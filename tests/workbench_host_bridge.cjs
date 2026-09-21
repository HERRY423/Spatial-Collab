// Real production UI + service calls behind a simulated MCP Apps host.
// This is NOT a ChatGPT/browser acceptance receipt.
const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const { webcrypto } = require('node:crypto');
const { createDOM } = require('./workbench_test_dom.cjs');
const [scriptPath, htmlPath, transferScript, baseURL, snapshotPath] = process.argv.slice(2);
const dom = createDOM(fs.readFileSync(htmlPath, 'utf8'));
dom.document.querySelector = () => ({content: ''});
dom.document.documentElement = {dataset: {}, scrollHeight: 900};
const messages = [], requests = [], callbacks = new Map(), intervals = new Map();
let intervalId = 0, listener, downloaded;
const parent = {postMessage(message) {
  messages.push(message);
  if (message.id === undefined) return;
  (async () => {
    let result;
    if (message.method === 'ui/initialize') result = {protocolVersion: '2026-01-26', hostCapabilities: {serverTools: {}, updateModelContext: {}}};
    else if (message.method === 'tools/call') {
      requests.push(message.params);
      const response = await fetch(`${baseURL}/${message.params.name}`, {method:'POST', headers:{'Content-Type':'application/json'}, body:JSON.stringify(message.params.arguments)});
      const payload = await response.json();
      result = payload.ok ? {structuredContent:payload.result} : {isError:true, structuredContent:{error:payload.error}};
    } else result = {};
    listener({source:parent, origin:'https://host.example', data:{jsonrpc:'2.0', id:message.id, result}});
  })().catch(error => { console.error(error); process.exitCode = 1; });
}};
const win = {parent, addEventListener(name, fn) {if(name==='message') listener=fn;}, devicePixelRatio:1};
const sandbox = {console, window:win, document:dom.document, crypto:webcrypto, Blob, Uint8Array, ArrayBuffer, btoa, atob,
  URL:{createObjectURL(blob){downloaded=blob;return 'blob:verified-test';},revokeObjectURL(){}},
  sessionStorage:{getItem:k=>callbacks.get(k),setItem:(k,v)=>callbacks.set(k,v),removeItem:k=>callbacks.delete(k)},
  ResizeObserver:class{observe(){}}, setTimeout, clearTimeout,
  setInterval(fn){intervals.set(++intervalId,fn);return intervalId;},clearInterval(id){intervals.delete(id);}};
const source = fs.readFileSync(scriptPath, 'utf8').replace('initialize().catch(error => status(error.message, "error"));', 'globalThis.ready = initialize();');
vm.createContext(sandbox); vm.runInContext(source, sandbox);
vm.runInContext(fs.readFileSync(transferScript, 'utf8'), sandbox);
const click = id => dom.get(id).listeners.click();
(async () => {
  await sandbox.ready;
  const api = win.SpatialCollab;
  const initial = api.currentProject();
  assert(initial);
  assert(messages.some(m=>m.method==='ui/notifications/initialized'));
  assert(messages.some(m=>m.method==='ui/update-model-context' && m.params.structuredContent.project_id===initial));
  assert.equal(intervals.size,1,'Embedded views must resync after model-side edits');
  listener({source:parent,origin:'https://attacker.example',data:{jsonrpc:'2.0',method:'ui/notifications/tool-result',params:{structuredContent:{project_id:'wrong',view:{}}}}});
  assert.equal(api.currentProject(),initial);
  await click('reloadProjects');
  const bytes = fs.readFileSync(snapshotPath);
  dom.get('uploadFile').files=[{name:'synthetic.json',size:bytes.length,arrayBuffer:async()=>bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength),slice:(start,end)=>new Blob([bytes.subarray(start,end)])}];
  await click('uploadAndImport');
  assert.match(dom.get('transferStatus').textContent,/新项目已创建/);
  assert.notEqual(api.currentProject(),initial);
  assert(requests.some(r=>r.name==='open_project'&&r.arguments.project_id===api.currentProject()));
  await click('downloadReview');
  assert.equal(dom.get('saveReviewFile').hidden,false,dom.get('transferStatus').textContent);
  assert(downloaded.size>0);
  assert.match(dom.get('transferStatus').textContent,/通过 SHA256/);
  const imported=api.currentProject();
  for(const poll of intervals.values()) await poll();
  assert.equal(api.currentProject(),imported);
  listener({source:parent,origin:'https://host.example',data:{jsonrpc:'2.0',id:'stop',method:'ui/resource-teardown'}});
  assert.equal(intervals.size,0);
  await assert.rejects(api.call('get_context'),/关闭/);
  fs.writeFileSync(snapshotPath+'.download.zip',Buffer.from(await downloaded.arrayBuffer()));
  console.log(JSON.stringify({handshake:true,origin_isolation:true,project_pinning:true,upload_import:true,download_checksum:true,embedded_resync:true,teardown:true,real_chatgpt:false}));
})().catch(error=>{console.error(error);process.exitCode=1;});
