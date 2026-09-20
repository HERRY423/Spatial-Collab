const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const { createDOM } = require('./workbench_test_dom.cjs');
const [scriptPath, htmlPath, baseURL] = process.argv.slice(2);
const dom = createDOM(fs.readFileSync(htmlPath, 'utf8')), { get } = dom;
const requests = [];
const sandbox = { console, setTimeout, clearTimeout, setInterval, clearInterval, document: dom.document,
  window: { parent: { postMessage() {} }, addEventListener() {}, devicePixelRatio: 1 },
  ResizeObserver: class { observe() {} },
  fetch: (url, options) => { requests.push({ name: url.split('/').pop(), args: JSON.parse(options.body) }); return fetch(baseURL + url, options); } };
let source = fs.readFileSync(scriptPath, 'utf8');
source = source.replace('initialize().catch(error => status(error.message, "error"));',
  'globalThis.h = { refresh, getState: () => state, getProtein: () => proteinViewData, getRun: () => proteinRun };');
vm.runInNewContext(source, sandbox, { filename: scriptPath });
const click = async id => { await get(id).listeners.click(); assert(!get('status').className?.includes('error'), get('status').textContent); };
(async () => {
  await sandbox.h.refresh();
  assert.equal(get('homeView').hidden, false);
  assert(get('homeMetrics').textContent.includes('RNA'));
  const selection = sandbox.h.getState().selection.selection_id;
  get('planRationale').value = 'unsaved human plan'; await get('planRationale').listeners.input();
  await click('navProtein'); assert.equal(get('rnaView').hidden, true); assert.equal(get('proteinView').hidden, false);
  await click('proteinInspect');
  assert.equal(sandbox.h.getProtein().summary.measured_count, 3);
  assert.equal(sandbox.h.getProtein().summary.zero_count, 1);
  assert(get('proteinStats').textContent.includes('缺失'));
  await click('navRna'); await click('navHome'); await click('openProtein');
  assert.equal(sandbox.h.getState().selection.selection_id, selection);
  assert.equal(get('planRationale').value, 'unsaved human plan');
  get('pairedRnaGene').value = 'R'; await click('proteinCompareAll');
  assert.equal(sandbox.h.getRun().parameters.features.length, 2);
  assert(get('proteinResults').textContent.includes('不可判定')); // constant panel-normalized RNA
  assert(get('proteinResults').textContent.includes('feature_id:P'));
  const stored = sandbox.h.getRun().run_id;
  get('proteinScale').value = 'log1p'; await get('proteinScale').listeners.change();
  assert.equal(sandbox.h.getProtein(), null); assert.equal(sandbox.h.getRun().run_id, stored);
  await click('proteinInspect'); assert.equal(sandbox.h.getProtein().scale, 'log1p');
  await click('homeQcAction'); assert.equal(get('rnaView').hidden, false);
  assert(!requests.some(r => ['apply_revision', 'propose_revision', 'set_selection'].includes(r.name)));
  console.log(JSON.stringify({ navigation_preserves_selection_and_draft: true, zero_and_unknown_retained: true,
    full_protein_comparison: true, stale_display_cleared: true, no_annotation_mutation: true }));
})().catch(error => { console.error(error); process.exitCode = 1; });
