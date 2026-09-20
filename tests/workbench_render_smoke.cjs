// DOM unit smoke: executes the real workbench event/render code with server fixtures.
// It does not claim browser layout, MCP Apps host or visible GUI acceptance.
const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const [scriptPath, htmlPath, payloadPath] = process.argv.slice(2);
const payload = JSON.parse(fs.readFileSync(payloadPath, 'utf8'));
const html = fs.readFileSync(htmlPath, 'utf8');
const { createDOM } = require('./workbench_test_dom.cjs');
const { get, document } = createDOM(html);
const requests = [];
const sandbox = {
  console, setTimeout, clearTimeout, setInterval, clearInterval,
  window: { parent: { postMessage() {} }, addEventListener() {}, devicePixelRatio: 1 },
  document,
  ResizeObserver: class { observe() {} },
  fetch: async (url, options) => {
    requests.push({ url, args: JSON.parse(options.body) });
    assert.equal(url, '/api/tool/inspect_quality');
    return { ok: true, json: async () => ({ ok: true, result: payload.quality }) };
  },
};
let source = fs.readFileSync(scriptPath, 'utf8');
const startup = 'initialize().catch(error => status(error.message, "error"));';
assert(source.includes(startup));
source = source.replace(startup, 'globalThis.testHarness = { render, renderRun, renderQuality };');
vm.runInNewContext(source, sandbox, { filename: scriptPath });

(async () => {
  sandbox.testHarness.render(payload.project);
  assert.equal(get('graphScope').querySelector('option[value="whole_slice"]').textContent, '选区来源 + 全导入窗口邻居');
  assert.equal(get('comparisonUniverse').querySelector('option[value="all"]').textContent, '全导入窗口');
  assert.match(get('scopeNotice').textContent, /未导入原始整张切片/);
  assert.match(get('sourceIdentity').textContent, /declared-slice-id/);
  assert.match(get('sourceIdentity').textContent, /空间窗口/);
  assert.match(get('sourceFiles').textContent, /fixture\/cells<untrusted>\.csv/);
  assert.match(get('annotationNotice').textContent, /自动标记组合工作假设/);
  assert.match(get('annotationNotice').textContent, /未经细胞身份验证/);
  const provenance = JSON.parse(get('sourceDetails').textContent);
  assert.deepEqual(provenance.source_files, payload.project.metadata.source_files);
  assert.deepEqual(provenance.working_annotation, payload.project.metadata.working_annotation);
  assert.deepEqual(provenance.import_scope, payload.project.metadata.import_scope);
  get('qualityScope').value = 'all';
  await get('inspectQuality').listeners.click();
  assert.equal(requests.length, 1);
  assert.equal(requests[0].args.selection_id, null);
  assert.equal(requests[0].args.revision_id, payload.project.view.revision_id);
  assert.equal(get('qualityRecord').hidden, false);
  assert.match(get('qualitySummary').textContent, /缺失及 null 不按零/);
  assert.match(get('qualitySummary').textContent, /细胞核面积/);
  assert(!get('qualitySummary').textContent.includes('NaN'));
  assert.deepEqual(JSON.parse(get('qualityDetails').textContent), payload.quality);
  sandbox.testHarness.renderRun(payload.sensitivity);
  assert.match(get('comparisonStatus').textContent, /计算假设/);
  assert.match(get('resultMetrics').textContent, /未创建修订/);
  assert.match(get('resultMeaning').textContent, /没有创建注释修订或记录研究者批准/);
  assert.match(get('resultMeaning').textContent, /全导入窗口/);
  assert.match(get('resultDenominators').textContent, /不可判定/);
  assert(!get('resultDenominators').textContent.includes('undefined'));
  assert.deepEqual(JSON.parse(get('runDetails').textContent), payload.sensitivity);
  sandbox.testHarness.renderRun(payload.region);
  assert.match(get('comparisonStatus').textContent, /区域描述性/);
  sandbox.testHarness.renderRun(payload.neighborhood);
  assert.match(get('comparisonStatus').textContent, /描述性判断/);
  const changed = structuredClone(payload.project);
  changed.head_revision = 'new-revision'; changed.view.revision_id = 'new-revision';
  sandbox.testHarness.render(changed);
  assert.equal(get('qualityRecord').hidden, true);
  assert.match(get('qualitySummary').textContent, /重新读取/);
  process.stdout.write(JSON.stringify({ crop_scope: true, qc_read_only_request: true, qc_missingness: true,
    hypothesis_renderer: true, existing_renderers: true, stale_qc_cleared: true, source_identity: true }));
})().catch(error => { console.error(error); process.exitCode = 1; });
