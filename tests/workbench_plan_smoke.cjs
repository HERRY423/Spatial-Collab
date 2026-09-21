// Real service-backed production JS interactions through a DOM fixture.
const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const { createDOM } = require('./workbench_test_dom.cjs');
const [scriptPath, htmlPath, baseURL] = process.argv.slice(2);
const dom = createDOM(fs.readFileSync(htmlPath, 'utf8'));
const { get } = dom;
const requests = [];
const sandbox = { console, setTimeout, clearTimeout, setInterval, clearInterval,
  document: dom.document, window: { parent: { postMessage() {} }, addEventListener() {}, devicePixelRatio: 1 },
  ResizeObserver: class { observe() {} },
  fetch: async (url, options) => { requests.push({ name: url.split('/').pop(), args: JSON.parse(options.body) }); return fetch(baseURL + url, options); } };
let source = fs.readFileSync(scriptPath, 'utf8');
const startup = 'initialize().catch(error => status(error.message, "error"));';
assert(source.includes(startup));
source = source.replace(startup, 'globalThis.testHarness = { refresh, render, renderRun, readPlanSpec, getEditors: () => variantEditors, getPlan: () => plan, getBinding: () => planBinding, getState: () => state };');
vm.runInNewContext(source, sandbox, { filename: scriptPath });
const h = sandbox.testHarness;
const click = async id => { assert(get(id).listeners.click, `Missing click handler ${id}`); await get(id).listeners.click(); };
const input = async (node, value) => { node.value = value; await node.listeners[node.tag === 'select' ? 'change' : 'input']?.(); };
const plain = value => JSON.parse(JSON.stringify(value));
async function tool(name, args) {
  const response = await fetch(`${baseURL}/api/tool/${name}`, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(args) });
  const value = await response.json(); assert(value.ok, JSON.stringify(value)); return value.result;
}
function findButton(root, title) {
  const walk = node => node.tag === 'button' && node.textContent === title ? node : (node.children || []).map(walk).find(Boolean);
  const button = walk(root); assert(button, `Missing dynamic button ${title}`); return button;
}

(async () => {
  await h.refresh();
  const original = plain(h.getState()), originalSelection = original.selection.selection_id;
  await input(get('sourceLabel'), 'T'); await input(get('targetLabel'), 'B');
  await input(get('planRadii'), '2, 100'); await input(get('planScopes'), 'whole_slice');
  await input(get('planName'), 'All declared branches'); await input(get('planRationale'), 'Synthetic explicit UI plan');
  const first = h.getEditors()[0];
  await input(first.name, 'known target relabel'); await input(first.rationale, 'Known source measurement predicate');
  await input(first.label, 'C');
  await input(first.clauses[0].field, 'label'); await input(first.clauses[0].value, 'B');
  await findButton(first.row, '添加 AND 条件').listeners.click();
  await input(first.clauses[1].field, 'attributes.qc'); await input(first.clauses[1].op, 'ge');
  await input(first.clauses[1].valueType, 'number'); await input(first.clauses[1].value, '1');
  await input(get('featureQuery'), 'G'); await click('searchFeatures');
  assert.match(get('featureCatalog').textContent, /同符号对应多个 ID/);
  const featureTable = get('featureCatalog').children[0];
  const featureRow = featureTable.children.find(row => row.children?.[0]?.textContent === 'gene A,1');
  assert(featureRow);
  await findButton(featureRow, '加入标记检查').listeners.click();
  assert(JSON.parse(get('genes').value).includes('feature_id:gene A,1'));
  await findButton(featureRow, '加入末个假设的表达条件').listeners.click();
  assert.equal(first.clauses[2].field.value, 'counts.gene A,1');
  await click('inspect');
  assert.equal(requests.filter(r => r.name === 'inspect_selection').at(-1).args.genes.includes('feature_id:gene A,1'), true);
  await input(get('featureQuery'), ''); await click('searchFeatures');
  assert(!get('nextFeatures').disabled); await click('nextFeatures');
  assert.match(get('featureStatus').textContent, /51–54/);
  await click('previousFeatures'); assert.match(get('featureStatus').textContent, /1–50/);

  await click('addVariant'); const second = h.getEditors()[1];
  await input(second.name, 'explicit source IDs'); await input(second.rationale, 'Preserved exact IDs');
  await input(second.mode, 'ids'); await input(second.included, 'false');
  await findButton(second.row, '复制当前共享选区的确切 IDs').listeners.click();
  assert.deepEqual(JSON.parse(second.ids.value), ['s']);
  await click('addVariant'); const third = h.getEditors()[2];
  await input(third.name, 'missing measurement'); await input(third.rationale, 'Unknown measured field must remain unknown');
  await input(third.label, 'Unknown hypothesis');
  await input(third.clauses[0].field, 'attributes.unrecorded_qc'); await input(third.clauses[0].op, 'ge');
  await input(third.clauses[0].valueType, 'number'); await input(third.clauses[0].value, '1');
  const spec = plain(h.readPlanSpec());
  assert.equal(spec.variants[0].selector.predicate.all.length, 3);
  assert.equal(spec.selection_id, originalSelection);
  await click('savePlan'); assert.equal(h.getPlan()?.version, 1, get('status').textContent);
  assert(!get('freezePlan').disabled && get('runPlan').disabled);
  await click('freezePlan'); const frozen = plain(h.getPlan());
  assert.equal(frozen.version, 2, get('status').textContent); assert.equal(frozen.status, 'frozen');
  assert.match(get('planDetails').textContent, /未知 4/);
  assert.deepEqual(frozen.resolved.variants[0].cell_ids, ['t']);
  await click('runPlan');
  assert.match(get('comparisonStatus').textContent, /全部敏感性/);
  let run = JSON.parse(get('runDetails').textContent);
  assert.equal(run.results.length, 6); assert.equal(run.plan_sha256, frozen.plan_sha256);
  assert(run.status_counts.computed && run.status_counts.unknown && run.status_counts.failed);
  assert.match(get('resultDenominators').textContent, /失败（保留）/);
  assert.match(get('resultDenominators').textContent, /来源等权/);
  const table = get('resultDenominators').children.find(node => node.tag === 'table');
  assert.equal(table.children.length - 1, 6);
  assert.deepEqual(table.children.slice(1).map(row => row.children[0].textContent.split('\n')[0]), run.difference_order);
  assert(!get('resultDenominators').textContent.includes('undefined'));
  assert(!get('resultDenominators').textContent.includes('NaN'));

  // Another host's active ROI must not redirect the loaded frozen plan.
  await tool('set_selection', { expected_revision: original.head_revision, cell_ids: ['t'] });
  await h.refresh();
  assert.equal(h.getBinding().selection_id, originalSelection);
  assert.match(get('planBinding').textContent, /当前共享选区或视图不同/);
  assert.deepEqual(JSON.parse(h.getEditors()[1].ids.value), ['s']);
  await click('runPlan');
  assert.deepEqual(requests.filter(r => r.name === 'run_hypothesis_plan').at(-1).args, { plan_id: frozen.plan_id, version: 2, project_id: original.project_id });
  // A concurrently appended draft makes local edits stale. No overwrite or run.
  await input(get('planName'), 'unsaved local edit');
  assert(get('runPlan').disabled && get('freezePlan').disabled);
  const runsBefore = requests.filter(r => r.name === 'run_hypothesis_plan').length;
  await click('runPlan'); assert.equal(requests.filter(r => r.name === 'run_hypothesis_plan').length, runsBefore);
  const other = await tool('revise_hypothesis_plan', { plan_id: frozen.plan_id, expected_version: 2, spec: { ...frozen.spec, name: 'Other host draft' } });
  await click('savePlan'); assert.match(get('status').textContent, /Stale hypothesis plan/);
  assert.equal(h.getPlan().version, 2); assert.equal(get('planName').value, 'unsaved local edit');
  await click('loadPlan'); assert.match(get('status').textContent, /未保存修改/);
  await click('newPlan'); await click('refreshPlans');
  get('savedPlan').value = JSON.stringify({ plan_id: frozen.plan_id, version: other.version });
  await click('loadPlan'); assert.equal(h.getPlan().version, 3, get('status').textContent);
  await input(get('planName'), 'Reviewed next draft'); await click('savePlan');
  assert.equal(h.getPlan().version, 4); await click('freezePlan'); assert.equal(h.getPlan().version, 5);
  await click('runPlan'); assert.equal(requests.filter(r => r.name === 'run_hypothesis_plan').at(-1).args.version, 5);
  get('savedPlan').value = JSON.stringify({ plan_id: frozen.plan_id, version: 2 });
  await click('loadPlan'); assert.equal(h.getPlan().version, 2);
  assert.equal(h.getBinding().selection_id, originalSelection);
  assert.match(get('planStatus').textContent, /历史冻结版本仍可运行/);

  await click('refreshAssets'); assert.match(get('assetNotice').textContent, /图像 缺失/);
  await click('inspectAsset'); assert.match(get('assetSummary').textContent, /不同声明细胞/);
  assert.match(get('assetNotice').textContent, /仍不可判定/);
  assert.deepEqual(JSON.parse(get('assetDetails').textContent).observations.map(row => row.cell_id), ['t']);
  await click('export'); assert.equal(requests.filter(r => r.name === 'export_review_bundle').at(-1).args.compact, true);
  get('compactExport').checked = false; await click('export');
  assert.equal(requests.filter(r => r.name === 'export_review_bundle').at(-1).args.compact, false);
  await click('useCurrentPlanSelection'); assert.notEqual(h.getBinding().selection_id, originalSelection);
  assert(get('runPlan').disabled);
  const uncalibrated = plain(h.getState()); uncalibrated.metadata.units = 'array_index';
  h.render(uncalibrated); assert(get('compare').disabled);
  assert.match(get('coordinateNotice').textContent, /物理尺度未验证/);
  const after = await tool('open_project', {});
  assert.equal(after.head_revision, original.head_revision); assert.equal(after.revisions.length, 1);
  process.stdout.write(JSON.stringify({ explicit_feature_ids: true, feature_paging: true, and_predicates: true,
    exact_id_copy: true, save_freeze_exact_run: true, all_status_rows_kept: true, dual_methods: true,
    pinned_selection: true, stale_version_rejected: true, old_frozen_retrievable: true,
    asset_read_only_identity: true, compact_export: true, no_annotation_commit: true, coordinate_ui_gate: true }));
})().catch(error => { console.error(error); process.exitCode = 1; });
