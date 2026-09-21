"use strict";
(() => {
  const $ = id => document.getElementById(id);
  const localToken = document.querySelector('meta[name="spatial-csrf"]').content;
  const local = Boolean(localToken);
  const palette = ["#b9aceb", "#79c6b5", "#e2ad7d", "#83b4df", "#d797b1", "#c8ce8a", "#90a6b0"];
  let state = null, proposal = null, run = null, quality = null, mode = "rect", dragging = null, polygon = [];
  let viewBounds = null, transform = null, zoomTimer = null, busy = false;
  let hostCapabilities = {}, requestId = 0, hostOrigin = null;
  let activeProjectId = null, tornDown = false;
  let contextCursor = null, pollTimer = null, polling = false;
  let plan = null, planDirty = false, variantEditors = [];
  let planBinding = null, planInitialized = false, planList = [], planContextHash = null;
  let assetList = [], assetInspection = null, featurePage = null, featureQuery = "";
  let page = "home", overview = null, proteinViewData = null, proteinRun = null, proteinTransform = null;
  let colorMode = "label";
  let activeAssetRaster = null, assetRasterImg = null, assetOpacity = 0.7;
  let cellExpressionMap = new Map(), activeFeatureName = "";
  let integrationJob = null;
  let lastIntegrationViews = null, lastIntegrationFrame = null, lastDisagreementSet = new Set(), syncedHoverCellId = null;
  const coordinateUnit = () => ({micrometer:"µm",pixel:"像素（物理尺度未验证）",array_index:"阵列坐标（物理尺度未验证）",unknown:"未知单位"}[state?.metadata.units] || "未知单位");
  function enableActions() {
    document.querySelectorAll("button,input,select,textarea").forEach(control => control.disabled = busy);
    if (state && state.metadata.units !== "micrometer") $("compare").disabled = true;
    $("freezePlan").disabled = busy || !plan || planDirty || plan.status !== "draft";
    $("runPlan").disabled = busy || !plan || planDirty || plan.status !== "frozen";
    $("previousFeatures").disabled = busy || !featurePage || featurePage.offset === 0;
    $("nextFeatures").disabled = busy || featurePage?.next_offset == null;
    $("inspectAsset").disabled = busy || !assetList.length;
    $("planBackgroundId").disabled = busy || $("planBackground").value !== "selection";
    for (const id of ["proteinInspect", "proteinCompare", "proteinCompareAll"]) $(id).disabled = busy || !overview?.protein_assays.length;
    $("proteinSelectVisible").disabled = busy || !proteinViewData?.points.length;
    $("proteinExport").disabled = busy || !proteinRun;
  }
  const pending = new Map();
  const canvas = $("map"), ctx = canvas.getContext("2d");
  const status = (message, kind = "") => { $("status").textContent = message; $("status").className = kind; };
  const short = id => id ? (id.length > 18 ? id.slice(0, 18) + "…" : id) : "—";
  const pretty = value => JSON.stringify(value, null, 2);
  const show = (id, visible) => { $(id).hidden = !visible; };
  const fmt = value => value == null ? "不可判定" : Number(value).toFixed(3);
  const selected = () => state?.selection;
  const selectedIds = () => new Set(selected()?.cell_ids || []);
  const textNode = (tag, text, className) => { const node = document.createElement(tag); node.textContent = text; if (className) node.className = className; return node; };
  function viridisColor(t) {
    const stops = [
      [68, 1, 84],
      [59, 82, 139],
      [33, 145, 140],
      [94, 201, 98],
      [253, 231, 37]
    ];
    t = Math.max(0, Math.min(1, t));
    const idx = t * (stops.length - 1);
    const i = Math.floor(idx);
    const f = idx - i;
    if (i >= stops.length - 1) return `rgb(${stops[stops.length - 1].join(",")})`;
    const c0 = stops[i], c1 = stops[i + 1];
    return `rgb(${Math.round(c0[0] + f * (c1[0] - c0[0]))},${Math.round(c0[1] + f * (c1[1] - c0[1]))},${Math.round(c0[2] + f * (c1[2] - c0[2]))})`;
  }

  function rpc(method, params, timeout = 30_000) {
    return new Promise((resolve, reject) => {
      if (tornDown) { reject(new Error("此视图已经关闭，请重新打开项目。")); return; }
      const id = `spatial-${++requestId}`;
      const timer = setTimeout(() => { pending.delete(id); reject(new Error("宿主未响应，请刷新或使用本地工作台。")); }, timeout);
      pending.set(id, { resolve, reject, timer });
      window.parent.postMessage({ jsonrpc: "2.0", id, method, params }, hostOrigin || "*");
    });
  }
  function notify(method, params) {
    window.parent.postMessage({ jsonrpc: "2.0", method, params }, hostOrigin || "*");
  }
  window.addEventListener("message", event => {
    if (local || event.source !== window.parent || (hostOrigin && event.origin !== hostOrigin)) return;
    const msg = event.data;
    if (!msg || msg.jsonrpc !== "2.0") return;
    if (pending.has(msg.id)) {
      // Bind to the actual parent's origin after handshake, with opaque iframe origins handled by '*'.
      if (!hostOrigin && event.origin !== "null") hostOrigin = event.origin;
      const item = pending.get(msg.id); clearTimeout(item.timer); pending.delete(msg.id);
      if (msg.error) item.reject(new Error(msg.error.message || "宿主拒绝请求")); else item.resolve(msg.result);
      return;
    }
    if (msg.method === "ui/notifications/tool-result") {
      const result = msg.params?.structuredContent;
      if (result?.project_id && result?.view) activeProjectId = result.project_id;
      if (result?.view && result?.head_revision) render(result);
      else if (state && !busy) refresh().catch(error => status(error.message, "error"));
    }
    if (msg.method === "ui/resource-teardown" && msg.id !== undefined) {
      tornDown = true;
      clearTimeout(zoomTimer);
      clearInterval(pollTimer);
      for (const item of pending.values()) { clearTimeout(item.timer); item.reject(new Error("视图已关闭。")); }
      pending.clear();
      window.parent.postMessage({ jsonrpc: "2.0", id: msg.id, result: {} }, hostOrigin || "*");
    }
    if (msg.method === "ping" && msg.id !== undefined) {
      window.parent.postMessage({ jsonrpc: "2.0", id: msg.id, result: {} }, hostOrigin || "*");
    }
    if (msg.method === "ui/notifications/host-context-changed") {
      const context = msg.params || {};
      if (context.theme) document.documentElement.dataset.theme = context.theme;
      if (context.availableDisplayModes) hostCapabilities.availableDisplayModes = context.availableDisplayModes;
    }
    if (msg.method === "ui/notifications/tool-input" && msg.params?.arguments?.project_id) {
      activeProjectId = msg.params.arguments.project_id;
    }
  });

  async function call(name, args = {}) {
    if (activeProjectId && !("project_id" in args)) args = { ...args, project_id: activeProjectId };
    const requestedProject = args.project_id;
    const assertCurrentProject = () => {
      if (requestedProject && activeProjectId !== requestedProject) throw new Error("项目已切换，旧请求的结果未显示。请在当前项目重试。");
    };
    if (local) {
      const response = await fetch(`/api/tool/${encodeURIComponent(name)}`, {
        method: "POST", headers: { "Content-Type": "application/json", "X-Spatial-CSRF": localToken }, body: JSON.stringify(args),
      });
      const payload = await response.json();
      assertCurrentProject();
      if (!response.ok || !payload.ok) throw new Error(payload.result?.error || payload.error || "请求失败");
      return payload.result;
    }
    if (!hostCapabilities.serverTools) throw new Error("此宿主未提供工具调用能力。请通过 Agent 工具或本地工作台继续。");
    const result = await rpc("tools/call", { name, arguments: args }, ["run_comparison", "run_hypothesis_plan", "run_protein_comparison", "export_review_bundle", "import_uploaded_project", "prepare_download"].includes(name) ? 120_000 : 30_000);
    assertCurrentProject();
    if (result.isError) throw new Error(result.structuredContent?.error || result.content?.find(x => x.type === "text")?.text || "工具执行失败");
    if (result.structuredContent) return result.structuredContent;
    const text = result.content?.find(x => x.type === "text")?.text;
    if (text) return JSON.parse(text);
    throw new Error("工具未返回结构化数据。");
  }

  async function updateContext() {
    if (local || !hostCapabilities.updateModelContext) return;
    const selection = selected();
    const context = {
      project_id: state.project_id, slice_id: state.metadata.slice_id,
      coordinate_system: state.metadata.coordinate_system, units: state.metadata.units,
      head_revision: state.head_revision,
      selection: selection ? { selection_id: selection.selection_id, revision_id: selection.revision_id,
        cell_ids_preview: selection.cell_ids.slice(0, 100), observation_count: selection.cell_ids.length,
        exact_ids_tool: "get_selection", polygon: selection.polygon, stale: selection.stale } : null,
      observation_semantics: state.semantics,
      proposal_id: proposal?.proposal_id || null, comparison_run_id: run?.run_id || null,
      hypothesis_plan: plan ? { plan_id: plan.plan_id, version: plan.version, status: plan.status,
        plan_sha256: plan.plan_sha256, unsaved_edits: planDirty, revision_id: plan.spec.revision_id,
        selection_id: plan.spec.selection_id } : null,
      scientific_authorization: "NOT_ESTABLISHED",
      active_workspace: page,
      atlas_view: page==="atlas"&&atlasCurrent()?{atlas_id:atlasCurrent().object_id,atlas_sha256:atlasCurrent().object_sha256,bounds:atlasBounds,layer:$("atlasLayer").value,feature:$("atlasGene").value||null,image_id:$("atlasImage").value||null,display_mode:atlasData?.mode}:null,
      protein_view: proteinViewData ? { assay_id: proteinViewData.assay_id, assay_sha256: proteinViewData.assay_sha256,
        feature: proteinViewData.feature, revision_id: proteinViewData.revision_id, selection_id: proteinViewData.selection_id } : null,
    };
    try { await rpc("ui/update-model-context", { structuredContent: context, content: [{ type: "text", text: `空间工作台共享对象：${pretty(context)}` }] }); }
    catch (error) { status(`状态已保存到服务端；宿主上下文未更新：${error.message}`, "error"); }
  }

  function action(id, fn) {
    $(id).addEventListener("click", async () => {
      if (busy) return;
      busy = true; enableActions();
      try { await fn(); }
      catch (error) { status(error.message, "error"); }
      finally { busy = false; enableActions(); }
    });
  }

  async function refresh() {
    const args = viewBounds ? { bounds: viewBounds } : {};
    const next = await call("open_project", args);
    const changed = state && state.head_revision !== next.head_revision;
    if (changed) { proposal = null; show("proposal", false); $("confirmation").checked = false; }
    render(next);
    const context = await call("get_context");
    contextCursor = context.head_revision === state.head_revision && context.selection_id === (selected()?.selection_id || null) ? context.cursor : null;
    if (context.hypotheses?.state_sha256 && context.hypotheses.state_sha256 !== planContextHash) {
      planContextHash = context.hypotheses.state_sha256;
      await refreshPlans();
    }
    if (context.latest_proposal_id && context.latest_proposal_id !== proposal?.proposal_id) {
      const sharedProposal = await call("get_proposal", { proposal_id: context.latest_proposal_id });
      if (!sharedProposal.stale && !sharedProposal.applied_revision && sharedProposal.selection_id === selected()?.selection_id) {
        proposal = sharedProposal; renderProposal(proposal); show("proposal", true); $("confirmation").checked = false;
      }
    }
    if (run) { run = await call("get_run", { run_id: run.run_id }); renderRun(run); }
    overview = await call("get_overview"); renderOverview();
    if (changed) status("共享版本已更新。旧选区需重新选择后才能修订；历史选区仍可用于前后比较。", "success");
    await updateContext();
  }

  function render(next) {
    activeProjectId = next.project_id || activeProjectId;
    const sourceChanged = state && (state.project_id !== next.project_id || state.source_sha256 !== next.source_sha256);
    const objectChanged = state && (state.head_revision !== next.head_revision || state.selection?.selection_id !== next.selection?.selection_id);
    if (objectChanged) {
      clearProteinView("共享选区或版本已变化，请重新读取蛋白证据。");
      for (const id of ["integrationLeftMap", "integrationRightMap"]) {
        const c = $(id); if (c) c.getContext("2d").clearRect(0, 0, c.width, c.height);
      }
      $("integrationReviewRows").replaceChildren();
      $("analysisRows").replaceChildren(); $("analysisResultSummary").textContent="共享状态已变化，请重新读取科学结果。";
      $("analysisMap").getContext("2d").clearRect(0,0,$("analysisMap").width,$("analysisMap").height);
      $("integrationSummary").textContent = "共享选区或版本已变化，请重新对照结果与当前对象。";
      $("integrationLeftTitle").textContent = "待重新读取";
      $("integrationRightTitle").textContent = "待重新读取";
      proposal = null; show("proposal", false); $("confirmation").checked = false;
      $("markerTable").replaceChildren(textNode("p", "共享版本或选区已变更。请在当前版本重新选择细胞并检查标记证据。", "muted"));
      assetInspection = null; show("assetRecord", false);
      $("assetSummary").replaceChildren(textNode("p", "共享版本或选区已改变，请重新核对对应对象。", "muted"));
    }
    if (quality && (quality.revision_id !== next.view.revision_id || (quality.selection_id && quality.selection_id !== next.selection?.selection_id))) {
      quality = null;
      $("qualitySummary").replaceChildren(textNode("p", "版本或 QC 选区已变更，请重新读取对应对象的测量记录。", "muted"));
      show("qualityRecord", false);
    }
    state = next;
    if (sourceChanged) {
      plan = null; planBinding = null; planInitialized = false; planDirty = false; variantEditors = [];
      planList = []; planContextHash = null; assetList = []; assetInspection = null; featurePage = null;
      $("assetChoice").replaceChildren(); $("featureCatalog").replaceChildren();
      $("planVariants").replaceChildren(); $("savedPlan").replaceChildren();
      $("genes").value = "";
    }
    const meta = state.metadata;
    const unit = state.semantics?.observation_unit || "cell";
    const unitName = { cell: "细胞", spot: "spot", bin: "bin" }[unit];
    $("projectName").textContent = meta.name || "空间转录组项目";
    $("projectMeta").textContent = `${state.semantics?.platform || meta.source_kind} · 切片 ${meta.slice_id} · ${meta.coordinate_system} · ${coordinateUnit()} · ${state.cell_count.toLocaleString()} 个${unitName}` + (state.view.revision_id !== state.head_revision ? ` · 历史视图；当前版本 ${short(state.head_revision)}` : "");
    $("coordinateCaption").textContent = `质心坐标 · ${coordinateUnit()}`;
    $("coordinateNotice").textContent = meta.units === "micrometer" ? "使用导入时声明的微米坐标；仍需核对原始标定依据。" : "物理尺度未验证：可浏览、选区、看表达和修订工作标签；物理半径分析不可用，假设计划将保留为未知。不会猜测阵列间距或像素大小。";
    const isWindow = meta.import_scope?.kind === "spatial_window";
    const importedUniverse = isWindow ? "全导入窗口" : "全部导入对象";
    renderSource(meta);
    $("scopeNotice").textContent = `当前对象：${unitName}。此处展示质心与原始标记计数；登记图像与掩膜的视觉复核在 napari 中进行，转录本未在此处复核。` + (unit === "cell" ? "细胞身份是待检验的工作注释；质心邻近不证明接触或通讯。" : "spot / bin 标签描述捕获区域，不等同于纯细胞身份或去卷积结果。") + (isWindow ? ` 当前为声明的导入窗口 ${JSON.stringify(meta.import_scope.bounds)} ${coordinateUnit()}，未导入原始整张切片。邻居与背景只来自此窗口。` : " 结果限于当前导入数据的描述性变化。");
    $("fitView").textContent = isWindow ? "适应导入窗口" : "适应导入数据";
    $("graphScope").querySelector('option[value="whole_slice"]').textContent = `选区来源 + ${importedUniverse}邻居`;
    $("comparisonUniverse").querySelector('option[value="all"]').textContent = importedUniverse;
    $("qualityScope").querySelector('option[value="all"]').textContent = importedUniverse;
    $("planScopes").querySelector('option[value="whole_slice"]').textContent = `${importedUniverse}邻居`;
    $("selectAll").textContent = `选择可见${unitName}`;
    $("editField").options[0].textContent = unit === "cell" ? "细胞工作标签" : "区域工作标签";
    $("headRevision").textContent = state.view.revision_id;
    $("sourceBadge").textContent = meta.source_kind === "synthetic" ? "合成演示 · 非真实组织" : `输入来源：${meta.source_kind}`;
    $("pointCount").textContent = `${state.view.returned_cells.toLocaleString()} / ${state.view.matching_cells.toLocaleString()} 个视野对象`;
    show("emptyView", !state.view.complete);
    const selection = selected();
    $("selectionCount").textContent = selection ? `已共享 ${selection.cell_ids.length.toLocaleString()} 个${unitName}` : "尚未选择对象";
    $("selectionState").textContent = selection ? (selection.stale ? "历史选区 · 修订前需重选" : `选区 ${short(selection.selection_id)}`) : "与 Agent 共享确切对象";
    $("selectionDetails").textContent = selection ? pretty(selection) : "尚无共享选区";
    if (!$("genes").value && meta.panel_genes?.length) $("genes").value = JSON.stringify(meta.panel_genes.slice(0, 5).map(id => `feature_id:${id}`));
    const revisions = state.revisions || [];
    $("revisionCount").textContent = `${revisions.length} 个版本`;
    $("revisionHistory").replaceChildren();
    for (const revision of [...revisions].reverse()) {
      const item = textNode("div", "", "history-item");
      item.append(textNode("strong", short(revision.revision_id) + (revision.revision_id === state.head_revision ? " · 当前" : "")));
      item.append(textNode("p", revision.rationale || revision.kind || "输入版本"));
      if (revision.reviewer) item.append(textNode("span", `记录确认者：${revision.reviewer}`));
      $("revisionHistory").append(item);
    }
    for (const id of ["baseRevision", "targetRevision", "revertTarget"]) {
      const old = $(id).value;
      $(id).replaceChildren();
      revisions.forEach((revision, i) => {
        const option = textNode("option", `${i === 0 ? "原始" : `版本 ${i}`} · ${short(revision.revision_id)}`);
        option.value = revision.revision_id; $(id).append(option);
      });
      if (id === "targetRevision") $(id).value = state.head_revision;
      else if (revisions.some(r => r.revision_id === old)) $(id).value = old;
    }
    const previousRun = $("retainedRun").value;
    $("retainedRun").replaceChildren();
    if (!state.runs?.length) { const empty = textNode("option", "尚无保存记录"); empty.value = ""; $("retainedRun").append(empty); }
    for (const saved of [...(state.runs || [])].reverse()) {
      const option = textNode("option", `${short(saved.run_id)} · ${saved.stale ? "历史结果" : "当前版本"}`);
      option.value = saved.run_id; $("retainedRun").append(option);
    }
    if ((state.runs || []).some(saved => saved.run_id === previousRun)) $("retainedRun").value = previousRun;
    if (run) { run.stale = run.target_revision !== state.head_revision; renderRun(run); }
    const labels = [...new Set(state.view.cells.map(c => c.label))].sort();
    $("legend").replaceChildren();
    for (const [i, label] of labels.entries()) {
      const item = textNode("span", "", "legend-item");
      const dot = textNode("span", "", "swatch"); dot.style.backgroundColor = palette[i % palette.length];
      item.append(dot, document.createTextNode(label)); $("legend").append(item);
    }
    $("legend").append(textNode("span", "○ 已排除细胞", "legend-item"));
    draw();
    if (!planInitialized) startNewPlan();
    renderPlanStatus();
    enableActions();
  }

  function renderSource(meta) {
    const scope = meta.import_scope;
    const working = meta.working_annotation;
    const fields = [
      ["来源切片 ID", meta.slice_id],
      ["样本 ID", meta.sample_id || "未单独记录"],
      ["记录平台", state.semantics?.platform || meta.platform || "未记录"],
      ["输入格式", meta.source_kind],
      ["导入特征数", String(meta.panel_gene_count ?? meta.panel_genes?.length ?? 0)],
      ["坐标系 / 单位", `${meta.coordinate_system} / ${meta.units}`],
      ["导入范围", scope?.kind === "spatial_window" ? `空间窗口 ${JSON.stringify(scope.bounds)} ${scope.units || "单位未记录"}` : (scope?.kind || "未声明空间裁剪；不能据此认定整张切片完整")],
    ];
    $("sourceIdentity").replaceChildren();
    for (const [label, value] of fields) {
      const item = textNode("div", "", "source-field");
      item.append(textNode("span", label, "muted"), textNode("strong", value || "未记录"));
      $("sourceIdentity").append(item);
    }
    $("sourceFiles").replaceChildren(textNode("h3", "原始来源文件"));
    const sources = Array.isArray(meta.source_files) ? meta.source_files : [];
    if (!sources.length) $("sourceFiles").append(textNode("p", "未记录原始来源文件。", "muted"));
    for (const source of sources) $("sourceFiles").append(textNode("code", source.path || "路径未记录"));
    const kind = working?.kind;
    $("annotationNotice").textContent = kind === "automated_exploratory_marker_rule"
      ? `自动标记组合工作假设 · 未经细胞身份验证。研究者批准记录：${working.researcher_approval || "未记录"}。当前标签来自声明的标记组合规则；可展开查看具体规则、来源版本与计划记录。`
      : kind ? `工作注释类型：${kind}。导入记录不构成细胞身份验证；请展开核对注释方法与审核状态。`
        : "未记录独立的工作注释方法；当前标签仍需结合原始来源和标记证据核对。";
    $("sourceDetails").textContent = pretty({ slice_id: meta.slice_id, sample_id: meta.sample_id || null, platform: state.semantics?.platform || meta.platform || null,
      panel_gene_count: meta.panel_gene_count ?? meta.panel_genes?.length ?? 0, feature_catalog: meta.feature_catalog || null,
      source_kind: meta.source_kind, source_files: sources, import_scope: scope || null, working_annotation: working || null });
  }

  function normalizedBounds() {
    const b = [...(viewBounds || state?.view.bounds || [0, 0, 100, 100])];
    if (b[0] === b[2]) { b[0] -= 1; b[2] += 1; }
    if (b[1] === b[3]) { b[1] -= 1; b[3] += 1; }
    return b;
  }

  function draw() {
    if (!state) return;
    const box = canvas.getBoundingClientRect(), dpr = window.devicePixelRatio || 1;
    if (!box.width || !box.height) return;
    canvas.width = Math.round(box.width * dpr); canvas.height = Math.round(box.height * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, box.width, box.height);
    const b = normalizedBounds(), pad = 42;
    const scale = Math.min((box.width - pad * 2) / (b[2] - b[0]), (box.height - pad * 2) / (b[3] - b[1]));
    const ox = (box.width - (b[2] - b[0]) * scale) / 2, oy = (box.height - (b[3] - b[1]) * scale) / 2;
    transform = { scale, ox, oy, bounds: b };
    ctx.strokeStyle = "#202831"; ctx.lineWidth = 1; ctx.fillStyle = "#637383"; ctx.font = "9px Consolas, monospace";
    for (let i = 0; i <= 4; i++) {
      const px = ox + i * (b[2] - b[0]) * scale / 4, py = oy + i * (b[3] - b[1]) * scale / 4;
      ctx.beginPath(); ctx.moveTo(px, oy); ctx.lineTo(px, box.height - oy); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(ox, py); ctx.lineTo(box.width - ox, py); ctx.stroke();
      ctx.fillText((b[0] + i * (b[2] - b[0]) / 4).toFixed(0), px - 9, oy - 12);
      ctx.fillText((b[1] + i * (b[3] - b[1]) / 4).toFixed(0), 7, py + 3);
    }
    if (activeAssetRaster && assetRasterImg) {
      ctx.save();
      ctx.globalAlpha = assetOpacity;
      const p = activeAssetRaster.pixel_to_world;
      const sc = transform.scale;
      const a = sc * p[0][0];
      const b = sc * p[1][0];
      const c = sc * p[0][1];
      const d = sc * p[1][1];
      const e = sc * p[0][2] + transform.ox - transform.bounds[0] * sc;
      const f = sc * p[1][2] + transform.oy - transform.bounds[1] * sc;
      ctx.transform(a, b, c, d, e, f);
      ctx.drawImage(assetRasterImg, 0, 0);
      ctx.restore();
    }
    if (colorMode === "density") {
      const cells = state.view.cells;
      const b = transform.bounds;
      const nx = 70, ny = 70;
      const dx = (b[2] - b[0]) / nx, dy = (b[3] - b[1]) / ny;
      if (dx > 0 && dy > 0) {
        const grid = new Float32Array(nx * ny);
        let maxDensity = 0;
        for (const cell of cells) {
          if (!cell.included) continue;
          const ix = Math.floor((cell.x - b[0]) / dx);
          const iy = Math.floor((cell.y - b[1]) / dy);
          if (ix >= 0 && ix < nx && iy >= 0 && iy < ny) {
            const count = ++grid[iy * nx + ix];
            if (count > maxDensity) maxDensity = count;
          }
        }
        if (maxDensity > 0) {
          for (let iy = 0; iy < ny; iy++) {
            for (let ix = 0; ix < nx; ix++) {
              const count = grid[iy * nx + ix];
              if (count === 0) continue;
              const norm = Math.pow(count / maxDensity, 0.65);
              const p0 = worldToPixel([b[0] + ix * dx, b[1] + iy * dy]);
              const p1 = worldToPixel([b[0] + (ix + 1) * dx, b[1] + (iy + 1) * dy]);
              ctx.fillStyle = viridisColor(norm);
              ctx.fillRect(Math.min(p0[0], p1[0]), Math.min(p0[1], p1[1]), Math.abs(p1[0] - p0[0]) + 0.5, Math.abs(p1[1] - p0[1]) + 0.5);
            }
          }
        }
      }
    } else if (colorMode === "expression") {
      const sel = selectedIds(), radius = state.view.cells.length > 3000 ? 2.2 : 3.2;
      let minExpr = 0, maxExpr = 0;
      if (cellExpressionMap.size > 0) {
        const vals = [...cellExpressionMap.values()];
        minExpr = Math.min(...vals); maxExpr = Math.max(...vals);
      }
      for (const cell of state.view.cells) {
        const p = worldToPixel([cell.x, cell.y]);
        let color;
        if (cellExpressionMap.has(cell.cell_id)) {
          const val = cellExpressionMap.get(cell.cell_id);
          const ratio = maxExpr > minExpr ? (val - minExpr) / (maxExpr - minExpr) : 0.5;
          color = viridisColor(ratio);
          ctx.globalAlpha = sel.size && !sel.has(cell.cell_id) ? 0.38 : 0.95;
        } else {
          color = "#4a5568";
          ctx.globalAlpha = 0.25;
        }
        ctx.beginPath(); ctx.arc(p[0], p[1], radius, 0, Math.PI * 2);
        if (cell.included) { ctx.fillStyle = color; ctx.fill(); } else { ctx.strokeStyle = color; ctx.stroke(); }
        if (sel.has(cell.cell_id)) { ctx.strokeStyle = "#f0e9ff"; ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(p[0], p[1], radius + 2, 0, Math.PI * 2); ctx.stroke(); }
        if (lastDisagreementSet.has(cell.cell_id)) { ctx.strokeStyle = "#ff4d6d"; ctx.lineWidth = 1.2; ctx.beginPath(); ctx.arc(p[0], p[1], radius + 3, 0, Math.PI * 2); ctx.stroke(); }
      }
    } else {
      const labels = [...new Set(state.view.cells.map(c => c.label))].sort();
      const sel = selectedIds(), radius = state.view.cells.length > 3000 ? 2.2 : 3.2;
      for (const cell of state.view.cells) {
        const p = worldToPixel([cell.x, cell.y]);
        const color = palette[labels.indexOf(cell.label) % palette.length];
        ctx.globalAlpha = sel.size && !sel.has(cell.cell_id) ? 0.38 : 0.88;
        ctx.beginPath(); ctx.arc(p[0], p[1], radius, 0, Math.PI * 2);
        if (cell.included) { ctx.fillStyle = color; ctx.fill(); } else { ctx.strokeStyle = color; ctx.stroke(); }
        if (sel.has(cell.cell_id)) { ctx.strokeStyle = "#f0e9ff"; ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(p[0], p[1], radius + 2, 0, Math.PI * 2); ctx.stroke(); }
        if (lastDisagreementSet.has(cell.cell_id)) { ctx.strokeStyle = "#ff4d6d"; ctx.lineWidth = 1.2; ctx.beginPath(); ctx.arc(p[0], p[1], radius + 3, 0, Math.PI * 2); ctx.stroke(); }
      }
    }
    ctx.globalAlpha = 1;
    if (dragging && dragging.current && mode === "rect") {
      const a = dragging.start, b2 = dragging.current;
      ctx.fillStyle = "#b6a6ec22"; ctx.strokeStyle = "#b6a6ec";
      ctx.fillRect(a[0], a[1], b2[0] - a[0], b2[1] - a[1]); ctx.strokeRect(a[0], a[1], b2[0] - a[0], b2[1] - a[1]);
    }
    if (polygon.length) {
      ctx.strokeStyle = "#c2b4f4"; ctx.fillStyle = "#c2b4f422"; ctx.beginPath();
      polygon.forEach((p, i) => { const q = worldToPixel(p); if (!i) ctx.moveTo(...q); else ctx.lineTo(...q); });
      if (polygon.length > 2) ctx.closePath(); ctx.stroke(); ctx.fill();
    }
  }
  function worldToPixel(p) { return [transform.ox + (p[0] - transform.bounds[0]) * transform.scale, transform.oy + (p[1] - transform.bounds[1]) * transform.scale]; }
  function pixelToWorld(p) { return [transform.bounds[0] + (p[0] - transform.ox) / transform.scale, transform.bounds[1] + (p[1] - transform.oy) / transform.scale]; }
  function eventPoint(event) { const r = canvas.getBoundingClientRect(); return [event.clientX - r.left, event.clientY - r.top]; }
  function nearest(p) {
    let best = null, distance = 100;
    for (const cell of state.view.cells) { const q = worldToPixel([cell.x, cell.y]), d = (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2; if (d < distance) { distance = d; best = cell; } }
    return best;
  }
  async function select(args) {
    if (!state) return;
    await call("set_selection", { expected_revision: state.view.revision_id, name: "研究者共享选区", ...args });
    proposal = null; show("proposal", false); $("confirmation").checked = false;
    $("markerTable").replaceChildren(textNode("p", "选区已更新，请重新检查标记证据。", "muted"));
    await refresh(); status("选区已保存，研究者与 Agent 可以读取同一组确切细胞。", "success");
  }
  canvas.addEventListener("pointerdown", event => {
    if (!state || busy || !state.view.complete || !transform) return;
    const p = eventPoint(event);
    if (mode === "polygon") { polygon.push(pixelToWorld(p)); draw(); return; }
    dragging = { start: p, current: p }; canvas.setPointerCapture(event.pointerId);
  });
  canvas.addEventListener("pointermove", event => {
    if (!state || !transform) return;
    const p = eventPoint(event);
    if (dragging) { dragging.current = p; draw(); show("hoverCard", false); return; }
    const cell = nearest(p); show("hoverCard", Boolean(cell));
    if (cell) {
      let text = `${cell.cell_id}\n${cell.label}${cell.included ? "" : " · 已排除"}\nx ${cell.x.toFixed(2)} · y ${cell.y.toFixed(2)} ${coordinateUnit()}`;
      if (cellExpressionMap.has(cell.cell_id)) {
        text += `\n${activeFeatureName || "特征"}: ${cellExpressionMap.get(cell.cell_id)} 计数`;
      }
      if (lastDisagreementSet.has(cell.cell_id)) {
        text += "\n⚠️ 整合方法分歧位点";
      }
      $("hoverCard").textContent = text;
      $("hoverCard").style.left = `${Math.min(p[0] + 14, canvas.clientWidth - 225)}px`;
      $("hoverCard").style.top = `${Math.max(0, p[1] - 80)}px`;
    }
  });
  canvas.addEventListener("pointerleave", () => show("hoverCard", false));
  canvas.addEventListener("pointercancel", () => { dragging = null; draw(); });
  canvas.addEventListener("pointerup", async event => {
    if (!dragging || busy) return;
    const a = dragging.start, b = eventPoint(event); dragging = null; draw();
    try {
      if (Math.hypot(a[0] - b[0], a[1] - b[1]) < 5) { const cell = nearest(b); await select({ cell_ids: cell ? [cell.cell_id] : [] }); }
      else {
        const p = pixelToWorld(a), q = pixelToWorld(b);
        await select({ polygon: [[Math.min(p[0], q[0]), Math.min(p[1], q[1])], [Math.max(p[0], q[0]), Math.min(p[1], q[1])], [Math.max(p[0], q[0]), Math.max(p[1], q[1])], [Math.min(p[0], q[0]), Math.max(p[1], q[1])]] });
      }
    } catch (error) { status(error.message, "error"); }
  });
  canvas.addEventListener("wheel", event => {
    if (!state || busy || !transform) return;
    event.preventDefault(); const center = pixelToWorld(eventPoint(event)), b = normalizedBounds(), factor = event.deltaY > 0 ? 1.28 : 0.78;
    viewBounds = [center[0] + (b[0] - center[0]) * factor, center[1] + (b[1] - center[1]) * factor, center[0] + (b[2] - center[0]) * factor, center[1] + (b[3] - center[1]) * factor];
    draw(); clearTimeout(zoomTimer); zoomTimer = setTimeout(() => refresh().catch(error => status(error.message, "error")), 220);
  }, { passive: false });
  new ResizeObserver(draw).observe($("canvasWrap"));

  function renderMarkers(evidence) {
    const data = evidence.expression || evidence.expression_summaries || {};
    cellExpressionMap.clear();
    const measuredEntry = Object.entries(data).find(([_, v]) => v.measured);
    if (measuredEntry && evidence.cells?.length) {
      activeFeatureName = measuredEntry[0];
      const featId = (measuredEntry[1].feature_resolution || {}).feature_id || activeFeatureName;
      for (const c of evidence.cells) {
        if (c.counts && c.counts[featId] != null) {
          cellExpressionMap.set(c.cell_id, c.counts[featId]);
        }
      }
      const vals = [...cellExpressionMap.values()];
      if (vals.length) {
        const minVal = Math.min(...vals), maxVal = Math.max(...vals);
        $("colorbarTitle").textContent = `${activeFeatureName} 表达量`;
        $("colorbarMin").textContent = fmt(minVal);
        $("colorbarMid").textContent = fmt((minVal + maxVal) / 2);
        $("colorbarMax").textContent = fmt(maxVal);
      }
      if (colorMode === "expression") draw();
    }
    const wrapper = textNode("div", "", "table-scroll"), table = document.createElement("table"), head = document.createElement("tr");
    ["基因", "测量状态", "均值", "非零细胞", "零计数细胞"].forEach(x => head.append(textNode("th", x))); table.append(head);
    for (const [gene, value] of Object.entries(data)) {
      const row = document.createElement("tr");
      const available = value.measured && value.cell_count > 0;
      row.style.cursor = available ? "pointer" : "default";
      if (available) {
        row.title = "点击在切片上查看该特征空间表达分布";
        row.addEventListener("click", () => {
          activeFeatureName = gene;
          cellExpressionMap.clear();
          const featId = (value.feature_resolution || {}).feature_id || gene;
          if (evidence.cells) {
            for (const c of evidence.cells) {
              if (c.counts && c.counts[featId] != null) {
                cellExpressionMap.set(c.cell_id, c.counts[featId]);
              }
            }
          }
          const vals = [...cellExpressionMap.values()];
          if (vals.length) {
            const minVal = Math.min(...vals), maxVal = Math.max(...vals);
            $("colorbarTitle").textContent = `${activeFeatureName} 表达量`;
            $("colorbarMin").textContent = fmt(minVal);
            $("colorbarMid").textContent = fmt((minVal + maxVal) / 2);
            $("colorbarMax").textContent = fmt(maxVal);
          }
          setColorMode("expression");
        });
      }
      [gene, value.status === "ambiguous" ? `符号有歧义：${value.feature_resolution.candidate_feature_ids.join(", ")}` : !value.measured ? "未测量" : (value.cell_count > 0 ? "面板已测量" : "面板内 · 选区为空"), available ? fmt(value.mean) : "—", available ? value.nonzero_cells : "—", available ? value.zero_cells : "—"].forEach(x => row.append(textNode("td", String(x))));
      table.append(row);
    }
    wrapper.append(table); $("markerTable").replaceChildren(wrapper);
    $("markerTable").prepend(textNode("p", `证据版本 ${short(evidence.selection.revision_id)} · 选区 ${short(evidence.selection.selection_id)} · ${evidence.cell_count} 个细胞（包含已排除细胞）`, "small muted"));
    $("markerTable").append(textNode("p", `标签计数：${JSON.stringify(evidence.label_counts || {})}`, "small muted"));
  }

  function renderProposal(value) {
    const holder = $("proposalSummary"); holder.replaceChildren();
    holder.append(textNode("p", `选中 ${value.selected_cell_count} 个细胞；实际修改 ${value.changed_cell_count} 个。`));
    holder.append(textNode("p", `理由：${value.rationale}`, "small muted"));
    const table = document.createElement("table"), head = document.createElement("tr");
    ["细胞", "修改"].forEach(text => head.append(textNode("th", text))); table.append(head);
    const names = { label: "标签", region: "区域", included: "纳入" };
    const display = value => typeof value === "boolean" ? (value ? "纳入" : "排除") : String(value || "（空）");
    for (const delta of value.deltas.slice(0, 20)) {
      const row = document.createElement("tr"); row.append(textNode("td", delta.cell_id));
      row.append(textNode("td", Object.entries(delta.changes).map(([field, change]) => `${names[field]}：${display(change.before)} → ${display(change.after)}`).join("；")));
      table.append(row);
    }
    const scroll = textNode("div", "", "table-scroll"); scroll.append(table); holder.append(scroll);
    if (value.deltas.length > 20) holder.append(textNode("p", `此处展示前 20 个细胞；完整 ${value.deltas.length} 个修改见展开记录。`, "small muted"));
    holder.append(textNode("p", `${(value.affected_result_ids || []).length} 条已有结果将失效并保留。基线 ${short(value.base_revision)}。`, "small muted"));
    $("proposalDetails").textContent = pretty(value);
  }

  function renderRun(result) {
    if (result.analysis_schema === "spatial-collab.protein-region.v1") { show("comparisonResult", false); renderProteinRun(result); return; }
    show("comparisonResult", true);
    $("runBadge").textContent = result.stale ? "历史结果 · 相对当前版本已过期" : `已保存 ${short(result.run_id)}`;
    $("runBadge").className = result.stale ? "badge amber" : "badge";
    if (result.analysis_schema === "spatial-collab.hypothesis-plan-run.v2") { renderPlanRun(result); return; }
    $("resultMetrics").replaceChildren();
    if (result.analysis_schema === "spatial-collab.hypothesis-sensitivity.v1") {
      renderHypothesisRun(result);
      return;
    }
    const comparison = result.comparison, before = result.before, after = result.after;
    const names = { stable: "描述性判断保持稳定", changed: "描述性判断发生变化", indeterminate: "当前证据不足以比较判断" };
    const classes = { descriptive_supported: "达到预设描述性阈值", descriptive_not_supported: "未达到预设描述性阈值", inconclusive: "不可判定" };
    $("comparisonStatus").textContent = names[comparison.status] || comparison.status;
    if (result.analysis_schema === "spatial-collab.region-contrast.v1") {
      $("comparisonStatus").textContent = { unchanged: "区域描述性指标保持不变", descriptive_values_changed: "区域描述性指标发生变化", indeterminate: "区域比较的分母不足" }[comparison.status] || comparison.status;
      const table = document.createElement("table"), head = document.createElement("tr");
      ["度量 · 前景 − 背景", "修订前", "修订后"].forEach(x => head.append(textNode("th", x))); table.append(head);
      const rows = [
        ["前景纳入 / 选中对象", `${before.foreground.included_count} / ${before.foreground.selected_count}`, `${after.foreground.included_count} / ${after.foreground.selected_count}`],
        ["背景纳入 / 选中对象", `${before.background.included_count} / ${before.background.selected_count}`, `${after.background.included_count} / ${after.background.selected_count}`],
      ];
      const labels = [...new Set([...Object.keys(before.label_fraction_difference), ...Object.keys(after.label_fraction_difference)])].sort();
      for (const label of labels) rows.push([`${label} · 比例差`, fmt(before.label_fraction_difference[label] ?? (before.status === "descriptive" ? 0 : null)), fmt(after.label_fraction_difference[label] ?? (after.status === "descriptive" ? 0 : null))]);
      for (const gene of result.parameters.genes) {
        rows.push([`${gene} · 检出率差`, fmt(before.marker_contrast[gene].detection_fraction_difference), fmt(after.marker_contrast[gene].detection_fraction_difference)]);
        rows.push([`${gene} · 每万面板计数均值差`, fmt(before.marker_contrast[gene].mean_panel_counts_per_10000_difference), fmt(after.marker_contrast[gene].mean_panel_counts_per_10000_difference)]);
      }
      for (const values of rows) { const row = document.createElement("tr"); values.forEach(x => row.append(textNode("td", x))); table.append(row); }
      $("resultDenominators").replaceChildren(table);
      $("resultMeaning").textContent = `前景：${result.selection.foreground_name}；背景：${result.selection.background_name}。对象集合冻结且互不重叠。比例差和标记差为描述性指标，不能解释为显著富集或差异表达；未测量及分母为空显示不可判定。`;
      $("runDetails").textContent = pretty(result); return;
    }
    for (const [label, value, detail] of [["修订前 · 超背景比例", before.excess_over_abundance, classes[before.classification]], ["修订后 · 超背景比例", after.excess_over_abundance, classes[after.classification]], ["效应变化 Δ", comparison.effect_delta, `标签变化 ${comparison.label_changes} · 纳入状态变化 ${comparison.exclusion_changes}`]]) {
      const card = textNode("div", "", "metric"); card.append(textNode("span", label), textNode("strong", fmt(value)), textNode("small", detail)); $("resultMetrics").append(card);
    }
    const table = document.createElement("table"), head = document.createElement("tr");
    ["度量与分母", "修订前", "修订后"].forEach(text => head.append(textNode("th", text))); table.append(head);
    const rows = [
      ["选区纳入细胞", item => `${item.composition.selected_included} / ${item.composition.selected_total}`],
      ["来源细胞 / 图节点", item => `${item.graph.source_count} / ${item.graph.node_count}`],
      ["目标邻接边 / 全部来源邻接边", item => `${item.graph.target_edges} / ${item.graph.directed_source_edges}`],
      ["观测目标邻居比例", item => fmt(item.observed_fraction)],
      ["排除自身后的目标背景比例", item => fmt(item.null_fraction)],
      ["无邻居的来源细胞", item => String(item.graph.isolated_sources)],
    ];
    for (const [label, get] of rows) { const row = document.createElement("tr"); [label, get(before), get(after)].forEach(text => row.append(textNode("td", text))); table.append(row); }
    $("resultDenominators").replaceChildren(table);
    $("resultMeaning").textContent = `前后采用相同冻结对象与参数。来源→目标：${result.parameters.source_label} → ${result.parameters.target_label}；半径 ${result.parameters.radius_um} µm。判断仅描述本切片；缺少类型或邻接边时保留不可判定。`;
    $("runDetails").textContent = pretty(result);
  }

  function renderHypothesisRun(result) {
    $("comparisonStatus").textContent = "计算假设的敏感性 · 尚未写入注释";
    const params = result.parameters;
    const states = { stable: "描述性判断稳定", changed: "描述性判断改变", indeterminate: "不可判定" };
    const coverage = {
      complete_for_radius_within_declared_window: "按声明完整窗口，半径邻居完整",
      window_boundary_may_truncate_neighbors: "窗口边界可能截断邻居",
      ROI_induced_neighbors_intentionally_truncated: "选区内图，边界截断由定义决定",
      not_established: "邻居覆盖未确定",
    };
    const cards = [["真实锚定版本", short(result.base_revision), "所有假设从此版本独立开始"],
      ["假设 × 半径", `${params.variants.length} × ${params.radii_um.length}`, `${params.radii_um.join(", ")} µm · 全部声明值保留`],
      ["注释提交状态", result.revision_created === false ? "未创建修订" : "需要核查", "计算假设不构成研究者批准"]];
    for (const [label, value, detail] of cards) {
      const card = textNode("div", "", "metric");
      card.append(textNode("span", label), textNode("strong", value), textNode("small", detail)); $("resultMetrics").append(card);
    }
    const table = document.createElement("table"), head = document.createElement("tr");
    table.className = "hypothesis-table";
    ["独立假设", "半径 µm", "实际变化 · 选区 / 全导入", "超背景比例 · 基线 → 假设", "Δ", "邻居目标比例 · 基线 → 假设", "背景目标比例 · 基线 → 假设", "来源对象 · 基线 → 假设", "目标边 / 来源边 · 基线 → 假设", "描述性判断", "邻居覆盖"].forEach(value => head.append(textNode("th", value)));
    table.append(head);
    for (const row of result.results) {
      const edge = item => `${item.graph.target_edges} / ${item.graph.directed_source_edges}`;
      const values = [row.variant, row.radius_um, `${row.changed_in_analysis_roi} / ${row.actual_changed_count}`,
        `${fmt(row.before.excess_over_abundance)} → ${fmt(row.after.excess_over_abundance)}`, fmt(row.comparison.effect_delta),
        `${fmt(row.before.observed_fraction)} → ${fmt(row.after.observed_fraction)}`, `${fmt(row.before.null_fraction)} → ${fmt(row.after.null_fraction)}`,
        `${row.before.graph.source_count} → ${row.after.graph.source_count}`, `${edge(row.before)} → ${edge(row.after)}`,
        states[row.comparison.status] || row.comparison.status, coverage[row.neighbor_coverage] || row.neighbor_coverage];
      const tr = document.createElement("tr"); values.forEach(value => tr.append(textNode("td", String(value)))); table.append(tr);
    }
    const assumptions = textNode("div", "", "hypothesis-assumptions");
    for (const variant of params.variants) assumptions.append(textNode("p", `${variant.name}：假设对 ${variant.cell_ids.length} 个确切对象采用 ${pretty(variant.changes)}。理由：${variant.rationale}`, "small muted"));
    $("resultDenominators").replaceChildren(table, assumptions);
    const graph = params.graph_scope === "roi_induced" ? "选区内部邻域" : result.import_scope?.kind === "spatial_window" ? "全导入窗口邻域与背景" : "全导入数据邻域与背景";
    $("resultMeaning").textContent = `来源→目标：${params.source_label} → ${params.target_label}；${graph}；超背景阈值 ${params.min_effect}。各假设独立使用相同真实基线，不将前一个假设累积到后一个。保存的是计算结果，没有创建注释修订或记录研究者批准。邻居完整性仅针对声明窗口和指定半径，不证明原始整张切片无偏、细胞身份正确或生物学结论成立。`;
    if (result.revision_created !== false || result.researcher_approval !== "NOT_REQUESTED_HYPOTHESIS_ONLY") $("resultMeaning").textContent = "此记录的修订 / 确认状态与假设计算约定不一致，请核查完整记录。";
    $("runDetails").textContent = pretty(result);
  }

  function renderQuality(result) {
    const holder = $("qualitySummary"); holder.replaceChildren();
    holder.append(textNode("p", `版本 ${short(result.revision_id)} · ${result.selection_id ? `冻结选区 ${short(result.selection_id)}` : "全部导入对象"} · ${result.observation_count.toLocaleString()} 个对象，包含 ${result.excluded_count} 个已排除对象。`, "small muted"));
    const expression = result.expression;
    holder.append(textNode("p", `导入表达特征：${expression.imported_feature_count.toLocaleString()}；表达库计数为零：${expression.zero_library_observation_count == null ? "未知（未导入表达特征）" : expression.zero_library_observation_count}。这些是测量描述，不自动判定失败细胞。`, "small"));
    const labels = { transcript_counts: "厂家转录本数", total_counts: "厂家总计数", control_probe_counts: "对照探针计数", genomic_control_counts: "基因组对照计数", control_codeword_counts: "对照 codeword", unassigned_codeword_counts: "未分配 codeword", deprecated_codeword_counts: "弃用 codeword", cell_area: "细胞面积（源单位）", nucleus_area: "细胞核面积（源单位）", nucleus_count: "细胞核数量" };
    const table = document.createElement("table"), head = document.createElement("tr");
    ["记录字段", "有效 N", "缺失 / null", "零值 N", "中位数 [Q1, Q3]", "最小 / 最大", "无效 / 溢出 N"].forEach(value => head.append(textNode("th", value))); table.append(head);
    const metrics = [["导入表达库计数", expression.library_counts], ["导入特征检出数", expression.detected_features], ...Object.entries(result.recorded_numeric_fields).map(([key, value]) => [labels[key] || key, value])];
    for (const [name, metric] of metrics) {
      const values = [name, metric.valid_count, `${metric.absent_count} / ${metric.null_count}`,
        metric.zero_count == null ? "未知" : metric.zero_count,
        metric.valid_count ? `${fmt(metric.median)} [${fmt(metric.q25)}, ${fmt(metric.q75)}]` : "未记录 / 不可用",
        `${fmt(metric.min)} / ${fmt(metric.max)}`, metric.invalid_type_count + metric.invalid_range_count + metric.overflow_count];
      const row = document.createElement("tr"); values.forEach(value => row.append(textNode("td", String(value)))); table.append(row);
    }
    const scroll = textNode("div", "", "table-scroll"); scroll.append(table); holder.append(scroll);
    holder.append(textNode("p", "导入表达库只汇总本次导入的实测特征；厂家总计数可能包含其他特征类别。缺失及 null 不按零计算，已排除对象仍在检查分母中。", "small muted"));
    const methods = result.segmentation_method;
    holder.append(textNode("p", `分割方法（厂家记录，非正确性评分）：${methods.categories.map(item => `${item.value}：${item.count}`).join("；") || "未记录"}。缺失 ${methods.absent_count}，null ${methods.null_count}。${methods.omitted_category_count ? `另有 ${methods.omitted_category_count} 类、${methods.omitted_observation_count} 个对象未在摘要展开。` : ""}`, "small muted"));
    $("qualityDetails").textContent = pretty(result); show("qualityRecord", true);
  }

  function parseQueries(value) {
    const raw = value.trim();
    if (!raw) return [];
    if (raw.startsWith("[")) {
      const values = JSON.parse(raw);
      if (!Array.isArray(values) || values.some(x => typeof x !== "string" || !x.trim())) throw new Error("请使用非空字符串组成的 JSON 列表。");
      return values;
    }
    return raw.split(/[,\n]+/).map(x => x.trim()).filter(Boolean);
  }

  function makeControl(tag, value, label, onChange = markPlanDirty) {
    const control = document.createElement(tag);
    control.value = value; control.ariaLabel = label;
    control.addEventListener(tag === "select" ? "change" : "input", onChange);
    return control;
  }
  function makeSelect(options, value, label, onChange = markPlanDirty) {
    const select = makeControl("select", value, label, onChange);
    for (const [key, title] of options) { const option = textNode("option", title); option.value = key; select.append(option); }
    select.value = value;
    return select;
  }
  function smallButton(title, callback) {
    const button = textNode("button", title, "quiet"); button.type = "button";
    button.addEventListener("click", () => { if (busy) return; try { callback(); } catch (error) { status(error.message, "error"); } });
    return button;
  }
  function markPlanDirty() { planDirty = true; renderPlanStatus(); enableActions(); }

  function addClause(editor, clause = { field: "label", op: "eq", value: "" }, mark = true) {
    if (editor.clauses.length >= 8) throw new Error("每个假设最多保留 8 个 AND 条件。");
    const row = textNode("div", "", "clause-row");
    const field = makeControl("input", clause.field, "条件字段");
    field.placeholder = "label / included / attributes.nucleus_count / counts.特征ID";
    const op = makeSelect([["eq", "等于"], ["ne", "不等于"], ["gt", ">"], ["ge", "≥"], ["lt", "<"], ["le", "≤"], ["in", "属于列表"]], clause.op, "条件运算");
    const typeName = Array.isArray(clause.value) ? "list" : clause.value === null ? "null" : typeof clause.value === "number" ? "number" : typeof clause.value === "boolean" ? "boolean" : "text";
    const valueType = makeSelect([["text", "文本"], ["number", "数字"], ["boolean", "布尔"], ["null", "null（未知）"], ["list", "JSON 列表"]], typeName, "条件值类型");
    const value = makeControl("input", typeof clause.value === "string" ? clause.value : JSON.stringify(clause.value), "条件值");
    value.placeholder = "布尔填 true / false；列表填 [1,2] 或 [\"A\",\"B\"]";
    const item = { row, field, op, valueType, value };
    row.append(field, op, valueType, value, smallButton("移除条件", () => {
      editor.clauses = editor.clauses.filter(x => x !== item); row.remove(); markPlanDirty();
    }));
    editor.clauses.push(item); editor.clauseHolder.append(row);
    if (mark) markPlanDirty();
  }

  function addVariant(variant = null, mark = true) {
    if (variantEditors.length >= 8) throw new Error("每个计划最多保留 8 个独立假设。");
    const row = document.createElement("tr");
    const name = makeControl("input", variant?.name || "", "假设名称"); name.placeholder = "给这个独立假设命名";
    const rationale = makeControl("textarea", variant?.rationale || "", "假设理由"); rationale.rows = 3;
    const mode = makeSelect([["predicate", "满足全部 AND 条件"], ["ids", "明确列出的细胞 ID"]], variant?.selector?.cell_ids ? "ids" : "predicate", "假设对象选择方式");
    const ids = makeControl("textarea", JSON.stringify(variant?.selector?.cell_ids || []), "假设的确切细胞 ID"); ids.rows = 3;
    const clauseHolder = textNode("div", "", "clauses");
    const idHolder = textNode("div", "", "variant-ids");
    const label = makeControl("input", variant?.changes?.label || "", "假设修改后的标签"); label.placeholder = "留空保留原标签";
    const includedValue = variant?.changes?.included === true ? "true" : variant?.changes?.included === false ? "false" : "";
    const included = makeSelect([["", "保留纳入状态"], ["false", "假设排除"], ["true", "假设纳入"]], includedValue, "假设纳入状态");
    const editor = { row, name, rationale, mode, ids, clauseHolder, clauses: [], label, included };
    idHolder.append(ids, smallButton("复制当前共享选区的确切 IDs", () => {
      if (!selected()?.cell_ids.length) throw new Error("请先创建非空共享选区。");
      ids.value = JSON.stringify(selected().cell_ids); markPlanDirty();
      status(`已复制 ${selected().cell_ids.length} 个确切 ID；后续共享选区变化不会替换此列表。`);
    }));
    const selector = document.createElement("td"), changes = document.createElement("td");
    const conditionControls = textNode("div", "", "condition-controls");
    conditionControls.append(clauseHolder, smallButton("添加 AND 条件", () => addClause(editor)));
    selector.append(mode, conditionControls, idHolder);
    const switchMode = () => { conditionControls.hidden = mode.value !== "predicate"; idHolder.hidden = mode.value !== "ids"; };
    mode.addEventListener("change", switchMode); switchMode();
    changes.append(label, included);
    const nameCell = document.createElement("td"), reasonCell = document.createElement("td"), actions = document.createElement("td");
    nameCell.append(name); reasonCell.append(rationale);
    actions.append(smallButton("移除假设", () => { variantEditors = variantEditors.filter(x => x !== editor); row.remove(); markPlanDirty(); }));
    row.append(nameCell, selector, changes, reasonCell, actions); $("planVariants").append(row); variantEditors.push(editor);
    for (const clause of variant?.selector?.predicate?.all || [{ field: "label", op: "eq", value: "" }]) addClause(editor, clause, false);
    if (mark) markPlanDirty();
    return editor;
  }

  function adoptPlanSelection() {
    if (!selected()?.cell_ids.length) throw new Error("请先创建非空共享 ROI，再绑定分析对象。");
    if (selected().source_sha256 !== state.source_sha256) throw new Error("共享 ROI 与当前来源不一致，请刷新。");
    planBinding = { revision_id: selected().revision_id, selection_id: selected().selection_id,
      cell_count: selected().cell_ids.length, source_sha256: selected().source_sha256 };
    markPlanDirty();
  }

  function startNewPlan() {
    planInitialized = true; plan = null; planBinding = null; variantEditors = [];
    $("planVariants").replaceChildren(); $("planDetails").textContent = "";
    $("planName").value = "注释敏感性复核"; $("planRationale").value = "";
    $("planRadii").value = "15, 35, 75"; $("planScopes").value = "both";
    $("planBackground").value = "neighbor_universe"; $("planBackgroundId").value = "";
    addVariant(null, false);
    if (selected()?.cell_ids.length) planBinding = { revision_id: selected().revision_id, selection_id: selected().selection_id,
      cell_count: selected().cell_ids.length, source_sha256: selected().source_sha256 };
    planDirty = false; renderPlanStatus(); enableActions();
  }

  function renderPlanStatus() {
    const binding = planBinding;
    $("planBinding").textContent = binding ? `固定版本 ${binding.revision_id} · ROI ${binding.selection_id} · ${binding.cell_count ?? "待冻结核对"} 个对象。` +
      (selected()?.selection_id !== binding.selection_id || state?.view.revision_id !== binding.revision_id ? " 当前共享选区或视图不同；计划仍保留原绑定，点击上方绑定按钮才会改用当前 ROI。" : "") : "尚未绑定分析版本与共享 ROI。先共享选区，再点击绑定按钮。";
    const latest = plan && planList.find(item => item.plan_id === plan.plan_id);
    $("planStatus").textContent = plan ? `${plan.spec.name} · ${plan.plan_id} · 版本 ${plan.version} · ${plan.status === "frozen" ? "已冻结" : "草稿"}` +
      (planDirty ? " · 有未保存编辑，需保存并重新冻结后才能运行。" : plan.status === "frozen" ? " · 可按此精确版本运行；编辑会生成新草稿。" : " · 已保存，尚未冻结。") +
      (latest && latest.version > plan.version ? ` 服务端最新为版本 ${latest.version}；历史冻结版本仍可运行，继续修改前需载入最新版本。` : "") :
      (planDirty ? "有未保存的计划编辑。保存草稿会记录本次对象和全部参数，不提交注释。" : "尚未保存计划。填写假设、理由并核对绑定对象后保存。");
  }

  function clauseValue(editor) {
    const raw = editor.value.value;
    if (editor.valueType.value === "text") return raw;
    if (editor.valueType.value === "null") return null;
    if (editor.valueType.value === "number") {
      if (!raw.trim() || !Number.isFinite(Number(raw))) throw new Error("数字条件必须填写有限数值，空值不会当成零。");
      return Number(raw);
    }
    if (editor.valueType.value === "boolean") {
      if (!["true", "false"].includes(raw.trim())) throw new Error("布尔条件只能填写 true 或 false。");
      return raw.trim() === "true";
    }
    const values = JSON.parse(raw);
    if (!Array.isArray(values) || !values.length || values.length > 32 || values.some(x => x !== null && !["string", "boolean", "number"].includes(typeof x))) throw new Error("列表条件需要 1–32 个 JSON 标量。");
    return values;
  }

  function readPlanSpec() {
    if (!planBinding) throw new Error("请先绑定非空共享 ROI 与版本。");
    if (planBinding.source_sha256 !== state.source_sha256) throw new Error("计划绑定了其他来源，需重新建立计划。");
    const radii = $("planRadii").value.split(/[,\s]+/).filter(Boolean).map(Number);
    if (!radii.length || radii.length > 6 || radii.some(x => !Number.isFinite(x) || x <= 0) || new Set(radii).size !== radii.length) throw new Error("请声明 1–6 个互不重复的正数半径。");
    const scopes = $("planScopes").value === "both" ? ["whole_slice", "roi_induced"] : [$("planScopes").value];
    if (!variantEditors.length || variantEditors.length * scopes.length * radii.length > 96) throw new Error("请保留 1–8 个假设，总比较组合不能超过 96。");
    const variants = variantEditors.map(editor => {
      const changes = {};
      if (editor.label.value.trim()) changes.label = editor.label.value.trim();
      if (editor.included.value !== "") changes.included = editor.included.value === "true";
      if (!Object.keys(changes).length) throw new Error("每个假设需要明确修改标签或纳入状态。");
      let selector;
      if (editor.mode.value === "ids") {
        const ids = parseQueries(editor.ids.value);
        if (!ids.length || new Set(ids).size !== ids.length) throw new Error("假设 ID 必须非空且互不重复。");
        selector = { cell_ids: [...ids].sort() };
      } else {
        if (!editor.clauses.length) throw new Error("条件选择至少需要一个 AND 条件。");
        selector = { predicate: { all: editor.clauses.map(clause => {
          const field = clause.field.value.trim(), op = clause.op.value, value = clauseValue(clause);
          if (!["label", "included", "region"].includes(field) && !/^attributes\..+|^counts\..+/.test(field)) throw new Error("条件字段需要 label、included、region、attributes.字段 或 counts.稳定特征ID。");
          if (["gt", "ge", "lt", "le"].includes(op) && typeof value !== "number") throw new Error("大小比较必须使用数字类型。");
          if (op === "in" && !Array.isArray(value)) throw new Error("属于列表运算必须使用 JSON 列表类型。");
          if (op !== "in" && Array.isArray(value)) throw new Error("JSON 列表只能用于属于列表运算。");
          return { field, op, value };
        }) } };
      }
      return { name: editor.name.value.trim(), rationale: editor.rationale.value.trim(), selector, changes };
    });
    const background = { kind: $("planBackground").value };
    if (background.kind === "selection") background.selection_id = $("planBackgroundId").value.trim();
    const minEffect = Number($("minEffect").value);
    if (!$("minEffect").value.trim() || !Number.isFinite(minEffect) || minEffect <= 0 || minEffect > 1) throw new Error("超背景阈值必须是大于 0 且不超过 1 的有限数字。");
    return { name: $("planName").value.trim(), rationale: $("planRationale").value.trim(),
      revision_id: planBinding.revision_id, selection_id: planBinding.selection_id,
      source_label: $("sourceLabel").value.trim(), target_label: $("targetLabel").value.trim(),
      min_effect: minEffect, radii_um: radii, graph_scopes: scopes, background, variants };
  }

  function applyPlan(record, selection = null) {
    plan = record; const spec = record.spec;
    planBinding = { revision_id: spec.revision_id, selection_id: spec.selection_id, source_sha256: state.source_sha256,
      cell_count: record.resolved?.selection?.cell_ids.length ?? selection?.cell_ids.length };
    $("planName").value = spec.name; $("planRationale").value = spec.rationale;
    $("planRadii").value = spec.radii_um.join(", "); $("planScopes").value = spec.graph_scopes.length === 2 ? "both" : spec.graph_scopes[0];
    $("planBackground").value = spec.background.kind; $("planBackgroundId").value = spec.background.selection_id || "";
    $("sourceLabel").value = spec.source_label; $("targetLabel").value = spec.target_label; $("minEffect").value = String(spec.min_effect);
    variantEditors = []; $("planVariants").replaceChildren();
    spec.variants.forEach(variant => addVariant(variant, false));
    planDirty = false; planInitialized = true;
    $("planDetails").textContent = pretty(record);
    if (record.resolved) {
      $("planDetails").textContent = `冻结解析：${record.resolved.variants.map(v => `${v.name} 命中 ${v.matched_count}，未知 ${v.unknown_count}`).join("；")}。未知条件会保留为未知结果。\n\n` + pretty(record);
    }
    renderPlanStatus(); enableActions();
  }

  async function refreshPlans() {
    const result = await call("list_hypothesis_plans"); planList = result.plans;
    const previous = $("savedPlan").value;
    $("savedPlan").replaceChildren();
    for (const item of [...planList].reverse()) for (const version of [...item.versions].reverse()) {
      const option = textNode("option", `${item.name} · v${version.version} · ${version.status === "frozen" ? "已冻结" : "草稿"}${version.version === item.version ? "（最新）" : "（历史）"}`);
      option.value = JSON.stringify({ plan_id: item.plan_id, version: version.version }); $("savedPlan").append(option);
    }
    if (!planList.length) { const empty = textNode("option", "尚无保存计划"); empty.value = ""; $("savedPlan").append(empty); }
    if ([...$("savedPlan").options].some(option => option.value === previous)) $("savedPlan").value = previous;
    renderPlanStatus();
  }

  function renderPlanRun(result) {
    const frozen = result.parameters.frozen_plan, spec = frozen.spec;
    $("comparisonStatus").textContent = "冻结计划的全部敏感性结果 · 未写入注释";
    const counts = result.status_counts || {};
    const metrics = $("resultMetrics"); metrics.replaceChildren();
    for (const [label, value, note] of [["完整比较组合", result.results.length, `计划 v${frozen.version} · ${spec.variants.length} 个独立假设`],
      ["已计算 / 未知 / 失败", `${counts.computed || 0} / ${counts.unknown || 0} / ${counts.failed || 0}`, "全部保留，未选择获胜假设"],
      ["真实注释版本", "未创建修订", `固定基线 ${short(result.base_revision)}`]]) {
      const card = textNode("div", "", "metric"); card.append(textNode("span", label), textNode("strong", String(value)), textNode("small", note)); metrics.append(card);
    }
    const holder = $("resultDenominators"); holder.replaceChildren();
    holder.append(textNode("p", "全部组合按两种估计量中最大的 |变化量| 排序；未知与失败保留在表尾。同一行并列展示按边计权和来源对象等权，差异不代表方法准确率。", "small muted"));
    const rowsById = new Map(result.results.map(row => [row.row_id, row]));
    const order = [...new Set(result.difference_order || [])].filter(id => rowsById.has(id));
    result.results.forEach(row => { if (!order.includes(row.row_id)) order.push(row.row_id); });
    const table = document.createElement("table"); table.className = "plan-results-table";
    const head = document.createElement("tr");
    ["组合 / 假设", "半径 / 范围", "状态 / 原因", "修改与未知对象", "按边计权：基线 → 假设", "按边变化 / 判定", "来源等权：基线 → 假设", "来源等权变化 / 判定", "等权 − 按边：基线 → 假设", "前后分母"].forEach(title => head.append(textNode("th", title)));
    table.append(head);
    const methodText = method => method?.status === "computed" ? `${fmt(method.excess_over_abundance)}（观察 ${fmt(method.observed_fraction)} / 背景 ${fmt(method.null_fraction)}；分母 ${method.denominator}）` : `不可判定：${method?.reason || method?.status || "未提供"}`;
    const comparisonText = item => `${fmt(item?.effect_delta)} · ${{ stable: "阈值判断稳定", changed: "跨越声明阈值", indeterminate: "不可判定" }[item?.status] || "不可判定"}`;
    const graphText = value => value?.graph ? `来源 ${value.graph.source_count}；有邻居来源 ${value.graph.sources_with_neighbors}；孤立 ${value.graph.isolated_sources}；边 ${value.graph.directed_source_edges}；目标边 ${value.graph.target_edges}；背景 ${value.graph.background_count}` : `未计算：${value?.reason || "未知"}`;
    for (const id of order) {
      const item = rowsById.get(id), row = document.createElement("tr"); row.className = `plan-result-${item.status}`;
      const reasons = [item.before?.reason, item.after?.reason, item.before?.error?.message, item.after?.error?.message].filter(Boolean);
      const scope = item.graph_scope === "roi_induced" ? "选区内部" : result.import_scope?.kind === "spatial_window" ? "全导入窗口邻居" : "全部导入邻居";
      const columns = [`${item.row_id}\n${item.variant}`, `${item.radius_um} µm\n${scope}`,
        `${{ computed: "已计算", unknown: "未知 / 不可判定", failed: "失败（保留）" }[item.status] || item.status}\n${[...new Set(reasons)].join("\n")}`,
        `命中 ${item.matched_count}；实际改 ${item.actual_changed_count}；ROI 内改 ${item.changed_in_analysis_roi}；条件未知 ${item.selector_unknown_count}`,
        `${methodText(item.before?.methods?.edge_weighted)}\n→ ${methodText(item.after?.methods?.edge_weighted)}`,
        comparisonText(item.comparison?.edge_weighted),
        `${methodText(item.before?.methods?.source_equal_weighted)}\n→ ${methodText(item.after?.methods?.source_equal_weighted)}`,
        comparisonText(item.comparison?.source_equal_weighted),
        `${fmt(item.reference_comparison?.before?.source_equal_minus_edge_weighted)} → ${fmt(item.reference_comparison?.after?.source_equal_minus_edge_weighted)}`,
        `${graphText(item.before)}\n→ ${graphText(item.after)}`];
      columns.forEach(value => row.append(textNode("td", value))); table.append(row);
    }
    holder.append(table);
    $("resultMeaning").textContent = `每个假设均从版本 ${result.base_revision} 和冻结 ROI ${spec.selection_id} 独立出发。背景：${{ neighbor_universe: "各邻域候选范围", imported_universe: "全部导入对象", selection: "指定冻结选区" }[spec.background.kind]}。来源等权排除孤立来源并保留其数量；两种估计量均不是显著性检验或独立验证。没有创建注释修订或记录研究者批准。` +
      (result.revision_created !== false || result.researcher_approval !== "NOT_REQUESTED_HYPOTHESIS_ONLY" ? " 记录的修订/批准字段与假设类型不一致，请核查完整记录。" : "");
    $("runDetails").textContent = pretty(result);
  }

  async function searchFeatures(offset = 0) {
    const source = state.source_sha256;
    const result = await call("get_feature_catalog", { query: featureQuery, limit: 50, offset });
    if (state.source_sha256 !== source || result.source_sha256 !== source) throw new Error("特征目录来源已改变，请刷新后重新搜索。");
    featurePage = result;
    $("featureStatus").textContent = `匹配 ${result.matching_count} 个特征；显示 ${result.features.length ? offset + 1 : 0}–${offset + result.features.length}。同名符号各自保留，不合并计数。`;
    const table = document.createElement("table"), header = document.createElement("tr");
    ["稳定特征 ID", "原符号", "符号是否歧义", "选择确切特征"].forEach(x => header.append(textNode("th", x))); table.append(header);
    for (const feature of result.features) {
      const row = document.createElement("tr"), actions = document.createElement("td");
      actions.append(smallButton("加入标记检查", () => {
        const queries = parseQueries($("genes").value), query = `feature_id:${feature.feature_id}`;
        if (!queries.includes(query)) queries.push(query);
        $("genes").value = JSON.stringify(queries); status(`已加入确切特征 ${feature.feature_id}；点击查看证据读取表达。`);
      }), smallButton("加入末个假设的表达条件", () => {
        const editor = variantEditors[variantEditors.length - 1] || addVariant();
        if (editor.mode.value !== "predicate") throw new Error("末个假设使用确切 IDs。请先选择 AND 条件方式，避免替换已有对象列表。");
        addClause(editor, { field: `counts.${feature.feature_id}`, op: "ge", value: 1 });
        status(`已添加 ${feature.feature_id} 原始计数 ≥ 1 的可编辑条件；保存前请核对阈值与假设理由。`);
      }));
      [feature.feature_id, feature.symbol ?? "未记录符号", feature.symbol_ambiguous ? "同符号对应多个 ID" : "未发现同名"].forEach(x => row.append(textNode("td", x)));
      row.append(actions); table.append(row);
    }
    $("featureCatalog").replaceChildren(table); enableActions();
  }

  async function refreshAssets() {
    const source = state.source_sha256, result = await call("list_assets");
    if (state.source_sha256 !== source || result.assets.some(a => a.source_sha256 !== source)) throw new Error("资产来源与当前项目不一致，请核对登记。");
    assetList = result.assets; const previous = $("assetChoice").value; $("assetChoice").replaceChildren();
    for (const asset of assetList) { const option = textNode("option", `${asset.kind === "image" ? "图像" : "分割"} · ${asset.name} · 映射 ${asset.mapping_status}`); option.value = asset.asset_id; $("assetChoice").append(option); }
    if (assetList.some(a => a.asset_id === previous)) $("assetChoice").value = previous;
    const prevMapAsset = $("mapAssetSelect").value;
    $("mapAssetSelect").replaceChildren();
    const noneOpt = textNode("option", "无底图"); noneOpt.value = ""; $("mapAssetSelect").append(noneOpt);
    for (const asset of assetList) {
      const option = textNode("option", `${asset.kind === "image" ? "图像" : "分割"} · ${asset.name}`);
      option.value = asset.asset_id;
      $("mapAssetSelect").append(option);
    }
    if (assetList.some(a => a.asset_id === prevMapAsset)) $("mapAssetSelect").value = prevMapAsset;
    $("assetNotice").textContent = `已登记 ${assetList.length} 个资产；图像 ${result.image_available ? "有登记" : "缺失"}，分割 ${result.segmentation_available ? "有登记" : "缺失"}。当前仅核对登记记录，读取对象对应时才校验文件内容；缺少证据保持未知，不判断合并分割或真实共表达。`;
    $("assetDetails").textContent = pretty(result); show("assetRecord", true);
    $("assetSummary").replaceChildren();
    for (const asset of assetList) $("assetSummary").append(textNode("p", `${asset.name} · ${asset.slice_id} · ${asset.coordinate_system} / ${asset.units} · ${asset.source.path} · 已映射 ${asset.mapped_label_count} / ${asset.label_count} 标签`, "small muted"));
    enableActions();
  }

  for (const id of ["planName", "planRadii", "planScopes", "planBackground", "planBackgroundId", "planRationale", "sourceLabel", "targetLabel", "minEffect"]) {
    $(id).addEventListener(["planScopes", "planBackground"].includes(id) ? "change" : "input", markPlanDirty);
  }
  action("addVariant", () => addVariant());
  action("newPlan", () => { startNewPlan(); status("已清空编辑区并建立新的未保存计划。已有版本仍保留。"); });
  action("useCurrentPlanSelection", () => { adoptPlanSelection(); status("计划已改用当前共享 ROI 及其版本；需保存新草稿并重新冻结。"); });
  action("refreshPlans", async () => { await refreshPlans(); status("已更新计划清单；当前编辑内容和所载版本保持不变。"); });
  action("loadPlan", async () => {
    if (planDirty) throw new Error("有未保存修改。请先保存草稿，或点击“另建计划（清空编辑）”明确清空后再载入。");
    if (!$("savedPlan").value) throw new Error("尚无可载入的计划。");
    const requested = JSON.parse($("savedPlan").value), source = state.source_sha256;
    const record = await call("get_hypothesis_plan", requested);
    const result = await call("get_selection", { selection_id: record.spec.selection_id });
    if (state.source_sha256 !== source || record.plan_id !== requested.plan_id || record.version !== requested.version ||
      result.selection?.source_sha256 !== source || result.selection?.revision_id !== record.spec.revision_id) throw new Error("计划返回的确切版本、选区或来源不匹配，未载入编辑区。");
    applyPlan(record, result.selection); await updateContext();
    status(`已载入计划版本 ${record.version}；分析对象固定为其原 ROI，当前共享选区未改变。`, "success");
  });
  action("savePlan", async () => {
    if (plan && !planDirty) { status("计划没有未保存修改；不额外创建相同草稿版本。"); return; }
    const spec = readPlanSpec(), source = state.source_sha256, previous = plan;
    const result = previous ? await call("revise_hypothesis_plan", { plan_id: previous.plan_id, expected_version: previous.version, spec }) : await call("create_hypothesis_plan", { spec });
    if (state.source_sha256 !== source) throw new Error("保存期间项目来源已变化。计划已保存在原项目，请重新载入核对。");
    applyPlan(result); await refreshPlans(); await updateContext();
    status(`已保存草稿版本 ${result.version}；未创建注释修订。请检查完整参数后冻结。`, "success");
  });
  action("freezePlan", async () => {
    if (!plan || planDirty || plan.status !== "draft") throw new Error("请先保存当前草稿，再冻结这个确切版本。");
    const source = state.source_sha256, expected = plan;
    const result = await call("freeze_hypothesis_plan", { plan_id: expected.plan_id, expected_version: expected.version });
    if (state.source_sha256 !== source) throw new Error("冻结期间项目来源已变化，请重新载入计划。");
    applyPlan(result); await refreshPlans(); await updateContext();
    status(`已冻结计划版本 ${result.version} 的全部参数与对象。未知条件已保留；冻结不是研究者科学批准。`, "success");
  });
  action("runPlan", async () => {
    if (!plan || planDirty || plan.status !== "frozen") throw new Error("只可运行已载入且没有未保存编辑的冻结版本。");
    const frozen = plan, source = state.source_sha256;
    status(`正在计算冻结计划 v${frozen.version} 的全部组合；未知和失败也会保存…`);
    const result = await call("run_hypothesis_plan", { plan_id: frozen.plan_id, version: frozen.version });
    if (state.source_sha256 !== source || result.plan_sha256 !== frozen.plan_sha256 || result.parameters?.frozen_plan?.version !== frozen.version) throw new Error("返回结果没有对应此精确冻结版本，请从记录清单重新核对。");
    run = result; renderRun(run); await refresh(); $("retainedRun").value = run.run_id; await updateContext();
    status(`冻结计划全部 ${run.results.length} 个组合已保存；注释版本保持不变。`, "success");
  });
  action("searchFeatures", () => { featureQuery = $("featureQuery").value.trim(); return searchFeatures(0); });
  action("previousFeatures", () => searchFeatures(Math.max(0, (featurePage?.offset || 0) - 50)));
  action("nextFeatures", () => { if (featurePage?.next_offset == null) throw new Error("已到目录末页。"); return searchFeatures(featurePage.next_offset); });
  action("refreshAssets", refreshAssets);
  $("assetChoice").addEventListener("change", () => { assetInspection = null; show("assetRecord", false); $("assetSummary").replaceChildren(); });
  action("inspectAsset", async () => {
    const assetId = $("assetChoice").value, selection = selected(), revision = state.view.revision_id;
    if (!assetId || !assetList.some(a => a.asset_id === assetId)) throw new Error("请先读取清单并选择已登记资产。");
    if (!selection?.cell_ids.length || selection.revision_id !== revision) throw new Error("请在当前视图版本创建非空共享选区，再核对资产对应。");
    const ids = selection.cell_ids.slice(0, 100), result = await call("inspect_asset", { asset_id: assetId, cell_ids: ids });
    if (selected()?.selection_id !== selection.selection_id || state.view.revision_id !== revision || result.revision_id !== revision ||
      result.asset.source_sha256 !== state.source_sha256 || result.asset.asset_id !== assetId || JSON.stringify(result.observations.map(o => o.cell_id)) !== JSON.stringify(ids)) throw new Error("资产核对期间版本、选区或确切 ID 已改变；请重新读取。");
    assetInspection = result;
    const holder = $("assetSummary"); holder.replaceChildren();
    holder.append(textNode("p", `选区 ${selection.selection_id} 共 ${selection.cell_ids.length} 个对象；本次明确检查前 ${ids.length} 个，其余 ${selection.cell_ids.length - ids.length} 个未检查。像素对应不是分割准确率。`, "small muted"));
    const table = document.createElement("table"), head = document.createElement("tr");
    ["确切对象 ID", "最近像素 (x, y)", "像素值 / mask 标签", "mask 声明的细胞 ID", "对应状态"].forEach(x => head.append(textNode("th", x))); table.append(head);
    const names = { same_declared_cell: "同一声明细胞", different_declared_cell: "不同声明细胞（需复核）", mapping_unknown: "映射缺失，未知", background_at_centroid: "质心处为背景", outside_image: "位于图像范围外", image_only: "仅图像，没有 mask 身份" };
    for (const item of result.observations) {
      const row = document.createElement("tr");
      [item.cell_id, JSON.stringify(item.nearest_pixel_xy), item.label_id ?? item.pixel_value ?? "未读取", item.mapped_cell_id ?? "未提供", names[item.correspondence] || item.correspondence].forEach(value => row.append(textNode("td", String(value)))); table.append(row);
    }
    const scroll = textNode("div", "", "table-scroll"); scroll.append(table); holder.append(scroll);
    $("assetNotice").textContent = "来源文件与像素内容已校验；细胞对应只说明登记记录和最近像素。是否合并分割、doublet 或真实共表达仍不可判定；视觉复核请在 napari 加载同一已登记资产。";
    $("assetDetails").textContent = pretty(result); show("assetRecord", true); status("已只读核对确切对象的资产对应；没有修改注释或共享选区。", "success");
  });

  action("refresh", async () => { await refresh(); status("已从服务端刷新共享状态。", "success"); });
  function setColorMode(mode) {
    colorMode = mode;
    $("colorModeLabel").classList.toggle("active", mode === "label");
    $("colorModeExpression").classList.toggle("active", mode === "expression");
    $("colorModeDensity").classList.toggle("active", mode === "density");
    show("expressionColorbar", mode === "expression");
    show("legend", mode === "label");
    draw();
  }
  action("colorModeLabel", () => setColorMode("label"));
  action("colorModeExpression", async () => {
    setColorMode("expression");
    if (cellExpressionMap.size === 0) {
      try {
        const genes = parseQueries($("genes").value);
        if (genes.length && selected()?.cell_ids.length) {
          const evidence = await call("inspect_selection", { genes });
          renderMarkers(evidence);
        }
      } catch (e) {
        status("请先在切片上选择细胞并输入标记特征，以读取连续表达分布。", "error");
      }
    }
  });
  action("colorModeDensity", () => setColorMode("density"));
  $("mapAssetSelect").addEventListener("change", async () => {
    const assetId = $("mapAssetSelect").value;
    if (!assetId) {
      activeAssetRaster = null;
      assetRasterImg = null;
      draw();
      return;
    }
    try {
      status("正在读取并对齐底图栅格…");
      const raster = await call("get_asset_raster", { asset_id: assetId });
      activeAssetRaster = raster;
      const img = new Image();
      img.onload = () => {
        assetRasterImg = img;
        draw();
        status(`已对齐加载底图 ${raster.name}（${raster.width}×${raster.height}，步长 ${raster.step}）。`, "success");
      };
      img.src = raster.data_url;
    } catch (err) {
      status(`底图加载失败: ${err.message}`, "error");
    }
  });
  $("mapAssetOpacity").addEventListener("input", e => {
    assetOpacity = Number(e.target.value) / 100;
    draw();
  });
  action("modeRect", () => { mode = "rect"; polygon = []; $("modeRect").classList.add("active"); $("modePolygon").classList.remove("active"); show("finishPolygon", false); draw(); });
  action("modePolygon", () => { mode = "polygon"; polygon = []; $("modePolygon").classList.add("active"); $("modeRect").classList.remove("active"); show("finishPolygon", true); status("依次点击多边形顶点，然后点击「完成多边形」。按质心是否落入区域选择细胞。"); });
  action("finishPolygon", async () => { if (polygon.length < 3) throw new Error("多边形至少需要三个顶点。"); await select({ polygon }); polygon = []; draw(); });
  action("fitView", async () => { viewBounds = null; polygon = []; await refresh(); });
  action("selectAll", async () => { if (!state.view.complete) throw new Error("请先缩小视野，不能将空视图作为完整导入数据选择。"); await select({ cell_ids: state.view.cells.map(c => c.cell_id) }); });
  action("selectIds", () => select({ cell_ids: $("cellIds").value.split(/[,\s]+/).filter(Boolean) }));
  action("inspect", async () => {
    const expectedSelection = selected()?.selection_id, expectedRevision = state.head_revision;
    const genes = parseQueries($("genes").value), evidence = await call("inspect_selection", { genes });
    if (evidence.selection.selection_id !== expectedSelection || evidence.selection.revision_id !== expectedRevision || selected()?.selection_id !== expectedSelection || state.head_revision !== expectedRevision) {
      await refresh(); throw new Error("另一个界面改变了共享选区或版本。请检查当前选区后重新读取标记。");
    }
    renderMarkers(evidence); status("标记证据已按共享选区读取。", "success");
  });
  action("inspectQuality", async () => {
    const revision = state.view.revision_id;
    const selectionId = $("qualityScope").value === "selection" ? selected()?.selection_id : null;
    if ($("qualityScope").value === "selection" && !selectionId) throw new Error("请先共享选区，或选择全部导入对象。");
    const result = await call("inspect_quality", { revision_id: revision, selection_id: selectionId });
    if (state.view.revision_id !== revision || result.revision_id !== revision || result.selection_id !== selectionId || (selectionId && selected()?.selection_id !== selectionId)) {
      await refresh(); throw new Error("QC 检查期间版本或选区已改变，请在当前对象上重新检查。");
    }
    quality = result; renderQuality(result); status("已读取测量与分割记录；没有改变注释、纳入状态或共享选区。", "success");
  });
  $("editField").addEventListener("change", () => { const included = $("editField").value === "included"; show("includedValue", included); show("editValue", !included); });
  action("preview", async () => {
    proposal = null; $("confirmation").checked = false; show("proposal", false);
    if (!selected()) throw new Error("请先创建共享选区。");
    const field = $("editField").value, value = field === "included" ? $("includedValue").value === "true" : $("editValue").value.trim();
    proposal = await call("propose_revision", { expected_revision: state.head_revision, selection_id: selected().selection_id,
      changes: { [field]: value }, rationale: $("rationale").value.trim(), actor: "researcher-workbench" });
    renderProposal(proposal); $("confirmation").checked = false; show("proposal", true);
    await updateContext(); status("修订预览已生成。请审阅确切修改后填写确认者并确认。", "success");
  });
  action("commit", async () => {
    if (!proposal) throw new Error("请先预览本次修订。");
    const confirmed = $("confirmation").checked; $("confirmation").checked = false;
    const result = await call("apply_revision", { proposal_id: proposal.proposal_id, expected_revision: proposal.base_revision,
      reviewer: $("reviewer").value.trim(), confirmation: confirmed });
    proposal = null; show("proposal", false); await refresh(); status(`已创建版本 ${result.revision_id}。请运行前后检验；旧结果保留。`, "success");
  });
  action("revert", async () => { await call("revert_revision", { target_revision: $("revertTarget").value, expected_revision: state.head_revision,
    reviewer: $("revertReviewer").value.trim(), confirmation: $("revertConfirmation").checked }); $("revertConfirmation").checked = false; await refresh(); status("历史覆盖层已作为新版本恢复。", "success"); });
  action("compare", async () => {
    const useSelection = $("comparisonUniverse").value === "selection";
    if (useSelection && !selected()) throw new Error("请先创建共享选区，或明确选择全部导入对象。");
    if (state.metadata.units !== "micrometer") throw new Error("物理半径比较需要明确的微米坐标；可以冻结计划并保留未标定的未知结果。");
    status("正在以相同参数重算修订前后指标…");
    run = await call("run_comparison", { base_revision: $("baseRevision").value, target_revision: $("targetRevision").value,
      selection_id: useSelection ? selected().selection_id : null, radius_um: Number($("radius").value), graph_scope: $("graphScope").value,
      source_label: $("sourceLabel").value.trim(), target_label: $("targetLabel").value.trim(), min_effect: Number($("minEffect").value) });
    renderRun(run); await refresh(); $("retainedRun").value = run.run_id; await updateContext(); status("修订前后检验已保存，可展开查看分母、邻接图变化与局限。", "success");
  });
  action("loadRun", async () => { if (!$("retainedRun").value) throw new Error("尚无保存的检验记录。"); run = await call("get_run", { run_id: $("retainedRun").value }); renderRun(run); await updateContext(); status("已读取原始参数与版本对应的保存记录。", "success"); });
  action("regionCompare", async () => {
    if (!selected()) throw new Error("请先创建前景共享选区。");
    run = await call("run_region_comparison", { base_revision: $("baseRevision").value, target_revision: $("targetRevision").value,
      selection_id: selected().selection_id, background_selection_id: $("backgroundSelection").value.trim() || null,
      genes: parseQueries($("genes").value) });
    renderRun(run); await refresh(); $("retainedRun").value = run.run_id;
    status("区域与背景的前后比较已保存；展开记录可检查每个分母与面板范围。", "success");
  });
  action("export", async () => { const result = await call("export_review_bundle", { run_id: run?.run_id || null, compact: $("compactExport").checked });
    $("exportResult").textContent = pretty(result); show("exportResult", true); status("审阅包已写入本项目的 exports 目录。校验和证明文件完整性，不证明生物学结论。", "success"); });

  async function refreshIntegrations() {
    const data = await call("list_integrations");
    for (const id of ["integrationLeft", "integrationRight"]) {
      const old = $(id).value; $(id).replaceChildren();
      for (const r of [...data.results].reverse()) {
        const o = textNode("option", `${r.method} · ${r.mode} · ${r.observation_count} 对象 · ${r.historical_view ? "历史" : "当前"} · ${short(r.result_id)}`);
        o.value = r.result_id; $(id).append(o);
      }
      if (data.results.some(r => r.result_id === old)) $(id).value = old;
    }
  }
  function integrationJobView(r) {
    integrationJob = r.id;
    $("integrationJobStatus").textContent = `任务 ${r.id} · ${r.status}${r.cache_hit ? " · 已复用相同输入结果" : ""}${r.error ? " · " + r.error : ""}`;
  }
  async function integrationPoll() {
    if (!integrationJob) throw new Error("尚无任务");
    const r = await call("get_integration_job", {job_id: integrationJob});
    integrationJobView(r);
    if (r.status === "succeeded") await refreshIntegrations();
  }
  action("integrationSubmit", async () => {
    if (!currentAssay()) throw new Error("请先在蛋白页选择配对蛋白层。");
    const spec = {
      revision_id: state.head_revision,
      assay_id: currentAssay().assay_id,
      backend: $("integrationBackend").value,
      feature_selection: $("integrationFeatureSelection").value,
      clusters: Number($("integrationClusters").value),
      components: Math.min(10, currentAssay().features.length),
      seed: Number($("integrationSeed").value),
      protein_transform: $("integrationProteinScale").value,
    };
    if ($("integrationScope").value === "roi") {
      if (!selected() || selected().stale) throw new Error("请先建立当前版本选区");
      spec.observation_ids = selected().cell_ids;
    }
    integrationJobView(await call("submit_integration", {spec}));
  });
  action("integrationPoll", integrationPoll);
  action("integrationCancel", async () => { if (!integrationJob) throw new Error("尚无任务"); integrationJobView(await call("cancel_integration_job", {job_id: integrationJob})); });
  action("integrationRetry", async () => { if (!integrationJob) throw new Error("尚无任务"); integrationJobView(await call("retry_integration_job", {job_id: integrationJob})); });
  action("integrationRefresh", refreshIntegrations);
  action("integrationFixed", async () => { const r = await call("filter_integration", {result_id: $("integrationLeft").value, revision_id: state.head_revision}); await refreshIntegrations(); $("integrationRight").value = r.result_id; status("已保留固定模型筛选结果；没有重新拟合。", "success"); });
  action("integrationRefit", async () => {
    const old = await call("get_integration", {result_id: $("integrationLeft").value});
    if (!old.method.parameters.backend) throw new Error("外部结果没有可运行后端；请在原环境重新拟合后导入。");
    const spec = {
      revision_id: state.head_revision,
      assay_id: old.input.protein_assay_id,
      rna_features: old.input.rna_features,
      protein_features: old.input.protein_features,
      components: old.method.parameters.components,
      clusters: old.method.parameters.clusters,
      backend: old.method.parameters.backend,
      feature_selection: old.method.parameters.feature_selection || "variance",
      seed: old.method.seed,
      protein_transform: old.preprocessing.protein,
    };
    if (old.fit_scope === "explicit_roi") { throw new Error("ROI 重新拟合请先在当前版本建立精确 ROI，再使用拟合按钮；不会静默改变拟合范围。"); }
    integrationJobView(await call("submit_integration", {spec}));
  });
  function renderIntegrationMaps() {
    if (!lastIntegrationViews || !lastIntegrationFrame) return;
    integrationMap("integrationLeftMap", lastIntegrationViews[0].points, lastIntegrationFrame, lastDisagreementSet, syncedHoverCellId);
    integrationMap("integrationRightMap", lastIntegrationViews[1].points, lastIntegrationFrame, lastDisagreementSet, syncedHoverCellId);
  }
  function integrationMap(id, points, frame, disagreements = new Set(), hoveredCellId = null) {
    const canvas = $(id), ratio = devicePixelRatio || 1; const w = canvas.clientWidth || 450, h = 340; canvas.width = w*ratio; canvas.height=h*ratio;
    const c = canvas.getContext("2d"); c.scale(ratio,ratio); c.clearRect(0,0,w,h);
    const colors = new Map([...new Set(points.map(p=>p.domain))].sort().map((d,i)=>[d,palette[i%palette.length]]));
    const scale = Math.min((w-24)/(frame[2]-frame[0]||1),(h-24)/(frame[3]-frame[1]||1));
    for(const p of points) {
      const px = 12+(p.x-frame[0])*scale, py = 12+(p.y-frame[1])*scale;
      c.fillStyle = p.included ? colors.get(p.domain) : "#666";
      c.beginPath(); c.arc(px, py, 2.5, 0, Math.PI*2); c.fill();
      if (disagreements.has(p.cell_id)) {
        c.strokeStyle = "#ff4d6d"; c.lineWidth = 1.2;
        c.beginPath(); c.arc(px, py, 4.5, 0, Math.PI*2); c.stroke();
      }
      if (hoveredCellId && p.cell_id === hoveredCellId) {
        c.strokeStyle = "#ffffff"; c.lineWidth = 2.0;
        c.beginPath(); c.arc(px, py, 6.0, 0, Math.PI*2); c.stroke();
      }
    }
  }
  action("integrationCompare", async () => {
    const args = {left_id: $("integrationLeft").value, right_id: $("integrationRight").value, left_partition: $("integrationLeftPartition").value, right_partition: $("integrationRightPartition").value, limit:50};
    const result = await call("compare_integrations", args);
    const views = await Promise.all([call("inspect_integration", {result_id:args.left_id,partition:args.left_partition,limit:10000}),call("inspect_integration", {result_id:args.right_id,partition:args.right_partition,limit:10000})]);
    lastIntegrationViews = views;
    const points = views.flatMap(v=>v.points);
    lastIntegrationFrame = [Math.min(...points.map(p=>p.x)),Math.min(...points.map(p=>p.y)),Math.max(...points.map(p=>p.x)),Math.max(...points.map(p=>p.y))];
    lastDisagreementSet = new Set(result.rows.map(r => r.cell_id));
    renderIntegrationMaps();
    $("integrationLeftTitle").textContent=`${args.left_partition} · ${views[0].mode} · ${views[0].historical_view ? "历史" : "当前"}`;
    $("integrationRightTitle").textContent=`${args.right_partition} · ${views[1].mode} · ${views[1].historical_view ? "历史" : "当前"}`;
    $("integrationSummary").textContent=`共同对象 ${result.common_count} · 左侧独有 ${result.left_only_count} · 右侧独有 ${result.right_only_count} · ARI ${fmt(result.adjusted_rand_index)} · 平均边界分歧 ${fmt(result.mean_boundary_disagreement)}。各图颜色独立编号；统计使用全部共同对象。${views.some(v=>v.next_offset!==null) ? "图仅展示前 10000 个对象，请通过 Agent 的视口分页检查其余对象。" : ""}`;
    $("integrationDetails").textContent=pretty({comparison:result,views:views.map(({points,...rest})=>rest)});
    const holder=$("integrationReviewRows"); holder.replaceChildren();
    for(const r of result.rows) {
      const row=textNode("div", "", "row disagreement-row");
      row.append(textNode("span", `${r.cell_id} · ${r.left_domain} / ${r.right_domain} · 边界分歧 ${fmt(r.boundary_disagreement)} · 当前标签 ${r.current_label}`));
      const spotBtn = textNode("button", "聚焦分歧点", "quiet");
      spotBtn.addEventListener("click", () => {
        syncedHoverCellId = r.cell_id;
        renderIntegrationMaps();
      });
      const button=textNode("button","检查原始证据 / 提出修订");
      button.addEventListener("click",async()=>{
        try {
          await select({cell_ids:[r.cell_id]});
          showPage("rna");
          renderMarkers(await call("inspect_selection", {genes:parseQueries($("genes").value)}));
        } catch(e){status(e.message,"error");}
      });
      row.append(spotBtn, button);
      holder.append(row);
    }
  });
  for (const [mid, idx] of [["integrationLeftMap", 0], ["integrationRightMap", 1]]) {
    $(mid).addEventListener("pointermove", e => {
      if (!lastIntegrationViews || !lastIntegrationFrame) return;
      const pts = lastIntegrationViews[idx].points;
      const rect = $(mid).getBoundingClientRect();
      const mx = e.clientX - rect.left, my = e.clientY - rect.top;
      const w = $(mid).clientWidth || 450, h = 340;
      const scale = Math.min((w-24)/(lastIntegrationFrame[2]-lastIntegrationFrame[0]||1),(h-24)/(lastIntegrationFrame[3]-lastIntegrationFrame[1]||1));
      let nearestCell = null, minDist = 120;
      for (const p of pts) {
        const px = 12+(p.x-lastIntegrationFrame[0])*scale, py = 12+(p.y-lastIntegrationFrame[1])*scale;
        const d = (px - mx)**2 + (py - my)**2;
        if (d < minDist) { minDist = d; nearestCell = p; }
      }
      if (nearestCell && syncedHoverCellId !== nearestCell.cell_id) {
        syncedHoverCellId = nearestCell.cell_id;
        renderIntegrationMaps();
      }
    });
    $(mid).addEventListener("pointerleave", () => {
      if (syncedHoverCellId) {
        syncedHoverCellId = null;
        renderIntegrationMaps();
      }
    });
  }
  action("integrationExport", async () => {const r=await call("export_review_bundle", {compact:true});$("integrationDetails").textContent=pretty(r);status("已导出源数据、修订和全部整合结果。");});
  action("multimodalContext",async()=>{$("multimodalDetails").textContent=pretty(await call("get_multimodal_context"));});

  let analysisJob = null, analysisPage = null, analysisCatalog = [], analysisPowerPlan = null;
  async function refreshAnalysis() {
    const [catalog, inputs, results, jobs] = await Promise.all([call("get_analysis_catalog"),call("list_analysis_inputs"),call("list_analysis_results"),call("list_analysis_jobs")]);
    function options(id, rows, value, label, blank=false) { const e=$(id), old=e.value; e.replaceChildren(); if(blank){const o=textNode("option","无第二输入");o.value="";e.append(o);} for(const row of rows){const o=textNode("option",label(row));o.value=value(row);e.append(o);} if([...e.options].some(o=>o.value===old))e.value=old; }
    analysisCatalog=catalog.methods;
    const methodSelect=$("analysisMethod"), oldMethod=methodSelect.value;methodSelect.replaceChildren();
    for(const [tier,label] of [["internal_numerical","内置数值方法"],["external_adapter","外部方法适配器"]]){const group=document.createElement("optgroup");group.label=label;for(const row of catalog.methods.filter(r=>r.method_contract.tier===tier)){const option=textNode("option",row.name+(row.available?"":" · 环境待安装"));option.value=row.method;group.append(option);}methodSelect.append(group);}if([...methodSelect.options].some(o=>o.value===oldMethod))methodSelect.value=oldMethod;
    options("analysisPrimary",inputs.inputs.filter(r=>r.kind==="spatial"),r=>r.object_id,r=>`${r.sample_id} · ${r.shape[0]} 对象 · ${r.object_id.slice(-8)}`);
    options("analysisSecondary",inputs.inputs,r=>r.object_id,r=>`${r.kind} · ${r.sample_id} · ${r.shape[0]} 对象`,true);
    options("analysisResult",results.results,r=>r.result_id,r=>`${r.method} · ${r.result_id.slice(-10)}`);
    options("analysisSavedJob",jobs.jobs,r=>r.id,r=>`${r.recipe.method} · ${r.status} · ${r.id.slice(-8)}`);
    updateAnalysisMethod();
  }
  function updateAnalysisMethod(){
    const method=$("analysisMethod").value, desc=analysisCatalog.find(x=>x.method===method);
    const fields={analysisSecondary:["nnls","cell2location","harmony","paste"],analysisFeatures:["mofa","mefisto","nnls","cell2location","harmony","paste","spagcn","spatial_spectral"],analysisComponents:["mofa","mefisto","harmony","spagcn","spatial_spectral"],analysisClusters:["spagcn","spatial_spectral"],analysisEpochs:["mofa","mefisto","cell2location","spagcn"],analysisPermutations:["moran_svg","spatial_lr"],analysisCellsPerSpot:["cell2location"]};
    for(const [id,methods] of Object.entries(fields))$(id).closest("label").hidden=!methods.includes(method);
    $("analysisDesign").hidden=!["harmony","paste"].includes(method);
    const hints={mofa:"比较 RNA 与蛋白的共同变化及各自贡献；本方法不使用空间坐标。",mefisto:"使用空间坐标和稀疏高斯过程拟合 RNA / 蛋白因子；需要同对象配对测量。",nnls:"将单细胞参考的表达签名拟合到空间位置，输出 RNA 贡献比例及残差。",cell2location:"根据单细胞参考和每位置细胞数先验估计丰度后验；请同时检查区间与训练记录。",harmony:"校正表达表示中的批次差异；请填写输入顺序对应的条件和批次。",paste:"估计完整重叠切片之间的空间对应；推断坐标与实测坐标分别保留。",spagcn:"使用官方 SpaGCN 的表达和空间邻接模型识别组织域。",spatial_spectral:"在稀疏空间图上结合表达相似性识别组织域，作为可复算基线。",moran_svg:"检验全部实测 RNA 基因的正空间自相关并统一校正；置换次数决定最小 p 值。",spatial_lr:"使用物种匹配的共识数据库，检验相邻位置的配体受体共表达；复合体要求全部亚基可测。",progeny:"使用 PROGENy 足迹与官方 ULM 推断 RNA 通路活性。当前使用学术资源范围；其他范围通过明确配方选择。"};
    $("analysisMethodHint").textContent=(hints[method]||"")+(desc&&!desc.available?" 当前分析环境尚缺该方法，请安装后重新读取方法列表。":"");
    $("analysisSubmit").disabled=!desc?.available;
    $("analysisMethodHint").textContent+=desc?` ${desc.method_contract.label}：${desc.method_contract.tier==="internal_numerical"?"核验指定数值定义与复算，不据此保证生物学正确性。":"负责输入冻结、调用与输出追踪；不为外部模型正确性背书。"}`:"";
    $("analysisPowerPanel").hidden=!["moran_svg","spatial_lr"].includes(method);
    $("powerModel").textContent=method==="moran_svg"?"生成模型：Gaussian SAR；rho 为邻接耦合参数，另报告所生成 Moran I 的分布。":"生成模型：lognormal neighbor coupling；rho 为邻域耦合参数，另报告相对随机摆放期望的 LR 超额分数。两者均未校准为真实组织的效应。";
  }
  function powerRecipe(){return {method:$("analysisMethod").value,input_ids:[$("analysisPrimary").value],parameters:{permutations:Number($("analysisPermutations").value)},seed:Number($("analysisSeed").value)};}
  function renderPowerCurve(r){const holder=$("analysisPowerCurve");holder.replaceChildren();if(!r.curve?.length)return;const table=document.createElement("table"),header=document.createElement("tr");for(const title of ["模型 rho","生成效应中位数","估计功效","功效 95% 区间"])header.append(textNode("th",title));table.append(header);for(const row of r.curve){const tr=document.createElement("tr");for(const value of [String(row.rho),row.generated_effect_median.toPrecision(4),`${(row.power*100).toFixed(1)}%`,row.power_ci95.map(v=>`${(v*100).toFixed(1)}%`).join("–")])tr.append(textNode("td",value));table.append(tr);}holder.append(table);}
  action("analysisPowerClear",async()=>{analysisPowerPlan=null;$("analysisPowerSummary").textContent="下次运行不关联声明；历史冻结记录仍保留。";$("analysisPowerCurve").replaceChildren();$("analysisPowerDetails").textContent="";});
  function powerText(r){const d=r.resolution;return `${d.observations} 个位置 / ${d.family_size} 项检验 / ${d.permutations} 次置换：最小 p=${d.minimum_p.toPrecision(3)}；BH 至少到第 ${d.minimum_possible_bh_rejection_rank} 位才有可能满足 q<${d.alpha}。假定排序第 ${d.assumed_bh_rank} 位至少需要 ${d.minimum_permutations_at_assumed_rank} 次置换。`+(r.status==="resolution_blocked_at_declared_rank"?" 此排序阈值受置换分辨率阻断，没有任何有限效应可以通过；其他信号使实际拒绝排序升高时需另行评估。":r.status==="simulated"?` 模型条件下的网格 MDE：${r.mde_rho===null?"当前网格尚未证明达到目标功效":`rho=${r.mde_rho}，生成效应中位数=${Number(r.mde_generated_effect_median).toPrecision(4)}`}。查看下方曲线与 95% 功效区间；这不是通用生物学最小效应。`:" 尚未指定模型并模拟最小可检出效应。")}
  action("analysisPowerCreate",async()=>{const spec=powerRecipe();const r=await call("create_analysis_power_plan",{spec,assumptions:{model:spec.method==="moran_svg"?"gaussian_sar":"lognormal_neighbor_coupling",bh_rank:Number($("powerRank").value),target_power:Number($("powerTarget").value),simulations:Number($("powerSimulations").value),effect_grid:$("powerGrid").value.split(",").map(Number),seed:spec.seed,rationale:$("powerRationale").value}});analysisPowerPlan=r;renderPowerCurve(r.result);$("analysisPowerSummary").textContent=powerText(r.result)+` 声明时点：${r.timing==="retrospective_sensitivity"?"本工作区已有相同输入的分析，属于回顾性敏感性评估":"在本工作区目标分析之前冻结；不等于外部预注册"}`;$("analysisPowerDetails").textContent=pretty(r);});
  $("analysisMethod").addEventListener("change",updateAnalysisMethod);
  action("analysisRefresh",refreshAnalysis);
  $("analysisSavedJob").addEventListener("change",async()=>{try{analysisJobStatus(await call("get_analysis_job",{job_id:$("analysisSavedJob").value}));}catch(e){status(e.message,"error");}});
  action("analysisSnapshot",async()=>{const args={revision_id:state.head_revision,species:$("analysisSpecies").value};if($("analysisScope").value==="roi"){const s=await call("get_selection");if(!s.selection||s.selection.stale||s.selection.revision_id!==state.head_revision)throw new Error("请在当前版本建立共享选区。");args.observation_ids=s.selection.cell_ids;}const r=await call("prepare_analysis_input",args);$("analysisStatus").textContent=`已冻结 ${r.shape[0]} 个对象、${r.shape[1]} 个 RNA 特征`;await refreshAnalysis();$("analysisPrimary").value=r.input_id;});
  function analysisJobStatus(r){analysisJob=r;$("analysisStatus").textContent=`${r.id} · ${r.status}${r.error?" · "+r.error:""}${r.cache_hit?" · 相同输入记录":""}`;}
  action("analysisSubmit",async()=>{const method=$("analysisMethod").value, input_ids=[$("analysisPrimary").value], p={};const number=id=>Number($(id).value);
    if(["nnls","cell2location","harmony","paste"].includes(method)){if(!$("analysisSecondary").value)throw new Error("请选择第二输入或单细胞参考。");input_ids.push($("analysisSecondary").value);}
    if(!["moran_svg","spatial_lr","progeny"].includes(method))p.n_features=number("analysisFeatures");
    if(["mofa","mefisto","harmony","spagcn","spatial_spectral"].includes(method))p.components=number("analysisComponents");
    if(["mofa","mefisto"].includes(method)){if(!currentAssay())throw new Error("请先在蛋白页选择配对蛋白层。");p.assay_id=currentAssay().assay_id;p.iterations=number("analysisEpochs");}
    if(["spagcn","spatial_spectral"].includes(method))p.clusters=number("analysisClusters");
    if(["spagcn","cell2location"].includes(method))p.max_epochs=number("analysisEpochs");
    if(method==="cell2location"){if(!$("analysisCellsPerSpot").value)throw new Error("请依据组织学提供每位置预期细胞数。");p.cells_per_location=number("analysisCellsPerSpot");}
    if(["moran_svg","spatial_lr"].includes(method))p.permutations=number("analysisPermutations");
    if(method==="paste"){if(!$("analysisFullOverlap").checked)throw new Error("请先检查并确认完整重叠假设。");p.overlap_assumption="full_overlap";}
    if(method==="harmony"){p.sample_conditions=$("analysisConditions").value.split(",").map(s=>s.trim());p.sample_batches=$("analysisBatches").value.split(",").map(s=>s.trim());}
    const spec={method,input_ids,parameters:p,seed:number("analysisSeed")};if(["moran_svg","spatial_lr"].includes(method)&&analysisPowerPlan)spec.power_plan_id=analysisPowerPlan.object_id;
    analysisJobStatus(await call("submit_analysis",{spec}));});
  action("analysisPoll",async()=>{if(!analysisJob&&$("analysisSavedJob").value)analysisJob={id:$("analysisSavedJob").value};if(!analysisJob)throw new Error("尚未提交分析任务。");analysisJobStatus(await call("get_analysis_job",{job_id:analysisJob.id}));if(analysisJob.status==="succeeded"){await refreshAnalysis();$("analysisResult").value=analysisJob.result_id;}});
  action("analysisCancel",async()=>{if(!analysisJob)throw new Error("尚无任务。");analysisJobStatus(await call("cancel_analysis_job",{job_id:analysisJob.id}));});
  action("analysisRetry",async()=>{if(!analysisJob)throw new Error("尚无任务。");analysisJobStatus(await call("retry_analysis_job",{job_id:analysisJob.id}));});
  async function inspectAnalysis(offset=0){const args={result_id:$("analysisResult").value,offset,limit:30};if($("analysisField").value)args.field=$("analysisField").value;const r=await call("inspect_analysis_result",args);analysisPage=r;$("analysisResultSummary").textContent=`${r.method} · ${r.total} 条结果 · 当前页 ${offset+1}–${offset+r.rows.length} · ${r.same_current_project_revision?"对应当前项目版本":"历史版本或外部样本，保留原始身份"}`;$("analysisDetails").textContent=pretty({...r,rows:undefined});const holder=$("analysisRows");holder.replaceChildren();for(const row of r.rows){const line=textNode("div","","row");line.append(textNode("span",row.cell_id?`${row.cell_id} · ${JSON.stringify(row.value)}`:JSON.stringify(row)));if(row.cell_id&&r.same_current_project_revision){const b=textNode("button","检查原始证据");b.addEventListener("click",async()=>{try{await select({cell_ids:[row.cell_id]});showPage("rna");renderMarkers(await call("inspect_selection",{genes:parseQueries($("genes").value)}));}catch(e){status(e.message,"error");}});line.append(b);}holder.append(line);}renderAnalysisMap(r);}
  function renderAnalysisMap(r,keepColumn=false){
    if(!keepColumn){$("analysisResultSummary").textContent+=` · ${r.method_contract?.label||"历史结果未分层"}`;const p=r.metadata.power_evidence;if(p){$("analysisResultSummary").textContent+=` · q<0.05：${r.metadata.discovery_summary?.discoveries??"未汇总"} 项。 `+powerText(p.plan?.result||{resolution:p.resolution});}}
    const rows=r.rows.filter(x=>x.xy), canvas=$("analysisMap"), ctx=canvas.getContext("2d");
    canvas.width=900;canvas.height=340;ctx.clearRect(0,0,900,340);
    if(!rows.length){$("analysisMapCaption").textContent="此结果按基因 / 分子关系汇总，不是逐位置地图。";return;}
    const width=Array.isArray(rows[0].value)?rows[0].value.length:1;
    if(!keepColumn){$("analysisValueColumn").replaceChildren();const names=r.metadata.cell_types||r.metadata.pathways||[];for(let i=0;i<width;i++){const o=document.createElement("option");o.value=i;o.textContent=names[i]||`分量 ${i+1}`;$("analysisValueColumn").append(o);}}
    const col=Number($("analysisValueColumn").value)||0, vals=rows.map(x=>Array.isArray(x.value)?x.value[col]:x.value);
    const numeric=vals.every(x=>typeof x==="number"&&Number.isFinite(x)), cats=[...new Set(vals)].sort();
    const xs=rows.map(x=>x.xy[0]),ys=rows.map(x=>x.xy[1]),xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);
    const low=numeric?Math.min(...vals):0,high=numeric?Math.max(...vals):Math.max(1,cats.length-1);
    const scale=Math.min(820/Math.max(xmax-xmin,1),270/Math.max(ymax-ymin,1));
    rows.forEach((x,i)=>{const v=numeric?vals[i]:cats.indexOf(vals[i]),f=(v-low)/Math.max(high-low,1e-12);ctx.fillStyle=`hsl(${240-220*f},70%,45%)`;ctx.beginPath();ctx.arc(40+(x.xy[0]-xmin)*scale,30+(x.xy[1]-ymin)*scale,5,0,2*Math.PI);ctx.fill();});
    $("analysisMapCaption").textContent=`仅本页 ${rows.length} / ${r.total} 个对象；实测输入坐标。${numeric?`蓝 ${low.toPrecision(3)} → 红 ${high.toPrecision(3)}`:`类别 ${cats.join(", ")}`}。完整范围请导出结果；原始证据通过下方对象按钮检查。`;
  }
  $("analysisValueColumn").addEventListener("change",()=>{if(analysisPage)renderAnalysisMap(analysisPage,true);});
  action("analysisInspect",()=>inspectAnalysis());
  action("analysisNext",()=>{if(!analysisPage||analysisPage.next_offset===null)throw new Error("没有下一页。");return inspectAnalysis(analysisPage.next_offset);});

  let atlases=[], atlasBounds=null, atlasData=null, atlasEdges=null, atlasImages=[], atlasTiles=[], atlasEpoch=0, atlasTimer=null, atlasDrag=null;
  const atlasCanvas=$("atlasMap"), atlasCtx=atlasCanvas.getContext("2d"), tileCache=new Map();
  const atlasCurrent=()=>atlases.find(r=>r.object_id===$("atlasChoice").value);
  function atlasTransform(){const b=atlasBounds,s=Math.min(atlasCanvas.width/(b[2]-b[0]),atlasCanvas.height/(b[3]-b[1]));return {s,x:(atlasCanvas.width-(b[2]-b[0])*s)/2,y:(atlasCanvas.height-(b[3]-b[1])*s)/2};}
  function atlasPoint(x,y){const b=atlasBounds,t=atlasTransform();return [t.x+(x-b[0])*t.s,t.y+(y-b[1])*t.s];}
  function atlasWorld(e){const r=atlasCanvas.getBoundingClientRect(),b=atlasBounds,t=atlasTransform();return [b[0]+((e.clientX-r.left)*atlasCanvas.width/r.width-t.x)/t.s,b[1]+((e.clientY-r.top)*atlasCanvas.height/r.height-t.y)/t.s];}
  function drawAtlas(){
    if(!atlasBounds)return;
    atlasCanvas.width=Math.max(400,atlasCanvas.clientWidth||900);atlasCanvas.height=560;
    const c=atlasCtx;c.clearRect(0,0,atlasCanvas.width,560);c.fillStyle="#10161b";c.fillRect(0,0,atlasCanvas.width,560);const frame=atlasTransform();c.save();c.beginPath();c.rect(frame.x,frame.y,(atlasBounds[2]-atlasBounds[0])*frame.s,(atlasBounds[3]-atlasBounds[1])*frame.s);c.clip();
    c.globalAlpha=Number($("atlasOpacity").value)/100;
    for(const tile of atlasTiles){const [p,q,r]=tile.world_corners.map(([x,y])=>atlasPoint(x,y));c.save();c.setTransform((q[0]-p[0])/tile.width,(q[1]-p[1])/tile.width,(r[0]-p[0])/tile.height,(r[1]-p[1])/tile.height,p[0],p[1]);c.drawImage(tile.image,0,0);c.restore();}
    c.globalAlpha=1;
    const records=atlasData?.records||[],maximum=Math.max(1,...records.map(r=>r.count||1));
    for(const r of records){const [x,y]=atlasPoint(r.x,r.y);if(atlasData.mode==="aggregate_density"){const [xx,yy]=atlasPoint(r.x+r.width,r.y+r.height);c.fillStyle=viridisColor(Math.log1p(r.count)/Math.log1p(maximum));c.globalAlpha=.8;c.fillRect(x,y,Math.max(1,xx-x),Math.max(1,yy-y));}else{c.fillStyle=atlasData.layer==="transcripts"?"#ffd37e":"#75c8bb";c.beginPath();c.arc(x,y,atlasData.layer==="transcripts"?1.8:2.5,0,Math.PI*2);c.fill();}}
    c.globalAlpha=1;
    for(const r of atlasEdges?.records||[]){c.strokeStyle=r.kind==="nucleus"?"#dbb0ff":"#86e9dc";c.beginPath();r.polygon.forEach(([x,y],i)=>{const p=atlasPoint(x,y);if(i)c.lineTo(...p);else c.moveTo(...p);});c.closePath();c.stroke();}c.restore();
  }
  async function readAtlasImages(){atlasImages=(await call("list_pyramids",{atlas_id:atlasCurrent().object_id})).images;const s=$("atlasImage");s.replaceChildren();const none=textNode("option","无图像");none.value="";s.append(none);for(const r of atlasImages){const o=textNode("option",`${r.levels.length} 层金字塔 · ${r.levels[0].width} × ${r.levels[0].height}`);o.value=r.object_id;s.append(o);}atlasTiles=[];tileCache.clear();}
  async function loadAtlasTiles(epoch){
    const img=atlasImages.find(r=>r.object_id===$("atlasImage").value);atlasTiles=[];if(!img)return;
    const a=img.pixel_to_world,b=atlasBounds,det=a[0][0]*a[1][1]-a[0][1]*a[1][0];
    const inv=(x,y)=>[(a[1][1]*(x-a[0][2])-a[0][1]*(y-a[1][2]))/det,(-a[1][0]*(x-a[0][2])+a[0][0]*(y-a[1][2]))/det];
    const corners=[[b[0],b[1]],[b[2],b[1]],[b[0],b[3]],[b[2],b[3]]].map(p=>inv(...p));
    const xs=corners.map(p=>p[0]),ys=corners.map(p=>p[1]),base=img.levels[0];let level=0;
    const desired=Math.max((Math.max(...xs)-Math.min(...xs))/atlasCanvas.width,(Math.max(...ys)-Math.min(...ys))/560);
    for(let i=1;i<img.levels.length;i++)if(base.width/img.levels[i].width<=desired)level=i;
    let tasks=[];
    function tasksAt(l){const v=img.levels[l],sx=v.width/base.width,sy=v.height/base.height,out=[];const x0=Math.max(0,Math.floor((Math.min(...xs)+.5)*sx/256)),x1=Math.min(Math.ceil(v.width/256)-1,Math.floor((Math.max(...xs)+.5)*sx/256)),y0=Math.max(0,Math.floor((Math.min(...ys)+.5)*sy/256)),y1=Math.min(Math.ceil(v.height/256)-1,Math.floor((Math.max(...ys)+.5)*sy/256));for(let y=y0;y<=y1;y++)for(let x=x0;x<=x1;x++)out.push({image_id:img.object_id,level:l,x,y});return out;}
    tasks=tasksAt(level);while(tasks.length>64&&level<img.levels.length-1)tasks=tasksAt(++level);
    if(tasks.length>64)throw new Error("该图像缺少足够粗的金字塔层；请放大视野，或生成更完整的图像金字塔。");
    let cursor=0;async function worker(){while(cursor<tasks.length&&epoch===atlasEpoch){const args=tasks[cursor++],key=JSON.stringify(args);let tile=tileCache.get(key);if(!tile){const r=await call("get_pyramid_tile",args),image=new Image();await new Promise((resolve,reject)=>{image.onload=resolve;image.onerror=()=>reject(new Error("图像瓦片无法解码"));image.src=r.data_url;});tile={...r,image};tileCache.set(key,tile);while(tileCache.size>64)tileCache.delete(tileCache.keys().next().value);}if(epoch===atlasEpoch){atlasTiles.push(tile);drawAtlas();}}}
    await Promise.all(Array.from({length:Math.min(4,tasks.length)},()=>worker()));
  }
  async function loadAtlas(){
    const current=atlasCurrent();if(!current)return;const epoch=++atlasEpoch;atlasBounds=atlasBounds||current.bounds.slice();atlasTiles=[];
    const args={atlas_id:current.object_id,bounds:atlasBounds,layer:$("atlasLayer").value};if(args.layer==="transcripts"&&$("atlasGene").value.trim())args.feature=$("atlasGene").value.trim();
    $("atlasStatus").textContent="正在按视野读取磁盘索引…";
    const [r,edges]=await Promise.all([call("get_atlas_view",args),$("atlasBoundaries").checked?call("get_atlas_view",{atlas_id:current.object_id,bounds:atlasBounds,layer:"boundaries"}):Promise.resolve(null)]);if(epoch!==atlasEpoch)return;
    atlasData=r;atlasEdges=edges;$("atlasStatus").textContent=`${current.metadata.name} · 当前视野 ${r.total.toLocaleString()} 个${r.layer==="transcripts"?"分子":"对象"} · ${r.mode==="aggregate_density"?"全部计入密度格，放大查看精确对象":"精确位置"}${edges?.mode==="zoom_required"?" · 分割边界过密，请继续放大":""}`;
    $("atlasDetails").textContent=pretty({dataset:current,view:{...r,records:undefined},boundaries:edges?{...edges,records:undefined}:null});drawAtlas();await call("share_atlas_view",{...args,image_id:$("atlasImage").value||null});if(epoch!==atlasEpoch)return;await updateContext();await loadAtlasTiles(epoch);
  }
  async function refreshAtlas(){atlases=(await call("list_atlases")).atlases;const s=$("atlasChoice"),old=s.value;s.replaceChildren();for(const r of atlases){const o=textNode("option",`${r.metadata.name} · ${r.totals.cells.toLocaleString()} 对象`);o.value=r.object_id;s.append(o);}if(atlases.some(r=>r.object_id===old))s.value=old;if(atlasCurrent()){atlasBounds=atlasCurrent().bounds.slice();await readAtlasImages();await loadAtlas();}else $("atlasStatus").textContent="尚未登记大型数据集。请按 large-data 操作说明导入；原有项目继续保留。";await refreshStudies();}
  const atlasAsync=fn=>()=>Promise.resolve().then(fn).catch(e=>{$("atlasStatus").textContent=e.message;});
  action("atlasRefresh",refreshAtlas);action("atlasApply",loadAtlas);action("atlasFit",async()=>{atlasBounds=atlasCurrent()?.bounds.slice();await loadAtlas();});
  $("atlasChoice").addEventListener("change",atlasAsync(async()=>{atlasBounds=atlasCurrent().bounds.slice();await readAtlasImages();await loadAtlas();}));
  for(const id of ["atlasLayer","atlasBoundaries","atlasImage"])$(id).addEventListener("change",atlasAsync(loadAtlas));$("atlasOpacity").addEventListener("input",drawAtlas);
  atlasCanvas.addEventListener("wheel",e=>{if(!atlasBounds)return;e.preventDefault();const p=atlasWorld(e),scale=e.deltaY>0?1.5:1/1.5;atlasBounds=atlasBounds.map((v,i)=>p[i%2]+(v-p[i%2])*scale);clearTimeout(atlasTimer);atlasTimer=setTimeout(atlasAsync(loadAtlas),180);},{passive:false});
  atlasCanvas.addEventListener("pointerdown",e=>{if(atlasBounds){atlasDrag={point:atlasWorld(e),bounds:atlasBounds.slice(),pan:e.shiftKey};atlasCanvas.setPointerCapture(e.pointerId);}});
  atlasCanvas.addEventListener("pointermove",e=>{if(!atlasBounds||atlasData?.mode!=="exact_points")return;const p=atlasWorld(e),t=atlasTransform();let nearest=null,distance=64;for(const row of atlasData.records){const d=((row.x-p[0])*t.s)**2+((row.y-p[1])*t.s)**2;if(d<distance){distance=d;nearest=row;}}atlasCanvas.title=nearest?`${nearest.feature||nearest.label||"对象"} · ${nearest.source_id||nearest.cell_id} · x=${nearest.x}, y=${nearest.y}${atlasData.layer==="transcripts"?` · z=${nearest.z??"未提供"} · qv=${nearest.qv??"未提供"} · 源归属=${nearest.cell_id??"未分配"}`:""}`:"";});
  atlasCanvas.addEventListener("pointerup",e=>{if(!atlasDrag)return;const p=atlasWorld(e),d=atlasDrag;atlasDrag=null;if(d.pan){atlasBounds=d.bounds.map((v,i)=>v+d.point[i%2]-p[i%2]);}else{if(Math.abs(p[0]-d.point[0])<1e-8||Math.abs(p[1]-d.point[1])<1e-8)return;atlasBounds=[Math.min(p[0],d.point[0]),Math.min(p[1],d.point[1]),Math.max(p[0],d.point[0]),Math.max(p[1],d.point[1])];}atlasAsync(loadAtlas)();});
  action("atlasFreeze",async()=>{if(!atlasCurrent()||!atlasBounds)throw new Error("请先读取数据集与视野。");const r=await call("freeze_atlas_roi",{atlas_id:atlasCurrent().object_id,bounds:atlasBounds});$("atlasFrozen").textContent=`已冻结 ${r.shape[0]} 对象 × ${r.shape[1]} 特征；可在整合与复核中选择输入 ${r.input_id}`;await refreshAnalysis();});
  async function refreshStudies(){const r=await call("list_study_analyses"),s=$("studyChoice"),old=s.value;s.replaceChildren();for(const d of r.study){const o=textNode("option",d.study_id);o.value=d.object_id;s.append(o);}if(r.study.some(d=>d.object_id===old))s.value=old;if(r.studyresult.length)renderStudy(r.studyresult.at(-1));}
  async function readStudyFile(id){const f=$(id).files[0];if(!f||f.size>2_000_000)throw new Error("请选择小于 2 MB 的 JSON 文件。");return JSON.parse(await f.text());}
  action("studyRegister",async()=>{const r=await call("register_study_design",{spec:await readStudyFile("studyDesignFile")});await refreshStudies();$("studyChoice").value=r.object_id;$("studyStatus").textContent="设计已冻结；登记本身不等于独立核实样本身份。";});
  function renderStudy(r){$("studyDetails").textContent=pretty(r);const holder=$("studyResults");holder.replaceChildren();for(const row of r.results){const box=textNode("div","","row");box.append(textNode("strong",row.metric),textNode("span",row.status==="tested"?`差异 ${fmt(row.effect)} · 95% CI ${row.confidence_interval.map(fmt).join(" 至 ")} · p=${Number(row.p_value).toPrecision(4)} · q=${Number(row.q_value).toPrecision(4)}`:`未检验：${row.reason}`));holder.append(box);}$("studyStatus").textContent=`以受试者为独立单位 · ${r.family_size} 项完整检验家族 · ${r.roi_inference==="exploratory_posthoc"?"事后探索 ROI":"已声明 ROI 设计"}`;}
  action("studyRun",async()=>{if(!$("studyScaleConfirmed").checked)throw new Error("请先核对指标尺度适用性。");const parse=id=>$(id).value.split(",").map(x=>x.trim()).filter(Boolean);const r=await call("run_study_inference",{study_id:$("studyChoice").value,records:await readStudyFile("studyRecordsFile"),plan:{method:$("studyMethod").value,control:$("studyControl").value,case:$("studyCase").value,roi_class:$("studyRoi").value,metrics:parse("studyMetrics"),covariates:parse("studyCovariates"),adjust_batch:true,missing_policy:$("studyMissing").value,outcome_scale:"continuous_section_summary",rationale:$("studyRationale").value}});renderStudy(r);});

  function showPage(next) {
    page = next;
    for (const [name, nav] of [["home", "navHome"], ["rna", "navRna"], ["protein", "navProtein"], ["integration", "navIntegration"], ["atlas", "navAtlas"]]) {
      show(`${name}View`, name === page); $(nav).className = name === page ? "active" : "";
      if ($(nav).setAttribute) $(nav).setAttribute("aria-current", name === page ? "page" : "false");
    }
    if (page === "rna") draw();
    if (page === "atlas") refreshAtlas().catch(e=>status(e.message,"error"));
    if (page === "protein") drawProtein();
    if (page === "integration") {refreshIntegrations().catch(e => status(e.message, "error"));refreshAnalysis().catch(e=>status(e.message,"error"));}
    if (state) updateContext();
  }
  function metric(holder, label, value) {
    const box = textNode("div", "", "metric"); box.append(textNode("span", label), textNode("strong", String(value))); holder.append(box);
  }
  function currentAssay() { return overview?.protein_assays.find(a => a.assay_id === $("proteinAssay").value); }
  function proteinFeatures() {
    const old = $("proteinFeature").value; $("proteinFeature").replaceChildren();
    for (const feature of currentAssay()?.features || []) {
      const option = textNode("option", `${feature.symbol} · ${feature.feature_id}`); option.value = `feature_id:${feature.feature_id}`; $("proteinFeature").append(option);
    }
    if ((currentAssay()?.features || []).some(f => `feature_id:${f.feature_id}` === old)) $("proteinFeature").value = old;
  }
  function renderOverview() {
    if (!overview || !state) return;
    $("homeProject").textContent = overview.name;
    $("homeScope").textContent = `${overview.sample_id || state.metadata.slice_id} · ${coordinateUnit()} · 当前版本 ${short(overview.head_revision)}`;
    const metrics = $("homeMetrics"); metrics.replaceChildren();
    metric(metrics, `空间对象 · ${overview.observation_unit}`, overview.observation_count);
    metric(metrics, "RNA 特征", overview.rna_feature_count);
    metric(metrics, "蛋白数据层", overview.protein_assays.length);
    metric(metrics, "共享选区对象", selected()?.cell_ids.length || 0);
    $("rnaSummary").textContent = `${overview.rna_feature_count} 个实测特征 · ${overview.revision_count} 个版本`;
    $("proteinSummary").textContent = overview.protein_assays.length ? overview.protein_assays.map(a => `${a.name} · ${a.features.length} 个通道`).join("；") : "尚未接入配对蛋白数据";
    const selectionText = selected() ? `共享 ROI：${selected().cell_ids.length} 个对象 · ${selected().selection_id} · ${selected().stale ? "历史版本，请重新选择后修订" : "当前版本"}` : "尚未选择 ROI；可进入转录组圈选区域，或在蛋白图中点击对象。";
    $("homeSelection").textContent = selectionText; $("proteinSelection").textContent = selectionText;
    $("proteinAvailability").textContent = overview.protein_assays.length ? `已接入 ${overview.protein_assays.length} 个蛋白数据层。配对依据为原始 ID、样本声明和精确坐标；${coordinateUnit()}。` : "当前项目没有蛋白测量。可通过 register-protein 接入同一样本的 H5AD 或宽表 CSV；需要原始对象 ID、坐标和测量类型。RNA 表达不能替代蛋白测量。";
    const old = $("proteinAssay").value; $("proteinAssay").replaceChildren();
    for (const assay of overview.protein_assays) { const option = textNode("option", assay.name); option.value = assay.assay_id; $("proteinAssay").append(option); }
    if (overview.protein_assays.some(a => a.assay_id === old)) $("proteinAssay").value = old;
    proteinFeatures();
    for (const id of ["proteinBase", "proteinTarget"]) {
      const previous = $(id).value; $(id).replaceChildren();
      for (const [i, revision] of (state.revisions || []).entries()) { const option = textNode("option", `${i ? `版本 ${i}` : "原始"} · ${short(revision.revision_id)}`); option.value = revision.revision_id; $(id).append(option); }
      if ((state.revisions || []).some(r => r.revision_id === previous)) $(id).value = previous;
      else if (id === "proteinTarget") $(id).value = state.head_revision;
    }
    $("homeRuns").replaceChildren();
    if (!overview.recent_runs.length) $("homeRuns").append(textNode("p", "尚无保存的分析。", "muted"));
    for (const saved of [...overview.recent_runs].reverse()) {
      const button = textNode("button", `${short(saved.run_id)} · ${saved.stale ? "历史版本" : "当前版本"}`, "quiet");
      button.addEventListener("click", async () => {
        if (busy) return; busy = true; enableActions();
        try { run = await call("get_run", { run_id: saved.run_id }); renderRun(run); showPage(run.analysis_schema === "spatial-collab.protein-region.v1" ? "protein" : "rna"); }
        catch (error) { status(error.message, "error"); } finally { busy = false; enableActions(); }
      }); $("homeRuns").append(button);
    }
    enableActions();
  }
  function clearProteinView(message = "显示条件已变化，请重新读取蛋白证据。") {
    proteinViewData = null; proteinTransform = null;
    $("proteinStats").replaceChildren(); $("proteinInterpretation").textContent = message;
    $("proteinMapTitle").textContent = "蛋白空间分布 · 尚未读取"; $("proteinRange").textContent = ""; $("proteinSource").textContent = "";
    drawProtein(); enableActions();
  }
  function drawProtein() {
    const map = $("proteinMap"), box = map.getBoundingClientRect();
    if (!box.width || !box.height) return;
    const dpr = window.devicePixelRatio || 1, context = map.getContext("2d");
    map.width = Math.round(box.width * dpr); map.height = Math.round(box.height * dpr); context.setTransform(dpr, 0, 0, dpr, 0, 0); context.clearRect(0, 0, box.width, box.height);
    const points = proteinViewData?.points || []; if (!points.length) return;
    const xs = points.map(p => p.x), ys = points.map(p => p.y), xmin = Math.min(...xs), xmax = Math.max(...xs), ymin = Math.min(...ys), ymax = Math.max(...ys);
    const scale = Math.min((box.width - 60) / (xmax - xmin || 1), (box.height - 60) / (ymax - ymin || 1));
    const ox = (box.width - (xmax - xmin) * scale) / 2, oy = (box.height - (ymax - ymin) * scale) / 2;
    proteinTransform = p => [ox + (p.x - xmin) * scale, oy + (p.y - ymin) * scale];
    const values = points.filter(p => p.included && p.display_value != null).map(p => p.display_value);
    const min = values.length ? Math.min(...values) : 0, max = values.length ? Math.max(...values) : 0, ids = selectedIds();
    $("proteinRange").textContent = values.length ? `${proteinViewData.scale} · ${fmt(min)} → ${fmt(max)}` : "没有有效信号";
    for (const point of points) {
      const [x, y] = proteinTransform(point), valid = point.included && point.display_value != null;
      const ratio = max === min ? 0.5 : (point.display_value - min) / (max - min);
      context.fillStyle = valid ? `hsl(${185 - ratio * 155},65%,${30 + ratio * 35}%)` : "#50565e";
      context.beginPath(); context.arc(x, y, points.length > 3000 ? 2 : 3, 0, Math.PI * 2); context.fill();
      if (ids.has(point.cell_id)) { context.strokeStyle = "#f6f5ff"; context.lineWidth = 0.7; context.stroke(); }
    }
  }
  async function readProtein() {
    const assay = currentAssay(); if (!assay) throw new Error("请先接入与当前样本对应的蛋白数据。");
    const revision = state.head_revision, selectionId = selected()?.selection_id || null;
    const result = await call("inspect_protein", { assay_id: assay.assay_id, revision_id: revision,
      feature: $("proteinFeature").value, selection_id: selectionId, scale: $("proteinScale").value, cofactor: Number($("proteinCofactor").value) });
    if (state.head_revision !== revision || selected()?.selection_id !== (selectionId || undefined) || result.assay_sha256 !== assay.assay_sha256) throw new Error("读取期间共享对象发生变化，请重新读取。");
    proteinViewData = result; const stats = result.summary, holder = $("proteinStats"); holder.replaceChildren();
    metric(holder, "已测 / 纳入", `${stats.measured_count} / ${stats.included_count}`); metric(holder, "缺失测量", stats.missing_count);
    metric(holder, "实测零信号", stats.zero_count); metric(holder, "原始均值", fmt(stats.mean_raw));
    $("proteinMapTitle").textContent = `${result.feature.symbol} · ${result.measurement_type === "antibody_count" ? "抗体计数" : "信号强度"}`;
    $("proteinInterpretation").textContent = `${selectionId ? "当前 ROI" : "全部导入对象"} · 原始中位数 ${fmt(stats.median_raw)} · ${result.scale} 均值 ${fmt(stats.mean_transformed)}。${result.view_complete ? "图显示全部导入对象；统计范围按共享 ROI。" : "超过 10,000 个对象，当前图省略；统计仍覆盖完整指定对象。"}显示变换不等于背景校正。`;
    $("proteinSource").textContent = pretty({ assay, inspection: { ...result, points: undefined } });
    drawProtein(); await updateContext(); status("蛋白信号已读取，缺失与实测零值分别计数。", "success");
  }
  function renderProteinRun(result) {
    proteinRun = result; const holder = $("proteinResults"); holder.replaceChildren();
    holder.append(textNode("p", `保存记录 ${result.run_id} · ${result.stale ? "历史版本" : "当前版本"} · ${result.parameters.scale} · 基线 ${short(result.base_revision)} → ${short(result.target_revision)}。以下是该次保存的参数与对象，修改输入不会改写结果。`, "notice"));
    const wrap = textNode("div", "", "table-scroll"), table = document.createElement("table"), header = document.createElement("tr");
    ["蛋白", "前景原始均值 前 → 后", "背景原始均值 前 → 后", "前景已测 / 缺失 前 → 后", "背景已测 / 缺失 前 → 后", "区域差（显示尺度）前 → 后", "区域差变化", "RNA 配对 n 前 → 后", "RNA 原始 ρ 前 → 后", "RNA 面板归一化 ρ 前 → 后"].forEach(x => header.append(textNode("th", x))); table.append(header);
    for (const feature of result.parameters.features) {
      const a = result.before[feature], b = result.after[feature], row = document.createElement("tr");
      const values = [feature, `${fmt(a.foreground.mean_raw)} → ${fmt(b.foreground.mean_raw)}`, `${fmt(a.background.mean_raw)} → ${fmt(b.background.mean_raw)}`,
        `${a.foreground.measured_count} / ${a.foreground.missing_count} → ${b.foreground.measured_count} / ${b.foreground.missing_count}`,
        `${a.background.measured_count} / ${a.background.missing_count} → ${b.background.measured_count} / ${b.background.missing_count}`,
        `${fmt(a.contrast)} → ${fmt(b.contrast)}`, fmt(result.comparison[feature].contrast_delta),
        a.paired_rna ? `${a.paired_rna.raw.n} → ${b.paired_rna.raw.n}` : "未指定 RNA",
        a.paired_rna ? `${fmt(a.paired_rna.raw.rho)} → ${fmt(b.paired_rna.raw.rho)}` : "未指定 RNA",
        a.paired_rna ? `${fmt(a.paired_rna.rna_panel_per_10000.rho)} (n=${a.paired_rna.rna_panel_per_10000.n}) → ${fmt(b.paired_rna.rna_panel_per_10000.rho)} (n=${b.paired_rna.rna_panel_per_10000.n})` : "未指定 RNA"];
      values.forEach(value => row.append(textNode("td", value))); table.append(row);
    }
    wrap.append(table); holder.append(wrap, textNode("p", `RNA：${result.parameters.rna_gene || "未指定"}。所有 ${result.parameters.features.length} 个通道均保留。注释名称改变可使组成变化而蛋白均值不变；完整记录包含标签分母。ρ 是前景内配对对象的描述性相关，无生物学重复推断。`, "small muted"));
    $("proteinRunDetails").textContent = pretty(result); enableActions();
  }
  async function compareProteins(all) {
    if (!selected()) throw new Error("请先在同一项目创建前景共享 ROI。");
    const assay = currentAssay(); if (!assay) throw new Error("尚未接入蛋白数据。");
    status("正在比较固定 ROI 的蛋白证据，并保存全部结果…");
    run = await call("run_protein_comparison", { assay_id: assay.assay_id, base_revision: $("proteinBase").value, target_revision: $("proteinTarget").value,
      selection_id: selected().selection_id, background_selection_id: $("proteinBackground").value.trim() || null,
      features: all ? assay.features.map(f => `feature_id:${f.feature_id}`) : [$("proteinFeature").value],
      scale: $("proteinScale").value, cofactor: Number($("proteinCofactor").value), rna_gene: $("pairedRnaGene").value.trim() || null });
    renderProteinRun(run); await refresh(); status("蛋白区域比较已保存，全部通道和未知结果均保留。", "success");
  }
  for (const [id, next] of [["navAtlas", "atlas"], ["navIntegration", "integration"], ["navHome", "home"], ["navRna", "rna"], ["navProtein", "protein"], ["openRna", "rna"], ["openProtein", "protein"], ["proteinToRna", "rna"], ["homeSelectionAction", "rna"]]) action(id, () => showPage(next));
  action("homeQcAction", async () => { showPage("rna"); quality = await call("inspect_quality", { revision_id: state.head_revision }); renderQuality(quality); status("已读取 RNA 测量质量。", "success"); });
  action("homeExport", async () => { const result = await call("export_review_bundle", { run_id: run?.run_id || null, compact: true }); $("exportResult").textContent = pretty(result); show("exportResult", true); showPage("rna"); status(`审阅记录已导出：${pretty(result)}`, "success"); });
  action("proteinInspect", readProtein);
  action("proteinCompare", () => compareProteins(false)); action("proteinCompareAll", () => compareProteins(true));
  action("proteinExport", async () => { if (!proteinRun) throw new Error("请先保存蛋白分析。"); const result = await call("export_review_bundle", { run_id: proteinRun.run_id, compact: true }); status(`蛋白审阅包已导出：${pretty(result)}`, "success"); });
  $("proteinAssay").addEventListener("change", () => { proteinFeatures(); clearProteinView(); });
  for (const id of ["proteinFeature", "proteinScale", "proteinCofactor"]) $(id).addEventListener("change", () => clearProteinView());
  action("proteinSelectVisible", async () => { if (!proteinViewData?.points.length) throw new Error("请先读取空间信号图。"); const ids = proteinViewData.points.map(p => p.cell_id); await select({ cell_ids: ids }); await readProtein(); });
  $("proteinMap").addEventListener("click", async event => {
    if (busy || !proteinTransform || !proteinViewData) return;
    const box = $("proteinMap").getBoundingClientRect(), x = event.clientX - box.left, y = event.clientY - box.top;
    let nearestPoint = null, distance = 100;
    for (const point of proteinViewData.points) { const p = proteinTransform(point), d = (p[0] - x) ** 2 + (p[1] - y) ** 2; if (d < distance) { distance = d; nearestPoint = point; } }
    if (!nearestPoint) return; busy = true; enableActions();
    try { await select({ cell_ids: [nearestPoint.cell_id] }); await readProtein(); }
    catch (error) { status(error.message, "error"); } finally { busy = false; enableActions(); }
  });
  window.addEventListener("resize", drawProtein);

  async function initialize() {
    if (local) $("connection").textContent = "● 本地共享项目";
    else {
      if (window.parent === window) throw new Error("请在兼容 MCP Apps 的宿主中打开，或启动本地工作台。");
      const init = await rpc("ui/initialize", { appInfo: { name: "Spatial Collab Workbench", version: "0.8.0-alpha.1" },
        appCapabilities: { availableDisplayModes: ["inline", "fullscreen"] }, protocolVersion: "2026-01-26" });
      if (init.protocolVersion !== "2026-01-26") throw new Error("宿主协商的 MCP Apps 版本尚未受此 Alpha 支持。可继续使用工具或本地工作台。");
      hostCapabilities = init.hostCapabilities || {};
      notify("ui/notifications/initialized", {}); $("connection").textContent = "● MCP Apps 共享项目";
      new ResizeObserver(() => notify("ui/notifications/size-changed", { height: document.documentElement.scrollHeight })).observe(document.body);
    }
    await refresh(); refreshAssets().catch(() => {}); status("项目已加载。先选择区域，再检查标记、提出修订并检验影响。", "success");
    if (state.metadata.source_kind === "atlas_workspace") {
      await showPage("atlas"); status("大数据工作区已加载。选择数据集后浏览全切片、分子与图像，或冻结区域进行分析。", "success");
    }
    pollTimer = setInterval(async () => {
      if (tornDown || busy || polling || document.hidden || dragging || polygon.length) return;
      polling = true;
      try { const context = await call("get_context", { after_cursor: contextCursor });
        if (context.changed) { await refresh(); status("已同步其他界面的共享选区、修订或分析记录。", "success"); }
      } catch (error) { status(`共享同步暂不可用：${error.message}`, "error"); }
      finally { polling = false; }
    }, local ? 2500 : 5000);
  }
  window.SpatialCollab = {
    call,
    currentProject: () => activeProjectId,
    openProject: async projectId => {
      if (busy || polling) throw new Error("共享状态正在更新，请稍后切换项目。");
      busy = true; enableActions();
      try {
      activeProjectId = projectId;
      state = null; proposal = null; run = null; viewBounds = null; contextCursor = null;
      plan = null; planDirty = false; planContextHash = null; proteinViewData = null; proteinRun = null;
      planBinding = null; planInitialized = false; planList = []; variantEditors = [];
      polygon = []; assetList = []; assetInspection = null; analysisPage = null;
      await refresh(); await refreshAssets(); await showPage("home");
      } finally { busy = false; enableActions(); }
    },
    fullscreen: () => rpc("ui/request-display-mode", { mode: "fullscreen" }),
  };
  initialize().catch(error => status(error.message, "error"));
})();
