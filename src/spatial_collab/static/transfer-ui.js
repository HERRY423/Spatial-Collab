"use strict";
(() => {
  const $ = id => document.getElementById(id);
  const bridge = window.SpatialCollab;
  let working = false, downloadUrl = null;
  const setStatus = text => { $("transferStatus").textContent = text; };
  const hex = bytes => Array.from(new Uint8Array(bytes), n => n.toString(16).padStart(2, "0")).join("");
  const digest = async bytes => hex(await crypto.subtle.digest("SHA-256", bytes));
  const encode = bytes => { let s = ""; for (let i = 0; i < bytes.length; i += 8192) s += String.fromCharCode(...bytes.subarray(i, i + 8192)); return btoa(s); };
  const decode = text => Uint8Array.from(atob(text), c => c.charCodeAt(0));
  const memory = {
    get: key => { try { return sessionStorage.getItem(key); } catch { return null; } },
    set: (key, value) => { try { sessionStorage.setItem(key, value); } catch { /* Host may disable iframe storage. */ } },
    remove: key => { try { sessionStorage.removeItem(key); } catch { /* Optional cache. */ } },
  };
  function action(id, fn) {
    $(id).addEventListener("click", async () => {
      if (working) return;
      working = true;
      for (const name of ["uploadAndImport", "downloadReview", "openPrivateProject", "reloadProjects"]) $(name).disabled = true;
      try { await fn(); } catch (error) { setStatus(error.message || String(error)); }
      finally { working = false; for (const name of ["uploadAndImport", "downloadReview", "openPrivateProject", "reloadProjects"]) $(name).disabled = false; }
    });
  }
  async function reloadProjects() {
    const data = await bridge.call("list_projects", { project_id: null });
    $("privateProjects").replaceChildren();
    for (const project of data.projects) {
      const option = document.createElement("option"); option.value = project.project_id;
      option.textContent = `${project.name} · ${project.cell_count} 个观测对象`;
      $("privateProjects").append(option);
    }
    if (bridge.currentProject()) $("privateProjects").value = bridge.currentProject();
  }
  action("reloadProjects", reloadProjects);
  action("openPrivateProject", async () => {
    const id = $("privateProjects").value;
    if (!id) throw new Error("请先刷新项目列表并选择项目。");
    await bridge.openProject(id); setStatus("已打开所选项目；后续操作均绑定该项目。");
  });
  action("fullScreenView", async () => { await bridge.fullscreen(); });
  action("uploadAndImport", async () => {
    const file = $("uploadFile").files?.[0];
    if (!file) throw new Error("请先选择文件。");
    if (!file.size || file.size > 128 * 1024 ** 2) throw new Error("请选择 1 byte–128 MiB 文件；更大数据请使用本地导入。");
    const options = file.name.endsWith(".h5ad") ? JSON.parse($("uploadOptions").value || "{}") : {};
    const sha256 = await digest(await file.arrayBuffer());
    const key = `spatial-upload:${sha256}:${file.size}`;
    let fileId = memory.get(key), record;
    if (fileId) {
      try { record = await bridge.call("get_transfer", { file_id: fileId, project_id: null }); }
      catch { fileId = null; }
    }
    if (!fileId) {
      record = await bridge.call("begin_upload", { filename: file.name, size: file.size, sha256, project_id: null });
      fileId = record.file_id; memory.set(key, fileId);
    }
    let offset = record.offset || 0;
    while (offset < file.size) {
      const bytes = new Uint8Array(await file.slice(offset, offset + 262144).arrayBuffer());
      record = await bridge.call("append_upload", { file_id: fileId, offset, data_base64: encode(bytes), project_id: null });
      offset = record.offset; setStatus(`已传输 ${Math.round(100 * offset / file.size)}%；尚未修改任何原项目。`);
    }
    await bridge.call("complete_upload", { file_id: fileId, project_id: null });
    const imported = await bridge.call("import_uploaded_project", { file_id: fileId, options, project_id: null });
    memory.remove(key); await reloadProjects(); await bridge.openProject(imported.project_id);
    $("privateProjects").value = imported.project_id;
    setStatus(`新项目已创建并打开。原文件 SHA256：${sha256}。`);
  });
  action("downloadReview", async () => {
    const projectId = bridge.currentProject();
    const file = await bridge.call("prepare_download", { compact: true, project_id: projectId });
    if (file.size > 128 * 1024 ** 2) throw new Error("审阅包超过当前浏览器 128 MiB 下载预算；完整原件已保存在私有项目 exports 目录，请在本地取回。");
    const chunks = []; let offset = 0;
    while (offset < file.size) {
      const part = await bridge.call("read_download", { file_id: file.file_id, offset, project_id: projectId });
      const bytes = decode(part.data_base64);
      if (!bytes.length || part.next_offset !== offset + bytes.length || part.sha256 !== file.sha256) throw new Error("下载分块不一致，请重新准备审阅包。");
      chunks.push(bytes); offset = part.next_offset;
      setStatus(`已取回 ${Math.round(offset * 100 / file.size)}%，正在核对内容。`);
    }
    const blob = new Blob(chunks, { type: "application/zip" });
    if (blob.size !== file.size || await digest(await blob.arrayBuffer()) !== file.sha256) throw new Error("下载校验失败；未提供保存链接。");
    if (downloadUrl) URL.revokeObjectURL(downloadUrl);
    downloadUrl = URL.createObjectURL(blob);
    const link = $("saveReviewFile"); link.href = downloadUrl; link.download = file.filename; link.hidden = false;
    setStatus(`文件已取回并通过 SHA256 校验：${file.sha256}。请点击“保存已校验的审阅包”。若宿主限制下载，可在本地工作台保存同一审阅包。`);
  });
})();
