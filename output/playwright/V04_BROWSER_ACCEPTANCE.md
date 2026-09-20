# v0.4 本机实际浏览器验收

2026-09-20，Playwright CLI 独立浏览器会话 `spatial-v04`。这些是实际页面的点击、输入、载入与渲染记录，不是 ChatGPT iframe 或 Claude 模型对话验收。

| 实际页面 | 操作与观察 | 记录 |
|---|---|---|
| 淋巴结 7,725 对象 | 载入原冻结 v3、读取并运行全部 24 格；18 computed / 6 unknown，差异排序、双估计量与分母可见 | `v04-loaded-plan.yml`、`v04-frozen-run.yml` |
| 淋巴结缺图 | 读取资产清单，0 个资产；页面明确缺图/分割，不判合并或真实共表达 | `v04-lymph-overview.png` |
| 新建计划 | 实际表格输入 `counts.CD3D eq 0` 数字条件、假设排除、15 µm/两范围；新草稿 v1→冻结 v2→执行 | `v04-draft-saved.yml`、`v04-unknown-plan-result.yml`、`v04-hypothesis-editor.png` |
| 新计划结果 | 2 格均 unknown，每格 7,725 个未知条件成员；没有注释提交 | 新计划 `plan_ffb19995fe384dd5bfea60c324db592c`；原分析计划 v3 保留 |
| HER2ST A1 | 实际查看 ERBB2/EPCAM/MS4A1 与 QC；MS4A1 显示未测，之后物理比较仍 disabled，区域比较 enabled | `v04-her2-marker-qc.yml` |
| Kidney tiny | 从原生登记清单选择 mask，输入 `mnlmjgpm-1,ipcgjdge-1` 建立选区；显示一个对应另一声明 cell、一个对应本 cell，二者均标签 221；不作生物学判决 | `v04-kidney-correspondence.yml` |
| 旧 Xenium 特征查询 | 页面默认明确 feature ID 查询，修复后五个特征均显示面板已测量（含正确的零） | 同上；新增核心回归 `test_legacy_source_accepts_explicit_catalog_queries_without_rewriting_source` |
| 安装版 SPOTS1 | 用最终独立 wheel 安装提供页面，2,568 spots/32,285 features；Ptp4a1 显示两个候选 ID 的歧义，两个明确 ID 分别是已测量零；目录保留两行 | `v04-spots-ambiguous-marker.yml`、`v04-spots-feature-catalog.yml`、`v04-spots-identity.png` |

淋巴结页面实际测得 `scrollWidth <= innerWidth`，运行/缺图检查时 console 0 errors。中途重启 kidney 本地服务后，旧网页会话出现预期的 CSRF 403；刷新页面取得新会话后恢复。原生 CUA/Node 浏览器工具曾因本机 kernel assets 路径不可用而失败，因此改用独立 Playwright CLI；这些失败不算通过。

浏览器验收增加了明确命名的计算计划、运行记录和 kidney 的精确共享选区，没有改动表达源、细胞注释或分割。没有声称两种方法中的任一种更准确，也没有科研人员效率或外部宿主采用证据。
