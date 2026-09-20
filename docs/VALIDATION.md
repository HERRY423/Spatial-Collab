# Spatial Collab 0.5.0-alpha.1 多组学扩展验收

日期：2026-09-20。真实数据结果及范围见 [蛋白复核报告](V05_PROTEIN_REVIEW.md)。

| 检查 | 结果 | 证据范围 |
|---|---|---|
| 主环境回归 | 442 passed / 2 skipped / 28 warnings | 可选原生测试另在 napari 环境运行 |
| 最终计算复用与显示分母改动复验 | 34 passed / 1 warning | 蛋白、离线重放、生产 JS 与协议；`projects/v05-final-targeted.xml` |
| 原生 napari 相关回归 | 47 passed / 43 warnings | 既有适配、图层、图像/分割功能；未新增蛋白 dock |
| 真实 RNA/ADT | 5,336 个 spots，2 × 21 个蛋白通道 | 所有 ID/XY、蛋白矩阵、区域均值、零值及相关性核对；2 个包 matched |
| 独立安装 | 25 个运行文件三方 SHA256 一致；31 个 MCP 工具 | 源码 / wheel / `projects/installed-v05-final`；Python `-I` 模块路径核验 |
| 安装版重放 | 3 个 matched | 两片 SPOTS 蛋白、旧 HER2ST RNA，兼容旧包 |
| 浏览器 | 主页、CD19 空间信号、全 21 通道结果、RNA 页面切换、RNA-only 缺失态 | Chromium 真实本地服务；不等于模型宿主对话验收 |

wheel：`dist/spatial_collab-0.5.0a1-py3-none-any.whl`，SHA256 `36baee72004aad17e6a8a0f39a994936282978d2c8c6e33a29879095b1028fa7`。独立安装收据：`projects/v05-audit/installed-delivery-receipt.json`；全部 25 个文件在执行前后均与源码、wheel 一致。

本轮最初沙箱创建测试目录被拒，随后在明确当前工作区范围内重跑通过；napari 环境缺少 MCP 的测试收集失败单独保留，相关测试在主环境通过。不存在部署、科学授权、独立研究者验收或多样本整合通过的声明。

---

# Spatial Collab 0.4.0-alpha.1 四项修复验收（历史）

2026-09-20。具体科学问题、真实数据数值与边界见 [四项复测报告](V04_REAL_DATA_REVIEW.md)。以下是包含最后旧项目特征前缀修复的交付检查；旧 candidate 与旧版收据继续保留。

| 检查 | 实际结果 | 范围 |
|---|---|---|
| 最终全部基础回归 | **431 passed，2 skipped，27 warnings；41.28 秒** | 跳过的两项要求真实 napari 运行时；其余包括身份/单位/图像核心、假设账本、数值、并发、重算、协议、实际服务驱动的 JS DOM 工作流 |
| 最终原生 napari 环境 | **47 passed，0 skipped；13.85 秒** | Python 3.12.13 / napari 0.8.0，原生隐藏 Image/Labels/Points/dock；包含最后核心兼容修复，执行前后源码一致 |
| 静态与插件清单 | Ruff、Codex plugin validator、Claude 原生插件校验、skill validator 通过 | 本地源码与清单格式；没有市场发布或外部账户安装 |
| 真实输入保真 | SPOTS 两片与 HER2ST A1 的完整源计数/坐标/身份一致；原文件 SHA 不变 | 两片 RNA 共 5,336 spots，各 32,285 features；HER2ST 346×12,562；不读封存答案 |
| 图像/掩膜真实格式 | Kidney tiny 358 对象，同坐标加载；357 同声明标签、1 不同声明标签，保留差异 | 检查对应关系，不是分割准确性或真实共表达证据 |
| 淋巴结声明假设 | 24 格全保留：18 computed、6 unknown；两种估计量前后均独立数值核对一致 | 单片描述性工作标签，非独立生物学验证；零阈值翻转不代表数值和分母未变 |
| 实际本机浏览器 | 淋巴结读取/运行 24 格；表格新建→冻结→执行 2 格 unknown；HER2ST marker/QC 后仍禁物理分析；kidney 错对应保留；安装版 SPOTS 重复符号与明确 ID 正确显示 | 独立 Playwright 浏览器；不是实际 ChatGPT iframe 或 Claude 模型任务对话 |
| 最终独立安装 | **24 个运行时文件**源码/wheel/安装逐个散列一致，`python -I` 从独立目录加载，版本 0.4.0a1，实际发现 **27 个 MCP 工具** | `projects/installed-v04-final`，运行前后文件未变 |
| 安装版真实重算 | **五份 matched**：v0.4 淋巴结 compact、HER2ST compact、两片 SPOTS compact、旧 v0.3 完整包 | 重建同一来源、注释、对象与参数；不授予科学授权 |

最终 wheel：`dist/spatial_collab-0.4.0a1-py3-none-any.whl`；SHA256 `93b8d7e03cc3a8e24331280693786b7978e70d85ef038af07897d584a50384aa`。首个 candidate 在源码/安装散列检查时发现与最后修复不同，主动停止，标为 superseded；未将它计作通过。

本机可审阅收据：

- `projects/v04-audit/final-source-checks-20260920T121050/receipt.json`，附完整测试日志与 XML、70 个源码相关文件前后散列。
- `projects/validation/napari-v04-final-legacy-20260920/acceptance.json`，附真实环境、关键核心/图像源码散列与完整日志。
- `projects/v04-audit/installed-delivery-receipt.json`，附五份真实 bundle 的重算结果、工具发现与安装路径核对。
- `output/playwright/V04_BROWSER_ACCEPTANCE.md`、相邻实际快照和 PNG。浏览器的独立验收计划/运行与共享选区均明确记录，没有提交注释。
- `projects/v04-audit/INDEPENDENT_REVIEW.md`，包括错误选区绑定、类型比较及修复复验；恶意 bundle 在开启或关闭重算时均拒绝。

失败和告警没有删除：初次 napari 只读缓冲区与渲染不兼容、旧特征前缀误报未测、原生浏览器工具 kernel assets 不可用以及本机服务重启后旧会话的 CSRF 403 均已记录并有后续可执行验证。基础 27 条 warning 主要为第三方弃用和刻意重复符号输入；原生环境另有 43 条第三方告警和既有无关插件缺模块日志。

没有外部发布、远端部署、完整组织性能、主淋巴结同源图像、HER2ST/SPOTS 物理标定、多患者推断或研究者效率对照。这些本地通过不算 SpatialBench 官方分数、真实模型宿主采用或科学结论成立。以下为历史状态。

---

# Spatial Collab 0.3.0-alpha.1 真实数据驱动验证

日期：2026-09-20。完整科学范围、失败与局限见 [真实数据分析报告](REAL_DATA_STUDY.md)，原始来源审查见 [输入审计](REAL_DATA_INPUT_AUDIT.md)。下面保留旧版记录，历史成功不算本版新证据。

| 本轮检查 | 实际结果 | 证据范围 |
|---|---|---|
| 全部基础环境回归 | **308 passed, 1 skipped, 4 warnings，36.96 秒** | 包括新窗口导入／缺失 QC／属性筛选／非提交假设／重放，及既有格式、协作、MCP；skip 为基础环境未装 napari |
| napari 原生环境 | **22 passed, 43 warnings，10.89 秒** | 真实 napari/npe2 隐藏窗口和既有跨客户端流程；不是科研人员体验试验 |
| 原始真实数据导入 | **7,725 个对象、4,624 个基因、2,140,795 个转录本** | 全窗口，无抽样；全部细胞表达总数与源 CSV 一致；原始哈希不变 |
| 独立原始列与坐标检查 | 全部 ID/XY、源 QC 与预先固定 10 个表达列一致 | 用原 H5/CSV 直接读出，不以插件结果生成期待值；更正了初次调用的错误组织 ID |
| 独立统计实现 | **36 条前后记录全部一致** | 直接坐标距离遍历，未调用插件分析函数/KDTree；不是外部生物学验证 |
| 真实导出离线复算 | 两种图范围的审阅包均 matched | 固定源、基线、选区、具体假设、参数和数值；不提供科学授权 |
| 实际本机浏览器 | 真实项目来源、3,383 对象共享选区、QC、未测量 marker、9 行假设表读取完成 | 页面显示未经验证的自动工作标签，无假设提交；没有横向页面溢出 |
| 原生清单与静态检查 | Codex、Claude、技能校验与 Ruff 通过 | 本地格式／代码检查；没有外部发布 |

本次真实分析产生于升级过程中的工作区源码，原始结果诚实记录当时的 source/distribution 版本字段为 0.2.0a1，不改写历史 provenance。交付版本为 0.3.0a1；安装包的独立复算验收另记，不能把旧安装版本字段当成当时运行了已发布 0.3 包的证明。

最终 UI 分母展示及协议定向回归为 **73 passed, 2 warnings，10.41 秒**。最终 wheel 在 `projects/installed-v03` 隔离安装，`python -I` 确认只从该安装目录导入：版本 0.3.0a1，**20 个运行时文件与当前源码及 wheel 逐个哈希一致**；18 个工具可发现，真实 QC 与属性筛选正确，两个真实敏感性 bundle 均重算 matched，项目状态没有变化。receipt 为 `projects/real-data-review/installed-delivery-receipt.json`。sdist 已包含新 DOM 回归所需的 `.cjs` 夹具；旧构建未包含该成员的问题已修正。

4 条基础警告与旧版相同类别：Starlette/httpx 弃用、刻意重复特征夹具、ome_zarr Scaler。napari 的 43 条为现有 npe2/Pydantic 弃用提示。局部验证不与完整回归相加。

未完成：整张高密度组织、多患者推断、形态／分割证据、标尺缺失的数据全面导入、科研人员效率对照、实际 ChatGPT/Claude 模型任务对话、远端部署和发布。本次真实数据运行不是 SpatialBench 官方评测。

---

# Spatial Collab 0.2.0-alpha.1 历史验证记录

日期：2026-09-20。本轮为本地研究原型扩展，下面保留 v0.1 历史验收，不能把两个版本的数字相加当作独立验证。

| 本轮层次 | 实际结果 | 范围 |
|---|---|---|
| 完整基础环境回归 | **215 passed, 1 skipped, 4 warnings，31.67 秒** | 格式、状态、分析、导出、重算和真实 MCP stdio/HTTP；skip 为基础环境未安装 napari |
| napari 原生环境 | **22 passed，43 上游告警，10.40 秒** | Python 3.12.13 / napari 0.8.0 / npe2 0.8.3，Windows 原生隐藏双窗口、真实 plugin dock、worker、Agent 提案交接 |
| 格式专项 | **61 passed** | AnnData、Xenium H5/MTX、10x MTX、标准 Visium、CosMx、MERSCOPE 合成格式夹具 |
| SpatialData 专项 | **10 passed** | 真实 SpatialData 0.8.0 序列化 Points/Shapes；非单位变换、乱序关联、源文件不变及错误输入拒绝 |
| 协作任务 | **8 / 8** | `scripted_reference_driver`，开发者编写合成题；没有运行 SpatialBench 或模型评测 |
| 原生清单/技能/静态检查 | 全部通过 | 官方 Codex validator、Claude plugin validate、skill validator、Ruff、JavaScript 语法 |
| 旧产物兼容 | v0.1 实际导出包重算 matched | 未重写旧快照或散列；兼容原指标与 bundle 格式 |
| 浏览器跨客户端验收 | 实际点击、填写和服务端交接完成 | 外部工具接口提出两对象提案→网页自动收到→测试署名确认→区域比较→保存→导出 |
| Python 包 | wheel + sdist 成功构建；独立目录安装检查通过 | `python -I` 仅插入独立安装路径，确认 7 格式、网页、napari manifest/entrypoint、MCP 和实际浏览器 bundle 重算 |

上述专项存在重叠，不相加报告。完整回归后针对最后的观测单位展示修改另跑核心/区域/任务/重算检查：**67 passed，9.78 秒**；实际浏览器导出包的 11 文件、9 组结果字段重算 matched。四条基础环境 warning 为 Starlette 测试客户端、httpx 请求体、刻意重复 var 名的 AnnData 夹具、ome_zarr Scaler 弃用提示。

跨客户端浏览器例子位于 `projects/demo-v02`：2 个合成对象从 T cell / Myeloid 改为 Uncertain，背景保持另外 478 个对象。标签比例差分别改变 −0.5，原始表达和 marker 对比保持不变；UNMEASURED 前后均不可判定。确认者 `Synthetic UI acceptance` 是自动化验收标识，不是真实研究者科学确认。导出目录为 `projects/demo-v02/exports/review_1d7c41fb55b54666b6d0b89910427a8f`。合成任务报告位于 `projects/benchmarks-v02/synthetic-collaboration-67de2b1be567/report.json`。

发现并修复：SpatialData 的 Parquet 分区目录被误当文件；浏览器 marker 响应选区与页面选区不同；失败的新预览留下旧确认意图；真实 napari shear 重置 API 不匹配；跨 FOV 身份/矩阵标识及原始计数边界。检查保留失败轨迹，不通过改答案掩盖错误。

napari `offscreen` 后端在本机不能创建 OpenGL，初次尝试失败；后改用 `viewer(show=False)` 的 Windows 原生后端。测试遇到系统临时目录权限问题，改用独立的新工作区临时目录后完成。已有旧 spatial-evidence-layer entrypoint 缺失模块提示未修改，详细告警与复现见 [NAPARI.md](NAPARI.md)。隐藏窗口验收不能写成可见界面人工体验或效率试验。

当前 16 个 MCP 工具，共享状态 token、分页游标、具体提案和结果均引用持久对象。状态 token 不是事件流，网页以轮询同步，napari 通过显式刷新同步；没有长任务调度、取消或自动恢复作业引擎。网页自动发现当前选区的最新提案，napari 可显式加载指定提案及其原始选区。

**未完成的验收**：真实组织的导入横评与科学结果、ChatGPT 实际会话/iframe、Claude 模型会话选择工具、公开远程服务、独立研究者体验与时间收益、全切片高密度性能、Visium HD 原生支持、组织图像/分割编辑、多患者推断、外部 CI 和发布。合成夹具、脚本任务、文件散列和隐藏窗口检查不能替代这些证据。

本轮适配器能力和具体格式边界见 [FORMATS.md](FORMATS.md)，公开基准来源与人机对照协议见 [BENCHMARK.md](BENCHMARK.md)。下面的 v0.1 限制是历史状态；涉及 SpatialData/napari 的旧“不支持”描述已由本轮有界实现更新。

---

# Spatial Collab 0.1.0-alpha.1 历史验证记录

日期：2026-09-20。本记录对应本地新建源码与安装包，不是远端 CI、正式发布、独立科学评审或真实用户试用。

## 本次通过的验证

| 层次 | 实际结果 | 能证明什么 |
|---|---|---|
| 自动化测试 | **119 passed, 3 warnings，12.21 秒** | 本机列出的输入、状态、计算与协议断言通过 |
| 静态检查 | `ruff check --no-cache src tests scripts` 通过 | 当前范围内无已启用规则问题 |
| Codex 插件格式 | 官方 `validate_plugin.py` 通过 | 插件清单与资源引用符合该校验器 |
| Claude 插件格式 | 本机 `claude plugin validate .` 通过 | Claude 插件清单通过原生验证 |
| 工作流技能 | 官方 `quick_validate.py` 通过 | 技能元数据格式通过 |
| MCP stdio | 真实子进程 + 官方 ClientSession 完成工作流 | 初始化、11 工具发现、选择、预览、确认、比较、导出及过期拒绝 |
| MCP Streamable HTTP | 真实本地服务 + 官方客户端完成同一工作流 | HTTP 协议可运行；不是云端宿主验收 |
| MCP Apps 资源 | 工具关联、资源读取、MIME 与元数据通过 | 标准 UI 资源可发现；不是实际 ChatGPT iframe 验收 |
| 本地浏览器 | 本机内嵌浏览器实际点击/输入/拖动 | 精确 ID、44 细胞矩形框选、18 细胞多边形圈选、标记、确认、前后结果 |
| 合成重放 | 导出后重新计算一致 | 输入、修订与指定描述性计算可重放 |
| Python 安装包 | wheel + sdist 构建成功，独立目标目录安装 | 从安装路径导入成功，包含静态页面，MCP 服务创建及离线重放成功 |

三条测试 warning 保留：Starlette 对测试客户端 httpx 的弃用提示、特意构造重复基因名的 AnnData 输入警告、httpx 测试请求的 raw-body 弃用提示。均未掩盖失败；原型没有以跳过检查代替通过。

## 可手算的合成例子

`python scripts/acceptance.py <新输出目录>` 创建 5 个细胞，其中一个来源细胞只有一个近邻。基线近邻标为 Myeloid，背景有 2/4 个目标细胞：观测比例 1，背景比例 0.5，超背景比例 **0.5**。

研究者测试记录把这个近邻改为 Stromal 后，观测比例 0，背景剩 1/4，超背景比例 **−0.25**。预设阈值 0.2，判断从达到变为未达到，差值 **−0.75**。所有结果只描述这个合成例子。

验收同时检查：未测量与零计数区别、预览不改 head、旧运行标为历史、过期提交被拒绝、原版本保留、恢复生成新版本、11 个导出成员验证通过、8 组结果字段重新计算一致。重放相对/绝对容差均为 `1e-12`。

本机验收输出：`projects/acceptance-20260920/acceptance.json`；审阅包：`projects/acceptance-20260920/exports/review_2d5a41dfbabc4fb3a7491886cd66a52f/`。源码原始输入散列：`0c62df0224f0601df12ba36fe884e336f2fab777f2fffaf94a4058da0c62fcbb`。

浏览器另验证了 4 个合成细胞改为 Uncertain 后，原始 ROI 的来源类型消失，返回 **indeterminate**；修改之外的另一个 ROI 在相同参数下结果稳定。它们分别检验不可判定与局部不受影响的情况，不是独立真实数据验证。

## 已覆盖的关键失败行为

- 两个客户端并发修改和过期 selection/proposal 不得静默覆盖；旧证据在界面明确失效。
- 缺失 panel 基因不当作零值；未知细胞 ID、重复 ID、混切片、非微米输入拒绝。
- Xenium 矩阵元数据尺寸不一致拒绝；条码顺序不同正确对齐；控制特征排除。
- 非整数、负数、非有限计数拒绝；稀疏无符号整数重复值不发生环绕溢出。
- ROI 外邻居的定义显式分开；排除后图、丰度与分母重算；没有来源/目标/边保留不可判定。
- 大邻接图在构造邻居索引前检查预算；视野超限不返回伪装成全体的抽样。
- 自交/退化多边形拒绝；大坐标平移不破坏正常面积判断。
- 导出路径约束、成员校验、版本链、提案与提交对应关系、重算差异检查。
- 浏览器跨 Origin/Host、缺失 CSRF、非 JSON、过大请求被拒绝；MCP 不接受字符串冒充布尔确认或额外任意路径参数。

## 证据边界与已知限制

**没有完成**：真实 Xenium 组织输入验收、ChatGPT 实际测试连接、Claude 实际模型对话/工具选择、MCP Apps 在真实宿主 iframe 中的运行、远端 CI、外部研究者/患者数据验证、生产力对照、发布或上架。

使用人工构造的 AnnData/Xenium 格式夹具测试适配器，不能写成“已在官方真实组织数据集上验收”。当前 480 细胞演示与 5 细胞手算案例全部是合成数据。

当前查看器只显示质心和计数摘要，没有组织切片、单转录本、细胞边界、分割编辑、SpatialData Zarr 直接适配或 Vitessce 集成。单切片指标没有空间置换检验、独立生物学重复、多患者比较或因果推断。

本次验证运行于 Windows/Python 3.13.9、MCP 1.28.1、NumPy 2.3.5、SciPy 1.16.3、Pydantic 2.13.4、Starlette 1.3.1；完整运行环境记录随比较结果导出。依赖声明范围的其他组合没有逐一验证。

具名确认只是归属记录；`scientific_authorization` 始终为 `NOT_ESTABLISHED`。散列验证与本软件重算不能证明第三方真实性或独立科学重复。

源码检查表见 `source-manifest.json`。可分发开发插件 ZIP 由 `scripts/package_plugin.py` 生成，排除 `projects/`、本地数据库和安装测试目录。它不是市场发布包已经获审的证明。
