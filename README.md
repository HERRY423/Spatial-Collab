# Spatial Collab · 空间多组学协作插件

**让研究者与 ChatGPT / Claude Code 在同一份空间数据上共同探索、修订注释，并检验修订对指定分析判断的影响。**

**v0.6 开发工作区新增：**正式整合结果接入、RNA-only / protein-only / 联合基线、官方 SMOPCA 小规模参考适配、“整合与复核”网页、固定模型筛选与新版本重新拟合、后台任务/缓存、分块蛋白矩阵与视口分页、显式模态/分子关系和个体层级研究设计。已加入第二种公开配对数据及独立用户试用包。详细操作、固定环境、方法边界和逐项交付见 [v0.6 实现说明](docs/V06_IMPLEMENTATION.md)。v0.5/v0.6 内容保留为历史基线说明。

**v0.7 新增可执行科学分析：**MOFA+/MEFISTO、参考去卷积、Harmony、PASTE、空间结构域、Moran 置换/FDR、数据库空间配体受体和 PROGENy 通路活性；附不可变稀疏输入、后台任务恢复、结果地图及基因编号映射。方法适用范围、安装及操作见 [v0.7 科学分析](docs/V07_SCIENTIFIC_WORKFLOWS.md)，实际运行证据见 [v0.7 验证](docs/V07_VALIDATION.md)。

**保持当前版本号的大数据与统计改造：**新增磁盘空间索引、完整密度金字塔、单分子与细胞/核边界、原生 tiled OME-TIFF 瓦片联动，以及受试者级 Welch、配对检验和含批次项的随机截距 LMM。新增“全切片与研究统计”页面；全景浏览与有界 ROI 分析分别执行。见 [操作与适用范围](docs/SCALABILITY_AND_STUDY.md) 和 [规模、像素、公开计数及统计验收](docs/SCALABILITY_VALIDATION.md)。

当前开发版本：**0.7.0-alpha.1，受限研究原型，尚未发布**。以下为 v0.5 基线能力：主页汇总当前样本、测量层、共享 ROI、版本和已保存分析；RNA 与配对蛋白使用同一组空间对象。原始文件保持不变，研究者和模型共享项目数据库中的选区、版本、具体修订提案和结果，不以聊天历史为数据状态。

新增蛋白 H5AD / 宽表 CSV 接入、单通道空间图、原始值 / log1p / asinh 显示、固定 ROI 的修订前后蛋白区域比较，以及配对 RNA 的原始 / 面板归一化 Spearman 对照。缺失值与实测零值分别保留。完整验证两片 SPOTS 的 5,336 个 spots、每片 21 个蛋白通道，真实矩阵逐值一致，离线包可复算。见 [蛋白组使用说明](docs/PROTEOMICS.md) 和 [v0.5 真实数据与验证](docs/V05_PROTEIN_REVIEW.md)。本版是同一样本的多模态协作基础，尚无跨样本整合、蛋白背景校正或自动联合注释。

本版针对真实使用中发现的四项不足，接入同坐标图像/分割证据，开放未知尺度下的表达与选区探索，保留跨样本对象身份和重复符号映射，并加入可编辑、可冻结的假设表与两种权重方法对照。已实际复测 Xenium 淋巴结 7,725 个细胞、HER2ST A1 346 个 spots、SPOTS 两片 5,336 个 spots，以及 Xenium 肾脏小样例的图像/分割对应。具体结果与剩余边界见 [四项修复验收](docs/V04_REAL_DATA_REVIEW.md)；原 [v0.3 真实数据报告](docs/REAL_DATA_STUDY.md) 保留为历史记录。

当前可用的闭环：

1. 导入 AnnData、Xenium H5/MTX、通用 10x MTX、标准 Visium spots、CosMx/MERSCOPE 指定 CSV 或 SpatialData Zarr 的表达表与 Points/Shapes。
2. 在网页或原生 napari 工作台查看质心、框选/圈选或输入确切对象 ID；Agent 用同一对象、版本和坐标系探索。明确区分 cell、spot、bin。
3. 查看所选细胞的标记基因计数，区分“未测量”和“测量为零”。
4. 提出重标注、标为待确认、排除/重新纳入或区域标签修订，预览后确认生成新版本。
5. 固定同一 ROI 与参数，重算细胞组成和半径邻域指标，展示前后差异、分母及判断是否跨越预设阈值。
6. 保留旧结果及其过期状态；恢复历史状态产生新版本；导出可校验和重放的审阅包。
7. 对冻结 ROI 与明确的、不重叠背景，比较修订前后的标签比例、marker 检出率、原始均值和每万面板计数均值，保留每个分母。
8. `inspect_quality` 检查来源 QC 与表达库；`query_observations` 按明确的属性条件分页得到对象 ID，缺失、null、零值分开；`run_sensitivity` 在同一版本独立计算多种标签/纳入假设与半径，保存可复算结果，**不创建注释版本或伪造研究者批准**。网页可读取这些已保存的假设结果。
9. 网页假设表和 MCP 共享草稿版本。研究者可编辑对象条件、修改内容、半径、图范围与背景；冻结后运行全部组合，按实际差值排序，展示按边加权与按来源对象等权的结果、分母和未知/失败原因。修改冻结方案会产生新草稿，不能改写旧记录。见 [假设计划](docs/HYPOTHESES.md)。
10. 通过本地显式登记或 napari 图层快照绑定图像/分割的来源、坐标变换与细胞标签映射；网页可检查对应证据，napari 可叠加图像、标签和质心。缺图或缺映射保持未知，任何质心偏差都不自动判为分割错误。

网页自动检查其他客户端的共享状态；Agent 的具体提案可在网页中审阅。MCP 新增分页筛选、轻量状态游标、提案读取与能力发现，避免反复把整个数据集塞进对话。napari 的图层、事件和后台 worker 适配见 [napari 使用指南](docs/NAPARI.md)。它们都连接同一个项目，原始输入不被覆盖。

## 本机立即试用

本机已有真实分析项目可直接打开（此前下载的原文件保留在原目录，项目不会重新生成）：

```powershell
cd C:\Spatial\spatial-collab
$env:PYTHONPATH = 'src'
python -m spatial_collab serve --project projects/v05-real-data/spots-spleen-1 --port 8773
```

打开 <http://127.0.0.1:8773/>。主页显示 2,568 个 spots、32,285 个 RNA 特征和 1 个含 21 通道的蛋白数据层。进入蛋白页选择 CD19 并读取信号；当前共享 ROI 是按源坐标中位数划分的 1,285 个 spots。可载入最近保存的全通道结果，也可指定 RNA 基因 Cd79a 重新比较。这些项目仅保存在本机，不包含在可分发 ZIP 中。已有 v0.4 RNA 项目仍可单独启动，见 [四项验收报告](docs/V04_REAL_DATA_REVIEW.md)。

创建独立的合成演示：

```powershell
cd C:\Spatial\spatial-collab
python -m pip install -e ".[import,dev]"
python -m spatial_collab demo projects\my-demo
python -m spatial_collab serve --project projects\my-demo --port 8765
```

打开 <http://127.0.0.1:8765/>。已有 `projects/demo` 演示时，可运行 `scripts/launch_demo.ps1`，或直接 `serve --project projects/demo`。不要对已有项目再次执行 `demo`；创建操作拒绝覆盖。

演示包含 480 个**合成细胞**，并特意构造部分标记与标签不一致的例子。它们用于软件和交互验证，不能成为肿瘤生物学结果或细胞注释准确性证据。

建议体验：先选择一个区域并检查 `CD3D, LYZ, EPCAM, UNMEASURED`；记录一次基线比较；把可疑细胞改为 `Uncertain` 或排除；确认后用相同区域重算；最后导出审阅包。若修改后来源或目标类型消失，结果会保留为“不可判定”，不会报告成阴性结论。

## 导入自己的数据

### 已注释 AnnData `.h5ad`

需要唯一的 `obs_names`、稳定且唯一的特征 ID、`obsm['spatial']` 的 N×2 坐标，以及显式指定的原始整数计数矩阵。标签可以指定列名，也可明确从 `Unannotated` 开始。特征 ID 默认来自 `var_names`，重复符号时使用通用 import 的 `feature_id_key` 指向原生 ID 列；原符号单独保存，不改名合并。坐标单位必须明确为 `micrometer`、`pixel`、`array_index` 或 `unknown`。未知尺度可选区、看表达和修订，物理半径分析仍要求微米框架；不会猜测换算。

```powershell
python -m spatial_collab import-h5ad C:\data\annotated.h5ad projects\study-a `
  --slice-id section-01 --coordinate-system xenium_native_xy --units micrometer `
  --label-key cell_type --counts-layer counts
```

若 `X` 确实是原始计数，可显式 `--counts-layer X`。存在 `slice_id`、`sample_id` 或 `library_id` 列时会检查混样；其他列名用 `--slice-key` 指定。未记录来源切片的文件需要导入者核对单切片条件；插件无法从坐标推断真实样本身份。AnnData 中未保留的基因按未测量处理。

### Xenium 输出

第一版读取 `cells.csv.gz`（或 `cells.csv`）、`cell_feature_matrix.h5`，以及研究者已有的注释 CSV：

```csv
cell_id,label
aaaaaaab-1,T cell
aaaaaaac-1,Myeloid
```

标签表必须与细胞表逐一对应；矩阵条码、特征数量和类型必须一致。只读取 `Gene Expression`，排除控制探针和蛋白特征。

```powershell
python -m spatial_collab import-xenium C:\data\xenium-outs C:\data\labels.csv projects\study-a `
  --slice-id section-01 --name "Xenium section 01"
```

Xenium 也支持已解压 `cell_feature_matrix/` MTX 目录。所有适配器仍只构建表达计数与质心的分析快照；不读取个体转录本，不验证分割准确性。SpatialData 适配器显式选择 table、Points/Shapes element 和目标坐标系，按 instance ID 关联并应用储存的二维变换；目标框架的微米尺度必须由导入者明确确认。

完整输入要求、支持的具体格式版本、标定和示例配置见 [多格式导入指南](docs/FORMATS.md)。目录探测只建议候选 reader，不自动猜测单位、样本或有歧义的格式：

```powershell
python -m spatial_collab formats
python -m spatial_collab probe C:\data\sample
python -m spatial_collab import visium C:\data\visium-outs projects\visium-a --options C:\data\visium-options.json
```

可选依赖：`pip install -e ".[import]"` 用于 H5AD/H5；`pip install -e ".[spatialdata]"` 用于 SpatialData；在已有 napari 环境中 `pip install -e ".[napari]"` 启用原生 dock。基础网页、MCP 和 CSV/MTX 工作流不要求安装 napari 或 SpatialData。当前 Visium reader 限标准 spots，尚未支持 Visium HD 原生输出；SpatialData 能导入声明为 bin 的二维点/形状，并不等于具有 HD 平台适配。

默认预算：100,000 个对象、20,000 个特征、2,000,000 个非零计数（稠密矩阵为元素数）、单文件 2 GiB。H5AD 可显式配置最高 100,000 特征/10,000,000 非零计数，预算随来源保存；真实 SPOTS 完整 RNA 使用 40,000/10,000,000，无截取基因。Xenium H5 可通过 `bounds` 显式选择连续窗口。其它超限输入需创建明确范围子集；不静默抽样。网页每个视野最多 10,000 个对象；大特征轴通过目录分页。紧凑审阅包只保存一次源矩阵加不可变修订，单成员限 256 MiB、全包 512 MiB，导出前检查。整张高密度组织的性能仍未验证。

## Claude Code

本地依赖安装后，在启动 Claude Code 的同一终端设置项目：

```powershell
$env:SPATIAL_COLLAB_PROJECT = 'C:\Spatial\spatial-collab\projects\demo'
claude --plugin-dir C:\Spatial\spatial-collab
```

插件包含 `.claude-plugin/plugin.json`、`.mcp.json` 和 `skills/spatial-review/SKILL.md`。它通过 stdio 访问同一项目。网页工作台需要另开 `serve` 命令；两端通过项目数据库共享状态。设置项目环境变量后要重新启动宿主服务，不能只在另一个终端临时赋值。

本地开发安装方式依据 [Claude Code 官方插件文档](https://code.claude.com/docs/en/plugins)。原生插件清单验证、MCP 协议调用与真实 Claude 对话中的工具选择是不同验收层；详见验证记录。

## ChatGPT / Codex

同一后端提供 MCP stdio、Streamable HTTP `/mcp`，以及标准 MCP Apps HTML 资源。`open_project` 关联交互视图；其余工具也可独立工作。本地 Codex 包使用 `.codex-plugin/plugin.json` 中的 `${PLUGIN_ROOT}` 配置；没有写入个人插件市场或修改宿主配置。

ChatGPT 开发连接需要 [官方支持的远程 HTTPS 端点或 Secure MCP Tunnel](https://developers.openai.com/plugins/deploy/connect-chatgpt)。对于当前本地原型，优先让 Secure MCP Tunnel 启动下面的 stdio 进程，并在该进程环境设置 `SPATIAL_COLLAB_PROJECT`：

```text
python C:\Spatial\spatial-collab\scripts\run_server.py
```

在 ChatGPT 中按官方指引添加该测试连接，然后调用 `open_project`。开发模式是否可用由账号与工作区策略决定。此交付**没有创建隧道、公开端点或市场发布**，没有记录真实 ChatGPT iframe 中的宿主验收。把 `localhost` URL 粘贴给云端 ChatGPT 不等于完成连接。

当前 HTTP 服务只绑定回环地址，并限制 Host/Origin。不要直接把它公开到互联网。正式远程部署还需要按研究环境实现身份、授权与传输控制；这部分没有作为已实现能力交付。

## 指标具体检验什么

检验的可执行问题是：“在给定半径和分析范围内，来源类型细胞的邻居中，目标类型的比例，是否超过其背景丰度至少 `min_effect`？”

- 邻居由二维质心欧氏距离 ≤ 半径定义，去除自身。
- `roi_induced`：节点只包含冻结 ROI 内仍纳入的细胞。
- `whole_slice`：来源细胞限于 ROI，候选邻居包含整张切片仍纳入的细胞。
- 指标是**来源邻接边加权**的目标比例减去非自身候选总体中的目标丰度；不是每个来源细胞等权。
- `stable/changed` 只表示该指标是否跨越预设阈值；指标本身的变化和分母变化同时报告。`indeterminate` 表示前后至少一侧证据不够。

该背景是简单丰度参照，**不是空间置换零分布**，没有控制组织结构、边界效应和空间自相关，也没有 p 值。空间邻近不能直接推出细胞通讯或因果机制。单切片、多个 ROI、很多细胞都不能替代独立患者重复。详见 [分析定义](docs/METHOD.md)。

新增 `run_region_comparison` 固定前景与背景对象，比较标签比例与标记特征；背景默认是同切片其余对象，也可指定另一个不重叠的已保存选区。面板归一化不是全转录组归一化，marker 对比不是差异表达检验。未测量基因与空分母保留不可判定；改变标签不会改变原始表达。所有 spot/bin 输出均声明其观测单位，不把区域标签解释为纯细胞身份。

## 从任务检验协作

[SpatialBench](https://benchmarks.bio/spatial) 的公开任务强调真实数据、具体产物、空间结构、背景选择及有限分辨率。这里借鉴其任务思想，加入 **8 个开发者编写的合成协作任务**，覆盖测量范围、精确选区、变更敏感性、邻域背景和可复算性。公开题目的数值评分不等于验证推理方法，因此测试还检查对象、操作轨迹和证据边界。

```powershell
python -m spatial_collab benchmark projects\benchmarks
```

每次生成独立目录，保留输入、答案、操作、失败、结果与校验清单。这是脚本驱动的软件验收，**不是 SpatialBench 成绩，也没有测量 Agent 准确率或研究者节省时间**。官方资料、具体启发和下一阶段人机对照协议见 [任务与评测设计](docs/BENCHMARK.md)。

## 结果、重放与测试

项目的 `project.sqlite3` 保存来源快照、选区、修订、提案与分析结果。导出物位于 `exports/review_*/`，包括 SHA256 清单、前后细胞表、参数、修订历史与计算环境。默认这些原始记录留在本地；通过 MCP 返回的摘要、坐标、标签及对象引用可能发送到宿主模型。

```powershell
python -m spatial_collab inspect --project projects\demo
python -m spatial_collab export --project projects\demo --run-id <run_id>
python -m spatial_collab verify-bundle <export_directory>
python -m pytest -q -p no:cacheprovider
python -m ruff check --no-cache src tests scripts
```

姓名是修订归属记录，不是身份认证、专家背书或科学授权；文件散列也不是独立科学验证。SQLite 的事务、比较并交换和内容校验用于避免状态错配，不能抵抗拥有完整本机写权限的恶意修改者。

完整交付范围、已验证项目和未完成验收见 [验证记录](docs/VALIDATION.md)；后续产品实验见 [下一轮真实任务验证](docs/PILOT.md)。
