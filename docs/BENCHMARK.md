# 从 SpatialBench 到研究者与 Agent 的协作任务

核查日期：2026-09-20。这里记录公开来源、采用的设计思想、当前本地验收的边界。Spatial Collab 与 SpatialBench 没有隶属或认证关系。

## 公开基准实际做了什么

SpatialBench 将真实流程切在某个分析步骤之前：交付当时可用的数据对象、科研问题、结构化答案要求和确定性评分器。它要求结论在合理方法选择下稳定，并通过人工检查和反捷径尝试剔除仅凭常识就能答对的题目。评分器检查最终答案，不能替代对过程和证据的审阅。[论文 v2，Methods 3.1–3.4](https://arxiv.org/html/2512.21907v2)

版本必须分开：论文 v2 描述 146 题；当前官方仓库说明 159 题，五个平台，七类任务；2026-09-20 实际打开的网页显示 Full 159 和 Verified 115 两种集合。五个平台为 Seeker、Visium、Xenium、AtlasXomics 和 MERFISH；七类为 QC、归一化、降维、聚类、细胞注释、差异分析和空间分析。仓库仅公开代表性题目，完整集合保留以减少污染。[官方仓库](https://github.com/latchbio/spatialbench)、[实时网页](https://benchmarks.bio/spatial)

当前仓库的方法说明记录每题运行三次、4 vCPU、32 GB RAM、100 GB 存储和六小时时限；一般不限制调用步数，环境允许网络，并记录内存耗尽重启。它与旧论文中的运行设置不能混写。评估必须锁定模型、宿主、工具版本、题目版本、预算及重试规则。[当前 METHODS](https://github.com/latchbio/spatialbench/blob/main/METHODS.md)

SpatialBench-Long 的公开仓库包含四个示例及提示、词表、轨迹和结果；网页上的全集数量与公开样例数量不是同一回事。长流程提醒我们检查前面决策对后续步骤的累积影响。[官方 Long 仓库](https://github.com/latchbio/spatialbench-long)

## 具体科研题目给本产品的启发

| 公开任务 | 必须理解的科研区别 | Spatial Collab 的协作要求 |
| --- | --- | --- |
| Seeker 卵巢卵泡计数 | 同一种表达状态的细胞可以分属不同空间结构；多个 beads 可以采到同一卵母细胞 | 同时保留观测单位、标记表达、坐标和结构成员集合；避免把标签计数称为结构计数 |
| MERFISH 内皮细胞双细胞评估 | 单个跨谱系标记共表达不能独自区分真实表达、混合或分割错误 | 研究者能够保留分歧、要求多标记和空间/形态补充证据；不把共表达自动转成排除 |
| 稀疏 Visium 骨组织谱系分辨率 | 参考图谱中存在某细胞类型，不代表本次测量能分辨它 | 同时显示面板覆盖、检测率和观测单位；允许粗粒度工作标签与未知标签 |
| Xenium CN7 构成与富集 | 最丰富的细胞类型可以不是相对背景最富集者；背景须来自同一条件 | 在结果旁保存分子、分母、条件和参考集合；修订后复算这些对象 |
| DBiT-seq 供者间性别差异 | 条码数不等于生物学重复数 | 样本/供者标识成为设计约束；缺乏独立重复时保留描述性结果，不能制造群体推断 |

上述案例来自实际公开题目：[卵泡](https://github.com/latchbio/spatialbench/blob/main/example_evals/cell_typing/curio_ovary_follicle_count_immature.json)、[MERFISH](https://github.com/latchbio/spatialbench/blob/main/example_evals/qc/qc_01_endothelial_doublet_assessment.json)、[Visium](https://github.com/latchbio/spatialbench/blob/main/example_evals/cell_typing/visium_bone_celltype_lineage_resolution_limit.json)、[CN7](https://github.com/latchbio/spatialbench/blob/main/example_evals/spatial_analysis/xenium_kidney_spatial_cn7_composition_day14.json)、[供者统计](https://github.com/latchbio/spatialbench/blob/main/example_evals/differential_expression/MOTIF04_neuronal_female_tf_family.json)。最后一题属于空间 ATAC，借鉴的是统计单位原则，并不表示本插件支持 ATAC 推断。

必须保留批判性判断：公开卵泡题实际评分仅检查计数是否位于容许区间，未验证是否使用空间信息；正确数字可能来自错误过程。MERFISH 题的目标值由题目作者的启发式证据推导，不能推广为所有组织的分割真值。Visium 题的检测阈值也不是跨平台通用判据。我们采用问题结构与失败模式，不复制答案、阈值或科学确定性。

## 本地可执行验收与外部协作研究分开

`spatial_collab.benchmark` 使用本项目自行构造的小型合成数据，生成科研问题、固定初始对象、运行记录、结构化观察值和确定性评分。它运行的是预先写好的参考操作流程，**不是对 ChatGPT、Claude、研究者效率或 SpatialBench 的评分**。它不下载原基准数据，也不调用外部模型或提交任何结果。

当前验收关注用户真正需要的行为：

1. 圈选必须返回精确 ID，包含坐标系和版本，而不是根据截图重猜。
2. 面板未测量与测量为零必须区分；无计数时仍能继续有限的几何/标签探索。
3. 对研究者指定的标签修改，固定比较对象和参数，检验指标是否改变；不要求每次修订都改变结果。
4. 缺少目标群体应为不可判定，不能偷换成“没有生物学效应”。
5. ROI 内部邻域与保留 ROI 外邻居的分析属于不同问题；背景集合必须可见。
6. 重开项目和重算导出包后，应恢复同一版本与数值；历史结果不能冒充当前结果。
7. 缺失空间单位应阻止物理距离解释，同时指出补充什么信息后可继续。
8. 同一份观测表声明为 spot 后，保存并恢复 spot 语义；不会因为兼容 API 沿用 `cell_id` 就当作单细胞数据。

本地结果单独报告对象一致性、科学数值检查、续接/重算状态和失败原因。任何失败保留原始记录；异常不计通过。合成输入和人工指定答案公开，适合回归开发，不能称为独立或留出测试。

开发者可调用 `run_benchmark(output_dir)`。每次在目标目录下创建新的运行子目录，包含 `tasks.json`、逐题输入、实际回答、操作轨迹、评分细项、`report.json` 和文件校验清单。已有项目不会被载入或覆盖。输入缺失、错误数值、重复 ID、布尔值冒充计数、不可判定结果缺少显式空值都不能获得通过。这个验收目录可离线保存和检查，不包含真实受试者数据。

## 后续真实协作评估协议

固定真实数据切片、任务起点和专家审核后的参考范围，按组织与研究来源划分开发和留出集合，不能只按同一切片中的细胞随机划分。独立维护答案与模型可见工作区。每个任务先由至少两种合理分析路线验证数值稳定性；有争议的生物学判断保留多种可接受结果和专家分歧。

对照应包括研究者原有流程、同模型的基础工具流程、同模型加插件流程。跨 ChatGPT、Claude Code、Codex、napari 比较时，模型版本、工具权限、运行预算和人工帮助都作为独立变量记录。可采用随机交叉顺序以降低熟悉题目造成的偏差。

每个任务保存：输入与任务版本、观测单位、空间单位、供者/切片/条件、数据缺失与资源限制、模型与宿主版本、人工修订次数、恢复点、工具调用与外部动作、耗时、费用、最终结构化结果、结果包及核对记录。只记录可观察的操作与简短理由，不索取隐藏思维链。

评估首先检查科学问题是否完成，再检查对象和证据是否正确，最后报告效率。分别汇总正确完成率、合理不可判定率、错误结论率、首次有效结果时间、研究者主动操作时间、修复时间、断点续接成功率；不能用少调用几次工具替代科研效率。按任务、平台、组织和样本独立性分层，展示失败轨迹及置信区间。

外部检索或上传不是默认为任务成功所必需的步骤。需要时记录具体目的、对象与用户授权；本地可完成的任务保留零外部动作对照。宿主轨迹必须实际采集，不能用模型自述“我未上传”作为证据。

## 与 napari 架构的对应

数据、状态与显示分离：不可变导入对象和版本是模型，Web/napari 是视图，MCP/用户界面操作是控制入口。格式 reader 只负责明确格式的解析；图层保留观测 ID、单位和版本引用。视觉点的行号不能替代稳定观测 ID。耗时任务不应阻塞观察，完成后必须核对输入版本再展示结果。napari 的图层、事件模型和 reader 贡献机制提供成熟参考；具体桥接实现及验收状态以平台文档为准。[napari 模型与事件](https://napari.org/stable/developers/architecture/napari_models.html)、[reader 插件](https://napari.org/stable/plugins/building_a_plugin/guides.html)
