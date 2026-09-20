# 描述性注释敏感性：v1 方法约定

## 分析对象

一个项目只包含一张切片和一个明确的二维微米坐标系，来源细胞 ID、坐标和计数不可变。版本覆盖层只改变 `label`、`included` 和 `region`。ROI 保存几何、坐标系、切片、来源散列、版本及具体细胞集合。多边形按质心包含关系选取，边界上的点计入；并不声称包含完整细胞形状。

比较固定两个版本及同一来源细胞集合。排除只改变其纳入状态，不能改写这个集合。区域标签修改不会重选 ROI；若科学问题要求改变区域边界，需要另外捕获选区。本版本未实现两套不同 ROI 的匹配比较。

## 邻域指标

令候选节点集合为 V，选区中来源标签的有效细胞为 S，目标标签为 T。半径 r 下记录有向邻边 `(s,v)`，其中 s∈S、v∈V、s≠v 且距离不超过 r。每个细胞 ID 单独计数，同坐标不同细胞仍是邻居。

目标邻边数 E_T 除以总来源邻边数 E，得到 `observed_fraction`。背景比例 B 为候选节点中目标细胞数量除以 |V|−1；当来源与目标标签相同，分子也减一。这保证排除自身后分子分母一致。

`excess_over_abundance = E_T/E − B`。给定 0 < τ ≤ 1，当该值 ≥ τ 时标为 `descriptive_supported`，否则 `descriptive_not_supported`。阈值应在分析前由研究者声明；调参后更好的结果不构成独立确认。

两侧均可计算时，是否跨越 τ 决定 `stable` 或 `changed`，同时报告效应变化、原始比例变化、背景变化、图节点数、来源细胞数及边分母变化。来源、目标、边或背景不可用时为 `inconclusive`，前后比较为 `indeterminate`。不把群体消失当作阴性证据。

这是描述性参照，尚未实现空间置换检验、病人层面推断、批次控制、组织隔离边界或生物学机制检验。

## 重算、失效和恢复

每次比较重新构建有效候选总体和半径邻域，采用 SciPy cKDTree。先查询邻居数量，再按单个来源细胞读取邻居；不建立 N×N 距离矩阵。预算限制来源有向边数和单个来源邻居数，超限给出失败，不能自动改变半径或抽样。

同一提案基于固定 head revision；提交使用事务中的比较并交换。两个宿主对同一版本并发修改时，最多一个提交成功，另一方必须重新预览。没有“最后写入者悄悄覆盖”路径。

结果保持不可变，`stale` 反映目标版本是否仍是当前版本。历史比较仍可解释为当时的结果，但不能代替当前版本。所有已有项目内结果保守处理；外部分析文件不在这个失效跟踪图中。

恢复历史覆盖层生成新 revision。若与历史数据状态相同也不能冒用旧 revision 的身份。导出记录当前状态与运行状态，重放不依赖聊天记录。

## 实现取舍

当前 working snapshot 使用 SQLite 和标准 JSON 审阅导出，仅索引本任务所需细胞与计数。原始 AnnData/Xenium 文件和散列保留为来源。没有替代 SpatialData 的多模态数据标准，也没有实现 SpatialData Zarr 直接读取。

为了先验收共享对象与修订闭环，界面采用有界质心画布。真实图像、分割与分块转录本应通过后续经过坐标验证的 Vitessce/SpatialData 适配补入，不能据质心视图声称已实现组织形态复核。

## 官方接口参考

- [10x Xenium 输出定义](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/outputs/xoa-output-understanding-outputs)：细胞表微米质心、表达矩阵及分割文件语义。
- [AnnData 文件读取](https://anndata.readthedocs.io/en/stable/generated/anndata.io.read_h5ad.html)：AnnData 输入格式。实现只读所需字段，不修改输入。
- [官方 MCP Python SDK](https://github.com/modelcontextprotocol/python-sdk)：传输与工具 schema。
- [MCP Apps 标准](https://modelcontextprotocol.io/extensions/apps/overview)、[OpenAI UI 接口](https://developers.openai.com/plugins/build/chatgpt-ui)：共享 HTML 资源和宿主消息。
# v0.2 扩展：观测单位与区域背景

兼容接口保留 `cell_id`、`cells` 字段名；科学含义由 `observation_unit`（cell/spot/bin）与 `label_semantics` 明确声明。spot/bin 的标签为区域工作注释，画成点或运行半径图均不能使其变为单细胞，也不能视为去卷积。旧项目不重写源快照；展示时以原版 cell 语义兼容。

`run_region_comparison` 保存前景 selection、背景 selection（或同切片补集）和两个修订。两个集合必须互不重叠，跨版本不改变成员。每个版本分别移除被排除对象，然后计算标签比例差、marker 原始均值差、检出率差及每万面板计数均值差。标签分母为该组纳入对象数；原始均值和检出率包含零库对象；归一化均值仅包含面板总计数大于零的对象，单独报告其数量。未测量基因不补零，空分母不输出支持性结论。

标签编辑不改变表达，排除操作会改变组内分母。该分析没有统计检验、独立重复建模、批次校正或差异表达宣称；面板归一化仅对应导入的特征范围。背景为同一切片不自动证明其为生物学匹配对照，仍需研究者确认取材和组织条件。

两种结果格式均可导出后离线重算；v0.1 的原导出包仍保留原格式和散列。状态游标只表示当前项目状态，不能用于声称具有完整事件流；没有异步长任务调度、取消或多步作业恢复系统。

