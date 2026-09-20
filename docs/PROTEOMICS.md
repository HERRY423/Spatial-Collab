# 空间蛋白组 · v0.5

同一样本的 RNA 与蛋白共用空间对象、选区、注释版本和纳入状态。RNA 源计数保持原样；蛋白存为独立、带来源散列的不可变数据层。这里的“整合”指对象配对与协作，不是联合嵌入、批次校正或自动细胞分类。

## 工作台

- **主页**：当前样本、RNA 特征数、蛋白数据层、共享 ROI、近期保存的分析。可以进入相应组学、检查 RNA 质量和导出记录。
- **空间转录组**：保留圈选、标记证据、修订预览/确认、图像对应、冻结假设和邻域分析。
- **空间蛋白组**：选择数据层和蛋白通道，读取全部导入对象的空间信号；选区白色描边。右侧统计是当前共享 ROI，没有 ROI 时统计全部对象。灰色表示缺失或排除，实测零值参与统计。可以点击一个对象建立精确选区；复杂圈选在 RNA 页完成。

切换页面不改变共享选区、不提交修订、不丢弃假设草稿。其他客户端更新选区或版本后，旧蛋白图/统计清空并提示重新读取；已保存分析仍绑定原来的参数和 ROI。网页当前最多显示 10,000 个对象，超出时明确省略图，完整范围的统计仍可用。

区域比较以当前共享 ROI 为前景，其余导入对象为默认背景，也可指定非空、不重叠的历史选区。基线和目标版本可分别选择。一次比较 1–32 个通道，每个通道保留均值、测量/缺失分母、零值、标签组成和前后区域差。修订标签未必改变全标签蛋白均值；排除对象会改变分母。只有研究者确认的原工作流才提交注释。

可选 RNA 基因提供同一前景内的配对 Spearman 相关，同时保留原始 RNA 与每万面板计数 RNA 两种视图。配对少于 3、RNA 或蛋白为常数时为 unknown；RNA 零库不进入归一化对照，分母完整记录。没有 p 值、细胞比例估计或患者层面的统计推断。

## 接入测量

当前支持配对 H5AD 和宽表 CSV。接入操作在本机显式执行，尚无浏览器文件上传器。H5AD 需要原始 obs ID、`obsm['spatial']` 和明确的 X 或 layer。CSV 格式如下，空值是缺失：

```csv
cell_id,x,y,CD3,CD19
original-barcode-1,5,7,0,12
original-barcode-2,6,8,,4
```

CSV 除 ID/x/y 外的列全部视为蛋白通道，勿混入元数据。H5AD 可指定 `feature_id_key`、`feature_symbol_key`；重复 symbol 可保留，稳定 ID 必须唯一，禁止静默聚合。默认使用原 var 名称，未猜测 antibody target 或同义词。

登记配置示例：将以下 JSON 保存为 `protein-options.json`，填写当前项目实际样本、源散列和坐标系（可从 `get_overview` / `open_project` 获取）。

```json
{
  "sample_id": "your-recorded-sample-id",
  "source_sha256": "the-current-project-source-sha256",
  "name": "Paired antibody panel",
  "measurement_type": "antibody_count",
  "coordinate_system": "your-recorded-coordinate-system",
  "units": "array_index",
  "registration_note": "Record the evidence that this is the same specimen and spatial coordinate frame",
  "format_id": "h5ad",
  "matrix": "X",
  "spatial_key": "spatial",
  "allow_partial": false
}
```

```powershell
python -m spatial_collab register-protein paired-adt.h5ad --project projects/my-study --options protein-options.json
python -m spatial_collab serve --project projects/my-study --port 8773
```

`antibody_count` 只接受非负整数或缺失；`intensity` 接受非负有限小数或缺失。负值校正矩阵暂不支持。H5AD 含 sample_id 时必须匹配。原始 ID 和 XY 必须逐对象精确对应，禁止按行序、最近邻或条码交集偷偷配对。全部对象默认必须匹配；`allow_partial: true` 可显式保留缺失覆盖，外来对象仍拒绝。未知尺度允许探索，不能声明 µm 半径。

限制：每项目最多 16 个蛋白层，每层最多 100,000 个对象、512 个通道、5,000,000 个矩阵元素；源文件 256 MiB、序列化层 32 MiB。适配任意平台前需提供对象和坐标对应，当前未对 CODEX、IMC、MIBI 原生文件做真实数据验收。

原始值、log1p 和 asinh(value/cofactor) 是显示/汇总尺度。它们不构成背景扣除、抗体特异性验证或批次校正；asinh 的系数完整保存。不同测量类型不合并为一种表达值。

## Agent 与可复算记录

MCP 增加 `get_overview`、`list_assays`、`inspect_protein`、`run_protein_comparison`，总共 31 个工具。现有共享选区/修订工具通用。网页将当前组学及已读蛋白数据层、特征、版本、选区写入兼容宿主的模型上下文。

`export_review_bundle` 保存蛋白矩阵快照、来源信息和散列，蛋白运行绑定确切 assay 散列。离线 `verify-bundle` 会重建同一注释/纳入状态及 ROI 并重新计算；旧 RNA-only 包保持兼容。完整性校验不是来源真实性或科学结论的第三方认证。

本版蛋白功能通过本地网页和 MCP 工具使用。napari 继续共用 RNA 选区、修订、图像/分割；尚未新增原生蛋白 dock。ChatGPT 的远端接入仍需用户部署和鉴权；本地协议通过不等于实际 ChatGPT / Claude 模型对话验收。
