# 大数据、分辨率与多样本推断改造

版本保持 **0.7.0a1 / 0.7.0-alpha.1**，本次不增加版本号。新增“全切片与研究统计”页面、对应 MCP 工具及本地导入命令。

## 数据路径发生了什么变化

原有 JSON 项目继续承担有版本的注释与有界分析。大数据使用新增的磁盘 atlas：按行/块读取数据，将观测、稀疏计数、分子和多边形写入独立 SQLite 数据库，建立空间 R-tree、特征索引和密度金字塔。不会把百万对象的全部计数转换成 Python 字典或发送给浏览器。多个 atlas 可以附在同一工作项目中，但每个数据集保留独立样本、坐标系、来源和对象身份，不冒充主项目的同一张切片。

全景按全部对象计算密度；放大后返回精确对象。密度格明确是计数汇总，没有伪造的细胞 ID。10,000 是一次精确点响应的预算，超过后自动返回完整密度层，不再出现整张图空白。框选视野可冻结精确 ROI 的稀疏计数，进入现有可执行科学分析。算法自己的内存预算继续生效，未宣称所有算法都能一次拟合百万对象。

原始转录本 ID、特征、x/y、可选 z、质量值和源分配细胞 ID 均保留。显示为明确的二维投影，不把低质量分子自动删除或把分配关系当作真实细胞归属。细胞/细胞核边界按确切原始 ID、顶点次序及类型登记，和细胞质心分开存储。全局范围覆盖分子和边界，避免只用质心范围裁掉边缘分子。

图像直接读取原生 tiled TIFF / OME-TIFF / BigTIFF 金字塔。每次只解码与一个 256 像素显示瓦片相交的源块，支持 YX、交错 RGB/RGBA，以及显式指定 T/Z/C 等额外轴的平面。使用全分辨率像素中心到世界坐标的 3×3 仿射；视图保持长宽比，底图、分子和边界同步缩放。非 uint8 图像要求声明统一对比度，不逐瓦片自动拉伸而产生接缝。

## 安装与本地导入

```powershell
.venv-v07\Scripts\python -m pip install ".[large-data,study]"
.venv-v07\Scripts\python -m spatial_collab large-data init --project projects/my-study --name "My spatial study"
.venv-v07\Scripts\python -m spatial_collab large-data import --project projects/my-study --spec atlas.json
.venv-v07\Scripts\python -m spatial_collab large-data image --project projects/my-study --spec image.json
.venv-v07\Scripts\python -m spatial_collab large-data verify --project projects/my-study --atlas-id atlas_ID
```

本机已验证的完整环境保存在 `requirements/scalability-py313-win64.lock.txt`。安装命令从源码安装的包版本仍是 0.7.0a1。`init` 创建没有虚构细胞、特征或坐标的空容器，浏览器直接进入全切片页面；已有工作项目则跳过此步。大型 atlas 的数据身份独立于项目原有数据。

`atlas.json` 的最小示例（路径与列名必须替换为生产者实际字段）：

```json
{
  "metadata": {"name":"Study slide 1", "sample_id":"S1", "species":"human", "platform":"Xenium", "observation_unit":"cell", "coordinate_system":"producer_xy", "units":"micrometer"},
  "cells": {"path":"C:/data/cells.parquet", "columns":{"id":"cell_id", "x":"x_centroid", "y":"y_centroid"}},
  "features": {"path":"C:/data/features.csv", "columns":{"id":"feature_id", "symbol":"symbol"}},
  "counts": {"path":"C:/data/cell_feature_matrix.h5", "format":"10x_h5", "feature_field":"id"},
  "transcripts": {"path":"C:/data/transcripts.parquet", "projection":"xy_projection", "columns":{"id":"transcript_id", "feature":"feature_name", "x":"x_location", "y":"y_location", "z":"z_location", "qv":"qv", "cell_id":"cell_id"}},
  "boundaries": {"path":"C:/data/boundaries.csv", "columns":{"id":"boundary_id", "cell_id":"cell_id", "kind":"kind", "vertex":"vertex_index", "x":"x", "y":"y"}}
}
```

- cells/features/counts 是必需输入；完整特征面板包括测量为零的特征。
- CSV、CSV.gz 与 Parquet 均支持流式读取。counts 也支持长表 `cell_id,feature,value`，或 `format: h5ad_csr` 加显式 `counts_layer` 和 `feature_field`。不静默把归一化 X 当作原始计数。
- 矩阵的观测和特征轴必须与声明表完全一致，包括零库对象与全零特征；重复 ID、重复 cell-feature 记录、外键错误会中止导入。
- 边界表需要显式 boundary_id、cell_id、kind（cell/nucleus）和顶点次序。原生平台若没有该列，需按生产者文件契约准备，不猜测顶点顺序。可通过 `constants: {"kind":"cell"}` 为单类型边界文件显式指定类型。
- 分子特征名不一定和计数矩阵的稳定基因 ID 相同；保留各自生产者标识，不猜测合并。图层筛选使用转录本表中实际的 feature 值。
- 当前没有专用 Stereo-seq GEM、CosMx 专有容器、Visium HD 全目录的一键识别器；可用声明后的 CSV/Parquet 或原始稀疏矩阵导入。未将通用格式支持冒称逐平台实测。

`image.json` 示例：

```json
{"path":"C:/data/morphology.ome.tif", "atlas_id":"atlas_ID", "series":0, "selectors":{"C":0}, "contrast":[0,5000], "pixel_to_world":[[0.2125,0,0],[0,0.2125,0],[0,0,1]]}
```

上述像素大小只是格式示例，必须来自该图像实际标定。axes 为 YX 时 selectors 使用 `{}`。不同通道可分别登记后选择。普通 strip TIFF 被明确拒绝，应先在图像工具中转换为 tiled pyramid；不会为显示而读取整张巨图。当前支持原生 TIFF 金字塔，不声称已实现远程 OME-Zarr、DZI 或云对象存储。

## 多样本统计现在可以直接执行

登记 Study → Subject → Sample → Section、条件/批次、ROI 来源、run 与 source 哈希。每个切片提供定义一致的连续结果，例如已明确变换的通路分数或组织指标。科学输入与其生成方法版本、单位绑定；直接使用原始 RNA counts、未变换比例、spot 数充当受试者数均不在模型合同内。

| 方法 | 统计单位与输出 | 约束 |
|---|---|---|
| Welch | 每个受试者内先平均切片；报告组均值差、t、自由度、95% CI、p/q | 独立受试者、单批次；两条件共用受试者时拒绝 |
| paired_t | 同一受试者两个条件的切片均值配对；报告均值差、CI、p/q | 两条件受试者完全对应；不默默删除未配对对象；不将 pair_id 自动当成同一人 |
| mixedlm | 切片连续值 ~ 条件 + 批次 + 可选数值协变量，受试者随机截距，REML | 明确是 subject 随机截距模型；不支持任意随机斜率、样本嵌套随机效应或空间残差结构 |

完全混杂或固定效应秩不足会拒绝推断。混合模型不收敛、奇异或随机效应落在边界时保留诊断、不给 p 值。MixedLM 的区间与 p 值为渐近正态 Wald，未实现 Satterthwaite/Kenward–Roger 小样本校正；少于 20 位受试者明确提示。最低样本数检查不是功效证明。

每次分析冻结完整指标家族，失败指标以 p=1 保留在 BH 分母中；不因为某项拟合失败而缩小家族。缺失默认不检验该指标；若研究者明确选择 complete_subjects，排除该指标中有缺失的整位受试者并列明名单。事后选择的 ROI 保持探索性标记。全部来源绑定仍属于研究者声明，并不独立核实实验身份。

界面可上传研究设计与逐切片指标 JSON，选择方法、比较条件、ROI 类别、协变量和缺失策略，直接查看效应、置信区间、FDR 和诊断。MCP 对应 `register_study_design`、`run_study_inference`、`list_study_analyses`。最小指标记录：

```json
{"section_id":"section-1", "source_sha256":"与设计相同的来源哈希", "run_id":"设计中登记的run", "metric":"pathway_score", "value":1.2, "units":"score", "method_version":"明确的生成方法版本"}
```

研究设计的每个 sections 条目需要 subject_id、sample_id、section_id、condition、batch、sample_source、roi_definition、roi_class、roi_origin、source_sha256、assay_ids、run_ids。数值协变量置于 `covariates` 对象。完整、可直接运行的**合成**控制模板由 `scripts/validate_study_models.py` 写入 `output/scalability/study-tmp/`；不能把这些文件当作真实参与者数据。

## 导出与证据边界

现有审阅包保存 atlas / pyramid 元数据和哈希，以及冻结 ROI 的实际稀疏计数、研究设计和统计结果。巨型数据库和原始图像不塞进原有 JSON 包；离线验证会明确报告外部资源仅有清单，迁移时应连同 atlas 数据库与图像另行归档。`large-data verify` 对磁盘数据库执行完整哈希校验；普通交互查询检查文件大小/时间绑定，不能把这当作抵御蓄意保留文件时间的篡改检验。

Welch、配对与 LMM 结果可从包内设计和指标复算。未招募独立用户、未制造多患者公开队列，也未将合成规模或统计控制当成生物学有效性验证。实际运行结果见 [验收记录](SCALABILITY_VALIDATION.md)。

实现依据：[tifffile 官方文档](https://www.cgohlke.com/docs/tifffile/)、[statsmodels MixedLM](https://www.statsmodels.org/stable/generated/statsmodels.regression.mixed_linear_model.MixedLM.html)、[SciPy Welch 检验](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html)。
