# 空间转录组导入格式与明确边界

本文描述 v0.5 保留的 RNA 读取器。导入创建新的工作项目，保存原始文件摘要、原始计数和质心坐标；后续注释修订只写入项目中的版本覆盖层。它不改写原始 AnnData、Xenium、CSV 或 Zarr 文件。新增配对蛋白 H5AD / 宽表 CSV 作为独立测量层登记，见 [蛋白接入说明](PROTEOMICS.md)。

一个项目只包含一张切片和一个明确声明单位的坐标系。H5AD 可保留未标定坐标：`array_index`、`pixel`、`unknown` 仍可探索表达与修订注释，但不能进行微米半径分析。`cell_id` 是兼容性的内部观察对象 ID；观察对象可以是 `cell`、`spot` 或 `bin`。点位、bin 的标签是区域或工作注释，不能解释为纯细胞身份或去卷积结果。

## 发现能力、探测文件、明确选择读取器

```powershell
python -m spatial_collab formats
python -m spatial_collab probe C:/data/example
python -m spatial_collab import tenx_mtx C:/data/example/matrix C:/Spatial/projects/example --options C:/data/example/options.json
```

`formats` 返回声明式能力清单，不导入 AnnData、SpatialData、HDF5、NumPy 或 SciPy。`probe` 只检查文件名及有界 CSV 表头，不加载计数矩阵；候选可能不止一个，它不会自动选择。Zarr 标记只产生候选，不能证明容器符合 SpatialData 模型。

`--options` 必须是 JSON 对象，文件最多 64 KiB。附属文件参数建议写绝对路径：相对路径按命令执行时的工作目录解析，不按 options.json 所在目录解析。导入项目的路径必须是新的项目；现有项目不会被覆盖。

下面 JSON 都是可修改的模板。所有 `0.5`、`0.12` 比例只是**假设的示例校准值**，不是对你的切片、仪器或图像的事实判断。必须替换为来自该数据集的可信校准。不能为了让程序运行而宣称像素坐标已经是微米。

## 依赖和共同限制

在项目根目录安装对应依赖：

```powershell
python -m pip install -e .
python -m pip install -e ".[import]"
python -m pip install -e ".[spatialdata]"
```

| 读取器 | 路径类型 | 额外依赖 | 已实现范围 |
| --- | --- | --- | --- |
| `anndata_h5ad` | `.h5ad` 文件 | `.[import]` | 指定原始计数层、N×2 空间坐标、明确的观察单位 |
| `xenium` | Xenium 输出目录 | H5 路径需要 `.[import]`；MTX 用基础依赖 | `cells.csv[.gz]` 加 H5 或已解包的 MTX 矩阵 |
| `tenx_mtx` | MTX 目录 | 基础依赖 | 10x feature×observation 矩阵、外部质心和注释 CSV |
| `visium` | 标准 Visium `outs` 目录 | 基础依赖 | 传统标准 spot，MTX 矩阵，旧或新 CSV 位置文件 |
| `cosmx_csv` | flat-file 目录 | 基础依赖 | RNA 宽表，原生全局像素质心，明确的 RNA 列 |
| `merscope_csv` | MERSCOPE CSV 目录 | 基础依赖 | 全局微米质心、cell-by-gene 原始计数 |
| `spatialdata_zarr` | SpatialData Zarr 目录 | `.[spatialdata]` | 精确 table→Points/Shapes 关联及目标坐标变换 |

共同的默认 alpha 上限是 100,000 个观察对象、20,000 个特征、2,000,000 个存储非零条目。H5AD 可显式声明 `max_features`、`max_nnz`，最高分别 100,000 和 10,000,000；预算随来源保存，不会静默放宽或裁减基因。稠密矩阵的同一条目预算按全部元素计。单源文件最多 2 GiB。CSV 计数读取还限制为 20,000,000 个被检查的值；新增文本读取器限制单行长度及解压后文本大小。SpatialData 另有 2 GiB/10,000 文件、数组元素及解码后元数据预算。仍超限时需要明确的子集或其他存储方案，不会悄悄采样。

仅接受有限、非负、可精确表示的原始整数计数。归一化、对数变换值、重复观察 ID、重复稳定特征 ID、错位条形码，以及跨切片混合数据都拒绝导入。重复基因符号必须用不同的稳定特征 ID 保留，不能自动加后缀或合并计数。0 是已测量的零；未出现在本次导入特征范围内的基因是未测量，不会补成零。

所有读取器默认要求注释。明确设置 `allow_unannotated: true` 后，可先以 `Unannotated` 占位进行数据探索；该占位不构成细胞类型、空间邻域组别或生物学结论。传入注释 CSV 时，ID 集合仍必须与实际导入集合完全一致，不能用此开关允许错误或部分 ID 连接。

`biological_replicates` 默认为 0，表示未知；单切片模式只接受 0 或 1。观察对象数量、FOV 数量和 ROI 数量都不能作为生物学重复数。

## AnnData：`anndata_h5ad`

```powershell
python -m spatial_collab import anndata_h5ad C:/data/example/data.h5ad C:/Spatial/projects/anndata-review --options C:/data/example/anndata-options.json
```

```json
{
  "slice_id": "slice-A",
  "sample_id": "sample-A",
  "coordinate_system": "native_xy_um",
  "units": "micrometer",
  "platform": "declared assay",
  "observation_unit": "cell",
  "label_key": "cell_type",
  "spatial_key": "spatial",
  "counts_layer": "counts",
  "allow_unannotated": false
}
```

读取 `obs.index` 作为唯一原始观察 ID、`obsm[spatial_key]` 作为 N×2 的 x/y 坐标。显式提供 `sample_id` 后，内部 ID 使用 `sample:<编码样本>::cell:<编码原始ID>`，顶层 `sample_id`、`source_cell_id` 保存原值；metadata 的 `identity_scope` 为 `sample_qualified`。这样不同切片的同名 barcode 不会成为同一对象。特殊分隔符和 Unicode 使用百分号编码，避免拼接碰撞。未提供 `sample_id` 的旧调用保留项目内原 ID，不伪造样本名，也不宣称跨项目 ID 唯一。若原始 `obs['sample_id']` 已有值，显式参数必须与之相同。

`feature_id_key` 指定 `var` 中的唯一稳定特征 ID 列，例如 `gene_ids`；省略时使用 `var.index` 并要求唯一。`feature_symbol_key` 指定符号列；省略时保留 `var.index` 的原始符号。`metadata.features` 保存完整 `{feature_id,symbol}` 映射，`counts` 与 `panel_genes` 都以稳定 ID 为键。不同 ID 对应同一 symbol 时仍保留各自计数；符号查询返回 `ambiguous` 和候选 ID，不求和。若某个符号恰好等于另一特征 ID，也会报告歧义。用 `feature_id:ENSMUSG...` 或 `symbol:Cd3d` 显式选择命名空间。

`counts_layer` 默认是 `layers['counts']`；只有明确指定 `"X"` 才读取 X，不会缺层后自动回退。

已有 `slice_id`、`sample_id`、`library_id` 列时会检查单一且非缺失的样本标识，也可以传 `slice_key` 增加指定列的验证；不能借此跳过已有 sample 列。`units` 必须显式为 `micrometer`、`pixel`、`array_index` 或 `unknown`；此读取器不猜测比例尺、不变换坐标。对于没有标签列的对象，可以明确设置 `allow_unannotated: true`。`observation_unit` 可设 `spot` 或 `bin`，并保留相应标签语义。

真实 SPOTS 小鼠脾脏 RNA 输入的完整导入参数（第二片只改 sample/slice 名称与输入文件）：

```json
{
  "slice_id": "spots-spleen-1",
  "sample_id": "spots-spleen-1",
  "coordinate_system": "source_spatial_array_xy",
  "units": "array_index",
  "platform": "SPOTS",
  "observation_unit": "spot",
  "counts_layer": "X",
  "feature_id_key": "gene_ids",
  "allow_unannotated": true,
  "max_features": 40000,
  "max_nnz": 10000000
}
```

本机两片原始 RNA 分别为 2568×32285、2768×32285，非零计数 8,736,679 与 9,360,397；原始 `gene_ids` 唯一，各有 38 个发生重复的 symbol，产生 40 条额外重复行。保留完整 RNA，不以截掉大量基因通过默认预算。源文件不提供物理标尺，`array_index` 不解释成 µm。此例不导入 ADT、不提供真实细胞注释，RNA/ADT 联合建模属于另一个明确工作范围。H5AD 导入还逐行预检序列化来源大小，超过 256 MiB 则拒绝，避免较长特征 ID 导致导入后无法重放；不会自动截断。大型项目导出应选紧凑来源+覆盖层包。

实际新建项目为 `projects/v04-real-data/spots-spleen-1` 和 `spots-spleen-2`；全部计数、坐标和身份逐项来源核对通过。直接 HDF5/SQLite 检查确认完整 ID/symbol 轴相同，两个样本共有 2158 个原条码，但没有相同的 qualified ID。完整来源导出为 221,608,338 与 237,268,121 字节，均通过离线完整性验证。两片按 x 中位数声明的几何选区分别包含 1285 和 1390 个 spot，表达区域结果均可离线重算；`Ptp4a1` 保留歧义，两个显式稳定 ID 分别显示已测量的零。这些是本地软件和输入保真验证，不证明细胞注释、生物学结论或外部宿主采用效果。

HER2ST prepared A1 也已以 `sample_id: "A1"`、`units: "unknown"`、`observation_unit: "spot"`、`counts_layer: "counts"` 完整导入到 `projects/v04-real-data/her2st-A1`：346×12562，796925 个非零条目。几何选区的表达比较和重放成功；微米半径计算明确拒绝。未读取封存真值，未推断比例尺或细胞身份，旧版导入拒绝收据仍保留。

## Xenium：`xenium`

目录需要 `cells.csv` 或 `cells.csv.gz`，其中使用 `cell_id,x_centroid,y_centroid`。原生 x/y 质心单位为微米。矩阵优先读取 `cell_feature_matrix.h5`；没有该文件时读取已解包的 `cell_feature_matrix/`，其中需要 `matrix.mtx[.gz]`、`features.tsv[.gz]`、`barcodes.tsv[.gz]`。不自动解压 tar 包。

```json
{
  "slice_id": "xenium-slice-A",
  "annotations": "C:/data/xenium/annotations.csv",
  "allow_unannotated": false
}
```

```powershell
python -m spatial_collab import xenium C:/data/xenium/outs C:/Spatial/projects/xenium-review --options C:/data/xenium/options.json
```

注释文件示例：

```csv
cell_id,label
aaaa-1,T cell
bbbb-1,Myeloid
```

不需要已有注释时，删除 `annotations` 并设 `allow_unannotated: true`。质心、矩阵和注释的 ID 必须精确一致，行顺序可以不同。只保留类型为 `Gene Expression` 的特征；Xenium MTX 必须提供三列特征类型，不能把未分类控制探针当基因。H5 与 MTX 同时存在时，本版本不额外比较它们的内容是否一致。[10x 官方输出说明](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/tutorials/outputs/xoa-output-understanding-outputs)

### 超出全量导入上限时：明确的 H5 空间窗口

可为 H5 输入提供 `bounds: [xmin, ymin, xmax, ymax]`，使用 Xenium 原生微米坐标。程序扫描完整质心 ID，验证与完整矩阵条形码的精确对应关系，再只读取窗口内细胞对应的 CSC 列。边界包含在内；窗口内所有细胞都会导入，不做随机采样，不根据基因表达挑细胞，不裁减已测量的 `Gene Expression` 面板。

```json
{
  "slice_id": "declared-xenium-slice",
  "bounds": [2700, 2700, 3300, 3300],
  "allow_unannotated": true
}
```

这组坐标只是已声明窗口的示例，不能无条件套用到其它数据。选择窗口应基于科研设计或明确的几何探索，并记录选择原因。若后续在内部 ROI 使用最大半径 r 的邻域统计，应保留至少 r 的外侧观测带；这仍不证明窗外组织没有其它结构。

窗口模式支持最多 2,000,000 个源细胞的有界身份扫描，并限制源 ID 字符总量；实际导入仍限制为 100,000 个细胞、2,000,000 个被选中矩阵存储条目（检查时包括控制特征）。已有全量导入上限保持不变。MTX 窗口读取不在当前实现范围。

`metadata.import_scope` 保存完整边界、原始细胞数、导入细胞数、包含边界的规则、源和所选条目数及无采样声明。窗口项目中的 `whole_slice` 图范围只指**当前导入的完整观察集合**，不能解释为原始整张切片；窗外邻居不可用。注释文件若提供，必须精确覆盖最终窗口内细胞。

已有的原生 QC 列会保存为不可变的观察属性，包括计数、控制探针、基因组控制、细胞面积、核面积、核数量及 `segmentation_method`。缺失字段保持缺失；可选 QC 中的空值、`NaN`、`NA`、`null`、`none` 保存为 null，实际 0 保持为 0。计数矩阵和坐标仍严格要求有限数值，QC 的 ±Infinity 也被拒绝。程序不会把这些字段转换成自行推断的置信度或自动排除规则。

## 通用 10x MTX：`tenx_mtx`

路径直接指向包含 `matrix.mtx[.gz]`、`features.tsv[.gz]` 或 `genes.tsv[.gz]`、`barcodes.tsv[.gz]` 的目录。矩阵方向必须是 feature×observation；只支持 general coordinate Matrix Market 格式。三列特征表按 `Gene Expression` 过滤；传统两列 gene-only 表可导入。相同类型的压缩与非压缩文件同时存在会报歧义。

```json
{
  "positions": "C:/data/mtx/positions.csv",
  "annotations": "C:/data/mtx/annotations.csv",
  "slice_id": "slice-A",
  "coordinate_system": "registered_xy_um",
  "units": "micrometer",
  "platform": "declared platform",
  "observation_unit": "cell",
  "allow_unannotated": false
}
```

位置文件和注释文件都使用 `cell_id` 这个兼容性字段，即使观察单位是 spot 或 bin：

```csv
cell_id,x,y
barcode-A,12.5,20.0
barcode-B,18.0,22.5
```

```csv
cell_id,label
barcode-A,T cell
barcode-B,Myeloid
```

位置和注释必须与矩阵条形码集合完全一致，不会做前缀猜测、数字化 ID、字符串裁剪或位置式连接。[10x 官方 MEX 格式](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-mex-matrices)

## 标准 Visium spots：`visium`

本版本支持传统标准 Visium，**不支持 Visium HD、HD bin、分割细胞或 parquet 位置文件**。需要 `spatial/tissue_positions.csv`，或无表头的旧版 `spatial/tissue_positions_list.csv`，以及 MTX 矩阵目录。矩阵 H5 读取不在此 profile 中。

```json
{
  "slice_id": "visium-slice-A",
  "microns_per_pixel": 0.5,
  "in_tissue_only": true,
  "annotations": "C:/data/visium/in-tissue-annotations.csv",
  "allow_unannotated": false
}
```

```powershell
python -m spatial_collab import visium C:/data/visium/outs C:/Spatial/projects/visium-review --options C:/data/visium/options.json
```

这里的 `0.5` 只是示例，必须替换为**原始全分辨率图像**的真实微米/像素值。计算明确使用 `x = pxl_col_in_fullres × microns_per_pixel`、`y = pxl_row_in_fullres × microns_per_pixel`。不会根据 spot 直径或低分辨率图像比例推断物理单位。[10x 官方空间输出说明](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/spatial-outputs)

`in_tissue_only` 必须明确设置。为 true 时优先 `filtered_feature_bc_matrix/`，缺失时尝试 `raw_feature_bc_matrix/`；为 false 时使用 raw。也可通过 `matrix_dir` 提供绝对路径。矩阵条形码必须精确等于全部位置，或明确请求的 in-tissue 子集；不会自动接受额外 QC 过滤后的不完整集合。

注释 CSV 的 `cell_id` 列放 spot barcode，标签可写区域或工作注释。它必须与**最终导入的 spot 集合**一致：若只导入 in-tissue，需提供对应子集注释。导入同时记录原始位置数、导入数量、明确排除的 off-tissue 数量与校准参数。保留全部 spots 时，通过 `region` 区分 `in_tissue` 和 `off_tissue`。

## CosMx RNA 全局 CSV：`cosmx_csv`

支持 Bruker/NanoString flat-file profile，使用元数据中的 `fov,cell_ID,CenterX_global_px,CenterY_global_px`。`CenterX_local_px`、`CenterY_local_px` 不能直接拼成共享坐标系；此读取器不会替你猜 FOV 配准。

```json
{
  "slice_id": "cosmx-slice-A",
  "metadata_path": "C:/data/cosmx/S0_metadata_file.csv.gz",
  "counts_path": "C:/data/cosmx/S0_exprMat_file.csv.gz",
  "annotations": "C:/data/cosmx/annotations.csv",
  "microns_per_pixel": 0.12,
  "gene_columns": ["CD3D", "LST1", "EPCAM"],
  "exclude_unassigned": true,
  "allow_unannotated": false
}
```

`0.12` 和基因列表都是示例，必须按所选数据集修改。`gene_columns` 是本次导入的明确 RNA 计数范围；其它计数列名称会记录在排除清单中。不能把归一化表、蛋白测量或控制探针列表作为 RNA 基因范围。未明确路径时仅允许唯一匹配的 `*_metadata_file.csv[.gz]` / `*_exprMat_file.csv[.gz]`，或 `metadata.csv` / `exprMat.csv`。

计数宽表需要 `fov,cell_ID` 和指定的基因列。由于 `cell_ID` 只在单个 FOV 内唯一，内部 ID 使用 `fov=<编码值>;cell=<编码值>`，避免不同 FOV 的同号细胞混淆。注释可直接沿用原始双键：

```csv
fov,cell_ID,label
1,1,T cell
2,1,Myeloid
```

也可提供带内部精确 ID 的 `cell_id,label`。`cell_ID=0` 是未分配转录本行；遇到它时默认拒绝，只有明确 `exclude_unassigned: true` 才排除并记录数量。若存在 slide/sample 标识会检查单切片；计数表带有相同上下文字段时必须与元数据一致。[Bruker/NanoString 官方团队格式说明](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/flat-file-exports/flat-files-compare.html)

## MERSCOPE CSV：`merscope_csv`

需要 `cell_metadata.csv[.gz]` 与 `cell_by_gene.csv[.gz]`。默认按 metadata 的 `EntityID` 与计数表的 `cell` 列精确连接；不会按行号拼接。`center_x,center_y` 是原生全局微米坐标。

```json
{
  "slice_id": "merscope-slice-A",
  "metadata_id": "EntityID",
  "counts_id": "cell",
  "annotations": "C:/data/merscope/annotations.csv",
  "allow_unannotated": false
}
```

对于明确采用其它 ID 表头的导出版本，可修改 `metadata_id`、`counts_id`。默认导入计数表中 ID 之外的全部基因列；若导出还包含非基因字段或只需明确的基因子集，可传 `gene_columns` 列表并记录排除列。注释仍使用 `cell_id,label`，其中 cell_id 是 EntityID 的原始字符串。

```csv
EntityID,fov,center_x,center_y
cell-A,1,12.5,20.0
cell-B,2,150.0,40.0
```

```csv
cell,CD3D,LST1
cell-B,0,8
cell-A,5,0
```

这两个文件可以有不同的行顺序，集合必须一致。导入的基因范围以实际 CSV 列为准，不会声称 post-processing 导出的缺失列是已测量的零。[Vizgen MERSCOPE Rev K 用户指南](https://vizgen.com/wp-content/uploads/2025/09/91600001_MERSCOPE_Instrument_User_Guide_RevK.pdf)

## SpatialData Zarr：`spatialdata_zarr`

```json
{
  "table_name": "table",
  "element_name": "cell_centroids",
  "coordinate_system": "global_um",
  "units": "micrometer",
  "slice_id": "slice-A",
  "observation_unit": "cell",
  "platform": "SpatialData",
  "label_key": "cell_type",
  "counts_layer": "counts",
  "allow_unannotated": false
}
```

```powershell
python -m spatial_collab import spatialdata_zarr C:/data/spatial/sample.zarr C:/Spatial/projects/spatialdata-review --options C:/data/spatial/options.json
```

读取明确的 table 和一个 Points 或 Shapes 元素。table 必须声明 `region_key`、`instance_key`；选中 region 的 instance ID 集合必须与全部所选几何对象 ID 精确一致，再按 ID 顺序关联。不会用观察行顺序猜测几何身份，也不会隐式丢弃不匹配的对象。

只应用该几何元素直接保存的、指向指定目标坐标系的有限可逆二维变换。目标坐标系必须已经得到物理校准；名称叫 `global` 本身不证明单位是微米。三维坐标、无效或空 Shapes、raster labels、图像、分割掩膜和转录本归属聚合均不在此导入范围内。几何对象被缩减为质心，不构成对形态或分割质量的验证。

`counts_layer` 可明确写 `"X"`；其它值需要存在于 table.layers。缺失标签列只有在 `allow_unannotated: true` 时允许。导入目标必须在原始 Zarr 目录之外；源文件及成员集合会在导入前后校验。

## 扩展读取器的接口

本地受信任 Python 扩展可以注册能力和延迟加载目标：

```python
from spatial_collab.import_registry import register_reader

register_reader(
    {
        "id": "my_local_reader",
        "description": "Explicit local assay profile",
        "platform": "My assay",
        "observation_units": ["cell"],
        "required_inputs": ["documented local inputs"],
        "required_options": ["slice_id", "coordinate_system"],
        "dependencies": ["my-local-package"],
        "ceilings": ["Document exact scope, limits and identity guarantees"],
    },
    "my_local_package.reader:import_project",
)
```

目标签名为 `reader(path, destination, **options) -> Project`。注册重复 ID 默认拒绝；受信任程序可明确使用 `replace=True`。这不是自动插件发现：数据文件不能要求执行任意模块，`probe` 也不会加载注册目标。读取器应先完成身份、坐标、计数、资源和来源完整性校验，再创建项目。
