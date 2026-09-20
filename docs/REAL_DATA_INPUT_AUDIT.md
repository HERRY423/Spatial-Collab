# 本机既有真实数据输入核查（2026-09-20）

这份记录区分实际文件检查、历史获取收据和官方资料。所有原始文件只读；未下载大型数据。局部数据成功读取或本地 SHA-256 一致不等于整套实验输出完整、独立注释真值或生物学有效性验证。机器可读明细在 `projects/real-data-review/discovery.json`。

## 1. 主分析输入：Xenium Prime 人反应性淋巴结

实际输入目录：`C:\Spatial Transcriptomics\histoweave\datasets_cache\raw_sources\xenium`。

| 输入 | 字节 | 本次 SHA-256 |
| --- | ---: | --- |
| `cell_feature_matrix.h5` | 206940548 | `21dd59f14a25fec60d7a41e2f1086c91a430c00a6b0b1878695fdba68f44f5a5` |
| `cells.csv.gz` | 39689688 | `165d132c5717c73f462d9f6af635121bad1caba13c13a900c3307a670f3243b1` |

本次直接读取 HDF5 结构和细胞表，矩阵为 11094 个特征 × 708983 个细胞、157120741 个非零项。708983 个细胞 ID 唯一，细胞表与矩阵 barcode 顺序完全相同。4624 个 `Gene Expression` 特征的名称与 ID 均唯一；其余类型为 3291 Deprecated Codeword、2509 Unassigned Codeword、609 Negative Control Codeword、40 Negative Control Probe、21 Genomic Control，不能混作基因表达的分母。

[官方数据页](https://www.10xgenomics.com/datasets/preview-data-xenium-prime-gene-expression)对应 1 位供者的人 FFPE 反应性淋巴结、708983 个细胞和 4624 个靶基因，是开发版试剂及分析流程的预览数据。官方同时明确标签迁移曾把一群可能的 MARCO+ 巨噬细胞误标为内皮细胞。此处记录这一已知问题作为复核动机，不能据此宣布本次窗口内某个细胞已被证实误标。

`x_centroid/y_centroid` 的原生单位由 [10x 输出规范](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs)定义为 µm。当前本地目录缺少完整的 `experiment.xenium`、形态图、分割和转录本输出，不支持图像配准、分割正确性或单转录本归属审查。旁边残缺 ZIP 的问题不自动使单独可读的矩阵和细胞表失效，但本次仅确认它们的结构、ID 对齐和本地字节指纹，未确认上游发布者校验和。

独立按原生坐标包含边界选择 `[2700,2700,3300,3300]`，得到 7725 个细胞。实际范围 x=2700.05419921875–3299.743408203125、y=2700.008544921875–3299.939697265625。窗口选择不是生物学代表性抽样，也不能生成新的独立供者。矩阵中有 `MARCO/PECAM1/MS4A1`，没有 `LYZ/VWF/CD3D`；应报告靶面板覆盖不足，不能把未测基因作为阴性证据。

### 实际导入检查与发现的来源身份错误

第一份项目 `projects/real-data-review/xenium-window-2700-3300` 的数值独立检查通过：全部 7725 个 ID、原生坐标、逐细胞 Gene Expression 总数、三类 QC 属性与原文件一致，4624 个基因完整；导入前固定的十个细胞（四个坐标极值和六个固定种子随机细胞）的 4624 维原始计数逐项完全相同。检查直接使用 CSV、HDF5 CSC 列和只读 SQLite，没有调用插件导入实现；未把十列抽查扩大表述为全矩阵逐项复核。详见 `independent-import-receipt.json`。

但该项目不可变 `metadata.slice_id` 实际为 `human-colon-cancer-xenium`，与真实反应性淋巴结来源不符。数值通过不能掩盖共享科研对象的样本身份错误；综合来源检查为 **FAIL_SOURCE_IDENTITY_METADATA**。另存 `independent-import-metadata-addendum.json`，保留原数值检查及错误项目，不修改其不可变来源。替代项目应使用正确切片 ID 并重复同一组核验。

随后新建的 `projects/real-data-review/xenium-lymph-node-window` 已独立复核通过，slice ID 为 `xenium-prime-reactive-lymph-node`。全部 ID、原生坐标、基因总计数、三类 QC 属性、4624 个基因及预先固定十列均再次匹配；原始文件哈希未变，新旧项目的完整 cell payload 精确相同。正确项目有独立来源快照和新 ID，旧项目及失败追加记录仍保留。最终收据为 `independent-import-receipt-corrected.json`；它只证明限定范围的输入保真与元数据更正，不证明注释或生物学结论正确。

只读方案复核保存在 `study-plan-independent-review.json`。B/T 表达程序是阈值工作分组，使用同一组定义标记的表达差异有循环性，不能证明身份。`min_effect=0.1` 仅为描述性阈值；需同时显示来源数、目标丰度和邻居分母。中心 ROI 的 100 µm 外围范围大于最大 75 µm 半径，只有保留该外围邻居、以中心 ROI 细胞为来源时才支持完整半径搜索。方案作者声明提前固定规则，本次复核不升级为独立预注册证明。

## 2. SPOTS 小鼠脾脏 RNA + ADT：真实、可配对，但没有物理标尺

目录：`C:\Plugin\ExoWarrant\output\sc-spatial-multiomics-20260911-01\inputs`。

| 输入 | 实际形状 | 字节 | 本次 SHA-256 |
| --- | --- | ---: | --- |
| `spleen1_adata_RNA.h5ad` | 2568 spots × 32285 genes | 73347838 | `afce90ec192aa55767a11f6639f6876055c72c5d97961b7f65c85e2389298322` |
| `spleen1_adata_ADT.h5ad` | 2568 spots × 21 proteins | 446952 | `76ab82fd830d421250bafeb704c3e686cc7ebdee3a854f10acb19bddd33a674a` |
| `spleen2_adata_RNA.h5ad` | 2768 spots × 32285 genes | 78422654 | `3981d22914374306823dad2167fa39fd365df4bfb1aa20f47fabcfde157ecc39` |
| `spleen2_adata_ADT.h5ad` | 2768 spots × 21 proteins | 470152 | `706b4542739d763246d54eac3616ef0385edf8562f027cc48f2b889cdde7a932` |

本次 4 个哈希均与既有 `input-manifest.json` 相同。原始 X 的非零值全部非负且为整数，RNA 是 CSR float32，ADT 是稠密 float32；原始 obs 无标签。两片各有 38 个发生重复的 gene symbol，产生 40 条额外重复行，但 `var['gene_ids']` 的 Ensembl ID 唯一。不能随意给重复符号加序号后当作新的生物特征，也不应静默合并不相同的 Ensembl 特征。

同一切片 RNA/ADT 的 barcode 与坐标逐项完全相等；两个切片之间有 2158 个相同 barcode，需结合 sample/slice ID 保持实体身份。`obsm['spatial']` 为供应方 array 坐标：第 1 片范围 `[0,3]–[66,127]`，第 2 片 `[2,0]–[77,127]`。文件无图像、比例尺、物理单位或配准变换；不得填入虚构的 µm/像素。可进行表达和 spot 注释探索；基于物理半径的邻域结论需要另行取得标尺。

[SpatialGlue 的 Zenodo 数据页](https://zenodo.org/records/10362607)明确包含 2 个 SPOTS 小鼠脾脏实验数据集；既有获取收据将其对应到 [GSE198353](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198353) 和 GSM5945505–GSM5945508。本次 NCBI 网页返回访问挑战，没有绕过，GEO 映射沿用本机 2026-09-11 元数据收据。下载来源 `https://zenodo.org/api/records/10362607/files/Data_SpatialGlue.zip/content`，所取为 ZIP 内的 4 个成员。整包发布 MD5 为 `b766948c1136bb3af2f3e6e85dc81cfe`，既有收据明确整包 MD5 **未核验**，记录的是成员 CRC32/Range 校验；本次不提升这一证据状态。

同包 `results/spleen{1,2}_integrated.h5ad` 保留 counts、gene_symbol、QC、RNA signature、protein_counts/CLR 及聚类。本次全矩阵比较证实 `layers['counts']` 与原始 X 精确一致，barcode 与 Ensembl 特征顺序亦一致。已有 `dominant_reference_signature` 是本机探索性映射，约 58.6% spot 标为模糊，不能作为独立真值或细胞比例。两个文件自身的 `analysis_boundary` 也明确限制这一解释。可用于旧注释复核和研究假说敏感性，不能将旧工具产物误称为人工审阅标签。

## 3. HER2ST 历史输入：真实切片，但标尺仍未确认

目录：`C:\Spatial Transcriptomics\histoweave-her2st-data\prepared\inputs`。主任务的并行检查记录在 `projects/real-data-review/her2st-inspection.json`；本审计另外核对了文件存在性、字节数以及 `prepared/source_manifest.json`。以下形状与新鲜哈希来自该并行检查，不虚报为本子任务独立重算。

| 样本 | spots × genes | 字节 | 核验 SHA-256 |
| --- | --- | ---: | --- |
| A1 | 346 × 12562 | 3865626 | `893a87c46211f76d367dbfcdc4e06e647c9c8c7a60b9faf128ddd9a445b83090` |
| B1 | 295 × 12682 | 2432539 | `eae1586e32ad7763126286beab374cf395c8fc0cee969a9a8897eabf1cccb54d` |
| C1 | 176 × 13416 | 3144982 | `66db06ea2ab15334c00f0de019daa3531db0d321b0df195e34ea8cdd8ff8496e` |
| D1 | 306 × 13814 | 5072508 | `37c57d124f5de7f37ca465877250b66925b2ba0c50e5372e27f2937126d198b1` |
| E1 | 587 × 12179 | 3141905 | `03dfdb9cc2371360d82e3ffc035156561ee838767afd715a86e1afb8b21e748f` |
| F1 | 691 × 12967 | 6207115 | `0e3afd7aa423eac11066d390e95d9bc817f94a985bdb04c4e6334e09d073c65f` |
| H1 | 613 × 12491 | 3095693 | `5aa91ff21bf7b363aa85abba53a3a40de680c2b435f180b75cf257d3fd4885ff` |

共 3014 个 spot，prepared manifest 对应 A/B/C/D/E/F/H 七个 donor。源仓库 remote 为 [almaan/her2st](https://github.com/almaan/her2st)，manifest 记录 commit `85df7411b987a7dec3feff9718b869a7a82091ec`；本地原作者 README 提供 [Zenodo3957257](https://zenodo.org/record/3957257) 为数据入口。这些是已有准备流程的溯源记录，本次没有重新下载或验证完整发布包。

检查确认 counts 为非负整数，spot/gene IDs 唯一；obs 有 sample_id、donor_id、array_row/array_col，未混入标签。没有物理单位或变换声明，A1 的空间坐标 `[9.852,12.991]` 与 array index `[10,13]` 相近也不能证明单位为 µm。因此当前导入拒绝保留为实际失败结果：`BLOCKED_UNVERIFIED_COORDINATE_CALIBRATION`。未填造标签、未猜测阵列间距、未打开封存的 private_truth，也没有宣称已完成此系列的科研分析。G2 的历史排除只沿用 manifest 的重复坐标原因，不重新开启其真值。

## 4. 旧空间项目、Xenium tiny 与排除项

`D:\Spatial Evidence Layer\acceptance-data\xenium_lymph_outs` 是历史 CRC/解压失败后的部分输出。旧项目 `docs/REAL_XENIUM_ACCEPTANCE_2026-08-12.md` 保留失败记录，不能使用其目录名推定完整性。

可读的 `D:\Spatial Evidence Layer\acceptance-data\xenium-kidney-v4-tiny\outs` 是官方人肾脏 tiny 格式测试输入。ZIP 20929485 字节，本次 SHA-256 `abd7e8f7fd047dcc6afdb1e9eece90d4533d3ead053c6f05c482be050bdf79d2`、MD5 `50c04dea5e751e1c7508ff24528242e8` 与官方公布值一致。矩阵为 358 cells × 520 features，其中 405 Gene Expression、4 Protein Expression，其他为对照；`experiment.xenium` 记录像素大小 0.2125 µm。BioNexus `data/flagship/xenium_spatial_truth` 另有相同文件的本地副本。

[官方说明](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/resources/xenium-example-data)限定 XOA v4 tiny 为跨两个 FOV 的三个人工裁剪小块，用于格式测试，不用于生物学结论。因此可测试多模态特征过滤和格式兼容，不能替代真实科研效能评估。

`C:\Plugin\BioNexus\data\flagship\citeseq_pbmc_sorted\pbmc_multimodal_zenodo_10213715.h5ad` 是 PBMC CITE-seq，来源 `https://zenodo.org/records/10213715/files/pbmc_multimodal.h5ad?download=1`。没有将它误计为空间数据，也未为本次任务加载其 1.945 GB 内容。

## 直接暴露的协作需求

- “同一份数据”必须锁定切片、barcode、原始来源和特征 ID；真实 SPOTS 已证明 barcode 跨切片重复。
- 接收未注释的真实数据，明确区分原始标签、旧工具假说、当前修改及人工确认，而不是要求先补出一列貌似确定的 cell type。
- 允许空间点在未标定坐标下进行明确受限的探索；物理距离方法应说明缺失标尺，不应阻断所有表达复核。
- 面板缺失、对照特征、多模态读出和原始计数来源必须进入可见上下文；弱证据不能靠标签名称变成细胞身份。
- 大切片可按精确坐标窗口读取；报告窗口边界和研究限制，并独立检查抽取的 ID、原生坐标及稀疏列计数。
