# v0.7 科学分析：从原始计数直接执行

本版针对六类算法缺口新增独立的科学分析工作流。旧版的结果登记接口和审阅对象仍可使用；新工作流会实际调用算法、保存不可变输入与结果，并记录失败。可用包、成功运行、统计正确性和生物学验证是不同状态。

## 六类任务及输出

| 任务 | 直接执行的方法 | 输出与适用边界 |
|---|---|---|
| 多模态表示 | 官方 MOFA+；官方 MEFISTO sparse GP | 同对象 RNA + 蛋白因子、载荷、单模态 PCA 对照。MOFA 不使用空间；MEFISTO 使用二维坐标。当前适配要求完整配对蛋白，拒绝隐式填补。 |
| 单细胞参考去卷积 | NNLS；官方 Cell2location | NNLS 的 RNA 贡献比例与残差；Cell2location 的丰度均值及 5%/95% 后验分位数。参考必须有明确物种、原始计数和来源标签，不能用空间聚类冒充。 |
| 跨样本 | 官方 Harmony；官方 PASTE | Harmony 只校正表达表示，保留校正前表示；PASTE 输出完整传输质量、熵和推断对应坐标，保留原坐标。两者用途不同。 |
| 空间结构域 | 官方 SpaGCN；稀疏空间图谱聚类基线 | 空间邻接与表达共同参与；SpaGCN 当前为坐标 + 表达模式，不假称使用了组织图像。簇不自动成为生物学注释。 |
| 空间变异基因 | Moran's I、随机置换、BH | 全面板库大小归一化后计算；所有接受检验的基因统一校正，再排序展示；低检出、常量和未检验保留。旧版 Moran 排序不等于此推断流程。 |
| 通讯与通路 | LIANA 物种特异共识库 + 明示空间置换统计；官方 decoupler ULM + PROGENy | 复合体要求全部亚基可测。LR 是空间 RNA 共表达假设，并非复现 CellChat/CellPhoneDB，也不证明信号传递。PROGENy 是独立的 RNA 通路活性推断，不证明特定 LR 的因果作用。 |

## 数据身份与可复算性

分析输入冻结原始稀疏计数、观察 ID 顺序、基因 ID、生产者提供的基因符号映射、物种、坐标、样本及来源。参考与空间输入使用相同稳定基因 ID 才能匹配。不会通过修改大小写、猜测 Ensembl 编号或取第一个重名基因配对。数据库分析排除歧义符号，保留不可检验状态。

计数以校验和绑定的压缩 CSR 矩阵保存，可随审阅包迁移。结果区分 RNA 贡献比例、模型丰度、计算配准、因子和统计量；不覆盖计数、坐标或正式注释。

`submit_analysis` 在独立进程执行并保存任务；`list_analysis_jobs` 可在重连后恢复，失败和取消记录保留，重试产生新任务。缓存绑定输入、方法参数、蛋白层、实现代码及配置环境。代码或环境配置在排队/运行过程中变化会导致任务拒收，不能假称同一次计算。

取消在算法检查点生效；官方单次训练调用尚不支持逐轮中断。两小时预算在检查点检查，不是操作系统级硬截止。大型生产任务仍需外部作业调度。

## 使用步骤

已验证环境为 Windows x64 / Python 3.13，完整版本快照在 `requirements/scientific-py313-win64.lock.txt`。建议新建独立环境，避免宿主环境中的 OpenMP/深度学习库冲突：

```powershell
py -3.13 -m venv .venv-v07
.venv-v07\Scripts\python -m pip install -r requirements/scientific-py313-win64.lock.txt
.venv-v07\Scripts\python -m pip install --no-deps .
.venv-v07\Scripts\python -m pip check
```

锁文件是本次 Windows 验证环境快照，不是跨系统通用保证。SpaGCN 的 `louvain==0.8.2` 在本机通过 MSVC 编译安装，缺少兼容 wheel 的 Windows 主机需要 C++ 构建工具。基础查看器不要求安装全部科学依赖；缺失方法在界面显示不可用。Harmony 固定使用官方 `0.0.10`。PASTE 1.4.0 与 POT 0.9.6.post1 的线搜索回调参数数量不同，本项目提供仅针对这两个版本的签名兼容桥，恢复原函数且不改目标函数，运行结果明确记录是否启用。

1. 使用下方验证环境安装后，为项目配置分析解释器。
2. 在“整合与复核 → 运行科学分析”中确认物种，冻结当前版本或精确 ROI 的 RNA 输入。
3. 需要参考或第二切片时，先在本地明确导入原始计数层、基因轴、样本及物种。
4. 选择方法、输入和参数。Cell2location 需要依据组织学提供每位置预期细胞数；Harmony 要填写每个样本的条件与批次；PASTE 要确认完整重叠假设。
5. 提交后查看状态，可取消、重试或恢复历史任务。结果分页和地图明确标注显示范围，点击对应当前版本的对象可回到原始 RNA 证据。
6. 导出审阅包；默认离线验证检查输入及结果结构。`verify-bundle --include-workflows` 才重新训练/计算，需相应依赖和计算时间。随机训练的复算不一定达到逐元素精确一致，失败须保留，不自动放宽容差。

```powershell
# 在仓库根目录，分析环境创建后执行
.venv-v07\Scripts\python -m spatial_collab analysis configure --project projects/my-study --python .venv-v07/Scripts/python.exe

# 生产者 var.gene_ids 是稳定 ID，var index 是符号，X 经来源确认确实是原始计数
.venv-v07\Scripts\python -m spatial_collab analysis import reference.h5ad --project projects/my-study --sample-id reference-1 --species mouse --kind reference --counts-layer X --label-key cell_types --feature-id-key gene_ids --feature-symbol-key _index

# 若原始计数在 layers['counts']，必须传 --counts-layer counts
.venv-v07\Scripts\python -m spatial_collab analysis run --project projects/my-study --spec recipe.json --background
.venv-v07\Scripts\python -m spatial_collab analysis result RESULT_ID --project projects/my-study --output full-result.json
```

配方示例（ID 必须替换为该项目登记的真实 ID）：

```json
{"method":"nnls","input_ids":["SPATIAL_INPUT_ID","REFERENCE_INPUT_ID"],"parameters":{"n_features":512,"min_reference_cells":5},"seed":0}
```

同一 `method/input_ids/parameters/seed` 对象可交给 MCP 的 `submit_analysis`。Agent 不接收任意系统路径或执行命令；文件导入和环境配置通过本地操作完成。

## 统计与规模边界

- 输入预算：每轴最多 100,000，最多 25M 非零计数，压缩矩阵最多 96 MiB。预算不是已验证的性能声明。
- 空间图 CSR 存储；已用 20,000 位置检查存储需求。表达降维和训练仍有独立预算，不能由此声称十万位置完整工作流已通过。
- 普通稠密表达矩阵最多 20M 项；MEFISTO 最多 10,000 位置且使用稀疏 GP；SpaGCN 官方稠密邻接当前限 6,000；PASTE 最多 4M 位置对。超限明确失败，研究者选择 ROI 后建立新输入，不默默抽样。
- PASTE 的完整重叠不适用于部分重叠切片；坐标单位必须一致。对应坐标不能当作测量事实。
- Harmony 条件与批次完全混杂时拒绝运行；混合度改善不是生物学保持的证据。跨患者显著性推断仍须独立研究设计。
- Cell2location 当前参考签名来自原始参考标签的平均计数，未拟合参考 NB 模型，也未传播参考签名不确定性。每位置细胞数先验、参考缺失类型和训练收敛必须另行评估。
- 置换 p 值最小为 `1/(B+1)`。19 次适于运行检查，不能支持大量关系的细粒度 FDR。默认 99 次同样需要结合检验数评估分辨率。
- PROGENy 网络由官方接口按物种和声明许可范围取得，结果保存完整权重网络及哈希。默认学术范围，商业使用需改为相应范围。网络下载失败不换成编造通路。
- 小鼠 PROGENy 来自官方 HCOP 人到鼠同源转换，不是小鼠扰动实验重新训练的网络。审阅包复算使用结果中冻结的网络，避免重新下载导致权重漂移。

## 官方实现依据

[MOFA+ / MEFISTO](https://biofam.github.io/MOFA2/MEFISTO.html)、[Cell2location](https://cell2location.readthedocs.io/)、[Harmony](https://github.com/slowkow/harmonypy)、[PASTE](https://paste-bio.readthedocs.io/)、[SpaGCN](https://github.com/jianhuupenn/SpaGCN)、[LIANA 共识资源](https://liana.readthedocs.io/en/stable/tutorials/notebooks/prior_knowledge.html)、[PROGENy 资源](https://decoupler.readthedocs.io/en/latest/api/generated/decoupler.op.progeny.html)、[ULM](https://decoupler.readthedocs.io/en/latest/api/generated/decoupler.mt.ulm.html)。算法接入不意味着在所有组织、平台或任务上达到最优水平。

运行证据、固定环境及已知限制见 [v0.7 验证记录](V07_VALIDATION.md)。历史 v0.6 独立用户试用包仍可用，尚无实际独立参与者，不宣称节省研究者时间已经得到验证。
