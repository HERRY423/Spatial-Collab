# v0.6 空间多组学整合与协作复核

本轮实现对应“空间组学插件建议”对话中的三项优先任务及六项开发建议。它仍是研究 Alpha。独立研究者试验尚未开展；试用包不能替代真实用户结果。

最终测试、公开数据核对、真实审阅包复算及未完成的独立验证见 [验收记录](V06_VALIDATION.md)。现成试用材料见 [中文试用指南](V06_PILOT_GUIDE.md)。

## 三项任务与六项建议的落点

| 对话要求 | 当前实现与入口 | 范围 |
|---|---|---|
| 正式整合结果对象 | `IntegrationResult`，CLI `import-integration`，MCP `register_integration` | 校验项目、源摘要、蛋白层、拟合版本、特征轴和完整输出 ID；外部文件不执行代码 |
| 联合分析完整任务 | `integrate` / `submit_integration`；网页“整合与复核” | RNA-only、protein-only、联合 PCA 基线；可选官方 SMOPCA 参考模型；共享空间图、分歧位置、原始证据和既有修订流程 |
| 第二平台与独立用户 | 公开人扁桃体空间 CITE-seq；`prepare_user_pilot.py` | 第二平台实际运行；独立用户 **0 人**，效率与科学准确性未确立 |
| 1. 结果接入再模型运行 | `integration.py`，`smopca_adapter.py` | 输入逐行对应；方法特异表示可存多个；不把簇编号当同一生物学身份 |
| 2. 模态与对应关系 | `multimodal.py`，`register_multimodal_object` | 同一对象、已知一对一、加权聚合、跨切片计算匹配、无配对关系。除严格同对象外只检查映射，不擅自允许配对相关或联合拟合 |
| 3. 预处理与分子关系 | 不可变派生层及分子关系对象 | 归一化/校正层可有负值；推断层要求模型来源；未测量用 null。校正声明要求提供声明的控制依据。关系分同一产物、复合体、标志物关联、自定义，均不自动核实生物学真实性 |
| 4. 方法与重新拟合敏感性 | `filter_integration`、后台任务、`run_integration_sensitivity` | 复用既有冻结假设计划，保留全部替代方案及未知/失败；无监督标签变化不改变拟合输入；假设性排除拟合到明确子集，不产生人工修订 |
| 5. 规模与环境 | 精确视口分页、分块 NPY 蛋白矩阵、按 ID 读取、完整环境锁 | >10,000 对象不再返回空图；统计使用完整 ROI；稀疏 H5 输入按块转密集，不整体稠密化。单切片对象上限仍存在，不能声称已验证整张高密度组织 |
| 6. 跨样本先研究设计 | `study.py`，`compare_study` | 研究→个体→样本→切片→assay→run；条件、批次、配对 ID、来源、ROI 定义及探索性。按个体汇总多个切片；完全混杂显式阻止条件效应解释，不执行跨患者模型或推断检验 |

## 可复现入口

从仓库根目录，在 Python 3.13 / Windows 环境执行：

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements/integration-py313-win64.lock.txt
.venv\Scripts\python -m pip install --no-deps -e .
.venv\Scripts\python scripts/reproduce_v06.py --workspace projects/v06-public --download
.venv\Scripts\python -m spatial_collab serve --project projects/v06-public/tonsil --port 8775
```

锁文件只对应本次验证的 Python 3.13 / Windows 64 位环境；其他平台应独立解析和验证，不照搬 `pywin32`。默认不需要 SMOPCA；其可选依赖为 scikit-learn，包含在上述锁中。

第二后端：

```powershell
.venv\Scripts\python scripts/reproduce_v06.py --workspace projects/v06-public --backend smopca
```

使用已有 SPOTS 项目时追加 `--additional-project <项目目录>`，无需个人原始数据路径。人扁桃体文件固定到官方提交，下载后必须匹配 SHA256 才可导入：

- 仓库提交：`9d59651d78149520661c57bab8a14d42446e6654`
- 文件：`data/RealDataSample/SpatialCITEseq/Humantonsil_filtered.h5`
- SHA256：`67ee4f78382c53580be49a2ba65e74a3c350834cf7ee08cc3a29e6bd8b11a336`
- 原始矩阵：2,491 × 28,417 RNA，2,491 × 283 ADT；不把 `normalized_*` 当作 raw counts。
- 公开文件共用 `cell` 轴；原始 ID、坐标和原始通道名保留。位置以未标定阵列坐标记录，未猜测微米换算。

## 如何复核一次分析

1. 进入“整合与复核”，选择后端、域数量、随机种子、处理方式和全体/共享 ROI。任务返回后可以继续查看数据；“更新任务状态”读取完成或具体失败。
2. 载入结果，选择 RNA-only、protein-only、joint 在同一组织坐标上对照。颜色在每个分区内独立编号。
3. 查看近邻共簇关系的分歧排序。ARI 和近邻关系比较不受簇编号置换影响，不把低一致性解释成标注错误概率。比较不同结果只使用其共同观测对象，并单独报告左右独有数。
4. 点击对象进入原始 RNA 证据与精确共享选区；蛋白页读取同一选区的测量。保留、暂不判断或在既有提案流程中提出修改。真正建立注释版本仍需研究者确认。
5. 修改后选择“固定模型筛选”或“重新拟合”。前者逐值保留父结果；后者绑定新版本和新拟合输入。ROI 重拟合要求明确当前 ROI，不能隐式重选。
6. 导出审阅包。`verify-bundle` 校验正式对象并重新计算内置后端结果；外部方法输出只校验身份和内容，不伪称复跑了外部算法。

## 数值方法与计算预算

简单基线默认按 raw-count variance 选最多 128 个 RNA 特征；这不是宣称最佳的 HVG/SVG 策略。可在 CLI/MCP 显式传入完整选择轴，最多 512 个 RNA 特征、总计 10M 个选定矩阵元素。RNA 先按完整导入 panel 库大小归一化到 10,000，再 log1p。蛋白默认 log1p 或 asinh/5。每个模态按特征标准化，再除以特征数平方根，连接后 SVD，最后固定种子的 k-means。单模态对照使用对应的预处理矩阵。简单基线不使用空间坐标训练，不声称是新的空间整合方法。

SMOPCA 使用官方未修改 `model.py`，附 MIT 许可，限制 3,000 个位置。其模型输入恢复为逐特征单位方差，不采用简单基线的 block balancing；蛋白先归一化库大小，模型内空间坐标以中位最近邻距离缩放。固定 Matérn ν=1.5、γ=1，不优化 γ；最多 20 次噪声迭代。全部参数和实际依赖版本进入结果。该设置是受限参考适配，不是论文全部流程或最佳参数复现。

首轮将模态平衡矩阵传入 SMOPCA 时，官方求根程序报 `f(a) and f(b) must have different signs`。保留此失败事实；根据官方单位方差前处理契约修正了适配尺度。没有修改官方求根逻辑、丢弃失败输出或用 PCA 冒充 SMOPCA。

任务存于项目 SQLite，单项目串行占用计算槽；成功结果按源、修订、assay、特征/处理选项、方法代码摘要、版本、环境和随机种子缓存。取消是阶段间协作取消，活动的数值计算可能先完成。故障保留错误；停止的 worker 可以取消后重新提交。30 分钟预算在阶段边界检查，并不是对 BLAS 调用的实时硬中断。

真实包复算发现 SMOPCA 对 BLAS 线程策略敏感：默认多线程环境曾出现输出不一致；显式单线程后跨安装环境的最大差值约 5.24e-13，两次相同环境拟合差值为 0，所有分区一致。现所有内置拟合及离线复算在运行时明确限制 BLAS 单线程，环境记录包含 NumPy/SciPy 构建信息及线程策略。失败记录保留，未放宽数值容差或修改官方模型；这也不保证任意平台均可逐值复现。

蛋白超过 100k 矩阵元素默认分块保存（每块最多 2,048 行）；缺失保留 NaN，API 返回 null。现有小型 JSON assay 保持兼容。按 ID 读取仅校验目标记录；JSON 的内容摘要缓存不跳过后续字节校验。数组也验证每块摘要，不能把文件名当可信校验。这里优化了存储和单 assay 访问，不声称 RNA 快照、所有历史结果或导出已经做到流式。

## 独立试用

```powershell
python scripts/prepare_user_pilot.py output/my-independent-pilot
python scripts/summarize_user_pilot.py output/my-independent-pilot/observations.csv
```

六位计划参与者覆盖常规 notebook+查看器、通用 Agent+工具、Agent+Spatial-Collab 的六种顺序。空白表记录实际操作时间、等待和总耗时、正确性、返工、遗漏、复算和资源用量。研究者须在试验前冻结三组不同 ROI、任务和独立评分依据。顺序平衡无法完全消除学习/迁移效应；小规模试用应作为预实验。没有真实参与者时不会生成效率结论。真实 ChatGPT/Claude Code 自然语言对话测试仍需相应宿主和实际试用者。

## 来源

- [SMOPCA 官方仓库](https://github.com/cmhimself/SMOPCA)
- [官方固定提交模型](https://github.com/cmhimself/SMOPCA/blob/9d59651d78149520661c57bab8a14d42446e6654/src/model.py)
- [官方空间 CITE-seq 教程](https://github.com/cmhimself/SMOPCA/blob/9d59651d78149520661c57bab8a14d42446e6654/tutorial/Tutorial-Spatial-CITE-seq.ipynb)
- [SMOPCA 论文及数据可用性](https://link.springer.com/article/10.1186/s13059-025-03576-9)

上述来源用于确认输入、模型和数据获取方式；本轮工程验收没有独立证明论文性能主张，也不重新确认对话内其他论文或期刊分区。
