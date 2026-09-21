# 功效声明与方法责任分层

版本保持 0.7.0a1 / 0.7.0-alpha.1。本次增加设计功能与责任契约，不宣称方法创新或独家竞品优势。

## 从零发现转向可检出性

对随机置换使用 `(b+1)/(B+1)`，最小 p 是 `1/(B+1)`。界面统一采用严格 `q < alpha`。BH 排序第 r 位的充分拒绝条件为 `p_(r) < alpha*r/m`，不是每个检验都必须通过 Bonferroni。

BH 的形式化 FDR 控制依赖独立性或适当的正相关条件；未宣称任意空间基因/共享亚基依赖结构都满足这些条件。这里的离散阈值推导是确定性算法性质，模拟功效也不额外证明真实家族的 FDR 控制。

对于 **300 spots、12,426 项检验、99 次置换、alpha=0.05**：

- 最小 p=0.01；至少到第 **2,486** 个排序位置才可能出现 BH 拒绝。
- 若声明按第 1 位阈值规划，则无论效应多强，都不能通过当前离散 p 值网格；该条件下功效为 0，MDE 为“无有限可达值”，不是效应等于零。
- 第 1 位满足严格阈值至少需要 **248,520** 次置换；当前执行器上限 9,999，因此报告明确标注预算超出现有支持。没有偷偷提高上限或把所需预算当成已运行。
- 若至少 2,486 个检验同时达到 p=0.01，BH 则可能检出它们。因此“不可能检测单个稀疏信号”不等于“所有信号均不可检测”。不能声称固定最小 q=m/(B+1)。

这些是精确的分辨率限制，尚不是效应量功效。仅由 n、m 和 B 无法推出唯一的最小 Moran I 或 LR 强度：空间图、信号结构、噪声、零值、复合体及假定信号稀疏程度都会影响功效。

## 两种显式生成模型

`create_analysis_power_plan` 在观察目标统计量之前，解析输入的完整检验家族、坐标和图，保存不可变 `powerdesign`，再运行模拟并保存 `powerplan`。

- **Moran / gaussian_sar**：`z=(I-rho*A)^(-1) epsilon`，A 为对称归一化空间邻接，epsilon 为独立标准正态。调用生产实现的 Moran 统计及相同单侧置换规则。rho 是模型参数；同时报告生成 Moran I 的中位数与 5%–95% 分位区间。
- **LR / lognormal_neighbor_coupling**：生成单位对数方差的正值配体场 L，以标准化 `W^T L` 构造受体对数均值，rho 控制耦合，残余正态噪声权重为 `sqrt(1-rho²)`；检验与产品相同的 `L^T W R/n`，仅置换 R。报告相对于随机摆放期望 `mean(L)*mean(R)` 的超额分数。

独立模拟重复给出经验功效及点态 95% Wilson 区间。网格 MDE 是“该点及所有更大已测网格点的功效区间下限均达到目标”的最小非零 rho；不插值，不假设未测区间单调，未达到时保留 unknown。点态区间不是多网格同时置信区间，因此 MDE 仍是规划摘要。

**这是以声明 BH 排序位置为条件的充分拒绝阈值功效，不是模拟完整 BH 检验家族的总体功效，也不是患者数功效分析。** 不将 rho 当成 logFC、生化结合强度或普适 Moran I 阈值。模型尚未校准真实计数分布、基因掉零、组织异质性、亚基测量误差和空间技术偏差。评估真实实验设计时需提供相符的先验模型/外部 pilot；当前实现不会以观察到的效应计算“事后功效”。

检验家族与生产路径一致：Moran 保留检测过滤且归一化后非常数的实测特征；LR 使用物种数据库和所有可唯一映射的必要亚基。家族、数据库、坐标和输入均绑定哈希，运行前不按目标 p 值筛选。局部结果历史已存在时，功效声明标为 `retrospective_sensitivity`；未存在只表示本工作区顺序，不等于外部预注册或未见数据。

## 使用

网页“整合与复核”的科学分析区选择 Moran/LR，展开“运行前：检验分辨率与条件功效声明”，填写排序假设、目标功效、模拟次数、rho 网格和理由，冻结后运行。更换输入、邻接、家族或置换数后，旧计划被拒绝绑定，需重新冻结。

```json
{
  "spec": {"method":"moran_svg","input_ids":["analysisinput_ID"],"parameters":{"permutations":99},"seed":0},
  "assumptions": {"model":"gaussian_sar","bh_rank":1,"target_power":0.8,"simulations":200,"effect_grid":[0,0.2,0.4,0.6,0.8,0.95],"seed":0,"rationale":"依据先验预计稀疏空间信号，采用最保守排序阈值规划"}
}
```

MCP：`create_analysis_power_plan(spec, assumptions)` → 返回 `powerplan_ID`；`get_analysis_power_plan(plan_id)` 读取；把 `power_plan_id` 加到 `submit_analysis` 的 spec。CLI：`python -m spatial_collab analysis power --project PROJECT --spec recipe.json --assumptions assumptions.json`。无需功效计划也能分析，但结果必须附精确分辨率和“未估计 MDE”，不会把零发现解释为无效应。审阅包保存声明、曲线、失败声明及目标结果，离线验证可重算功效。

MCP/CLI 使用同一配置科学环境的子进程评估功效，避免把核心宿主缺少 LIANA 误当成配置环境不可用。子进程超时为 60 秒；失败或超时不会覆盖已冻结声明。已提交、排队或失败的同输入/方法任务同样触发回顾性标签。实际数值与验收见 [验证报告](POWER_VALIDATION.md)。

## 方法层级与责任

| 层级 | 当前方法 | 项目负责什么 |
|---|---|---|
| 内置数值方法 | Moran 置换/FDR、参考 NNLS、空间 LR 自定义统计；另有 RNA/蛋白区域比较 | 指定统计量、预处理、置换与 FDR 规则、数值目标和独立数值复算；不背书生物学真值 |
| 外部方法适配器 | MOFA+、MEFISTO、Cell2location、Harmony、PASTE、SpaGCN、scikit-learn 谱聚类、decoupler PROGENy/ULM | 冻结输入、明确参数、真实调用、环境清单哈希、保存完整输出、轴/有限值/来源检查与解释边界；不为上游算法或模型的科学正确性背书 |

11 个科学方法条目保持不变；区域比较在原有工具中，不伪称第 12 个新方法。谱聚类虽然包含自定义亲和图，但核心聚类委托给 scikit-learn，保守归入适配层。NNLS 使用 SciPy 求解器，项目背书的是声明的非负最小二乘目标、残差和 RNA 贡献换算，不是自研求解器或绝对细胞数。LR 的 LIANA 数据库依赖与自定义检验分开说明。

目录、界面分组、新结果、分页结果和执行收据使用同一 `method_contract`。收据含环境清单哈希、输入哈希、输出哈希和执行实现哈希；验证器拒绝不一致。环境哈希覆盖记录清单，并不等于对每个安装二进制的校验。`external_runner.py` 的导入桥始终标为外部适配，即使调用者把结果命名为 Moran，也不会获得内置方法身份。历史结果保留为 `legacy_unclassified`，不追认新的责任或先验声明。

依据：[SciPy permutation_test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html)（随机置换 plus-one 规则）、[statsmodels 多重检验实现](https://github.com/statsmodels/statsmodels/blob/main/statsmodels/stats/multitest.py)（BH）。功效生成模型是本项目明确声明的规划模型，不冒充这些上游库提供的组织特异功效保证。
