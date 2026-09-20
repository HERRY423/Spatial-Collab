# v0.5 本地浏览器验收

2026-09-20，独立 Chromium + 实际本机服务。非 ChatGPT iframe 或 Claude 模型任务验收。

- SPOTS spleen 1：主页实际显示 2,568 spots、32,285 RNA 特征、1 个蛋白层、1,285 个共享 ROI 对象。`v05-home.yml`、`v05-home.png`。
- 顶部切入蛋白页，选择 `feature_id:CD19`，真实 `inspect_protein` 返回并显示已测 1,285 / 纳入 1,285、缺失 0、零 1、均值 79.918。源原始均值 79.91828793774319。
- 指定 Cd79a，经浏览器运行所有 21 通道，结果表 22 行（含表头），CD19 对应值 79.918 存在。`v05-protein-result.yml`、`v05-protein-results.png`。
- 切回 RNA：RNA 显示、蛋白页隐藏；共享 ROI 仍 1,285；未标定阵列坐标的物理分析仍禁用；蛋白 21 行保存结果仍在。没有提交注释。
- 从最终 wheel 独立安装启动 HER2ST A1 RNA-only 项目，进入蛋白页：蛋白层 0，读取和比较按钮 disabled，明确显示当前没有蛋白测量，不能用 RNA 替代。`v05-empty-home.yml` 和对应当次 `.playwright-cli` 快照。
- 最终主工作台改为已验证安装版，`http://127.0.0.1:8773/`，实际项目 `projects/v05-real-data/spots-spleen-1`。测试用 HER2ST 服务关闭；原有用户服务不受影响。

生产 JS 对真实 ToolService 的 DOM 测试另验证未保存假设草稿跨导航保留、共享 ROI 不变、常数相关 unknown、信号变换后旧图清空、保存记录不被改写及不隐式提交注释。DOM fixture 不算浏览器引擎验收。
