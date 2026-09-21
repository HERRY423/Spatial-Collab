# 本轮交付与验证文件

- `final-source-tests.txt`：最终源码回归，460 passed / 2 skipped。
- `installed-tests.txt`：最终独立安装包回归，449 passed / 4 skipped。
- `installed-protocol.txt`：安装包真实 MCP stdio 往返与正常退出。
- `validation-balanced_pca.json`、`validation-smopca.json`：首次成功的真实公开数据执行记录。
- `public-integrity-cross-build.json`：单线程策略修复前真实离线复算失败，予以保留。
- `smopca-repeatability.json`：单线程跨安装环境与同环境重复性诊断。
- `public-integrity-replay.json`：修复后原始数据逐值核对与 4 份真实结果完整离线复算。
- `performance.json`：开发者合成的 12,001 对象分页和单 assay 读取检查，不是科研人员效率数据。
- `browser-review-final.png`：实际浏览器双视图截图。
- `user-pilot/README.md`：可交给组织者的中文指南、任务卡、空白记录表与评分模板；实际参与者 0。

中间诊断文件 `clean-tests*.txt` 不是最终验收：其中一次运行继承了源码路径，另一次在科学库线程加载处中断。它们不能作为独立安装包通过的证据；最终结果仅以上述明确列出的文件为准。

完整解释及能力边界见 `docs/V06_VALIDATION.md`。安装包在项目 `dist` 目录，本轮未发布或推送。
