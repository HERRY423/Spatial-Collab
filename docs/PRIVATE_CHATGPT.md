# 私有 ChatGPT 接入与交付

版本：0.8.0-alpha.1。研究用途。

## 当前真实状态

用户在 2026-09-21 明确告知账号目前没有 ChatGPT 开发模式或隧道权限。因此云端 ChatGPT 接入状态为 **BLOCKED_ACCOUNT_PERMISSIONS**。没有创建隧道、录入 API 密钥或记录真实 ChatGPT iframe 验收。桌面插件安装、真实 MCP 协议与模拟宿主桥接验收不能替代这一步。

当前实现提供独立运行环境配置、空白私有项目、私有桌面插件启动、受校验的文件往返、显式项目路由、OAuth 资源服务器、部署与验收脚本。生产身份提供方、HTTPS 域名与云端账号授权由实际运营者提供，不生成虚构配置。

## Windows 首次安装

从解压的插件目录运行：

```powershell
python scripts/setup_private.py --create-empty --with-import
```

已存在项目时，将 `--create-empty` 换成 `--project 'C:\path\to\project'`。Python 3.11+ 是首次引导依赖，可直接用其完整路径执行脚本。脚本创建专用虚拟环境，不修改系统 Python 或 PowerShell 执行策略；首次安装需要访问 Python 包索引。可用 `--wheel <已核对的 wheel>` 安装确定的交付包。科学模型运行环境仍按现有方法文档单独配置，不默默安装 GPU 模型或承诺模型已经验证。

默认状态目录为 `%LOCALAPPDATA%\SpatialCollab`，包含运行环境、`runtime.json` 和明确为空的 private 项目。可以用 `--home` 指定其他目录；启动插件/隧道前，相应进程需设置 `SPATIAL_COLLAB_HOME`。不会将合成 demo 冒充用户数据。

Windows 插件通过 `scripts/launch_private.cmd` 读取独立环境和项目。`scripts\launch_private.cmd --doctor` 显示运行状态与缺失项，不读取或输出密钥内容。已有 `scripts/run_server.py` 和 Claude 的 stdio 入口仍保留。点击安装只完成插件登记；首次依赖与项目初始化须执行上述安装器，未初始化时会给出明确操作说明。

## 文件进入和取回

工作台新增“我的项目与文件”。上传来自用户主动选择的文件，不接受任意磁盘路径或任意 URL。JSON 文件必须含 `cells` 与 `metadata`；H5AD 需明确原始 counts 层、切片、坐标系、单位和 cell/spot/bin 观测单位。不会猜测 X 是原始计数，不会将 spot 当成纯细胞。

- 单次上传最大 128 MiB，JSON 解析最大 32 MiB；块大小 256 KiB，长度与完整 SHA256 校验。
- 上传采用顺序偏移，允许相同内容重试；重新选择同一文件可续传。宿主若禁用 sessionStorage，仍可使用 `get_transfer` 查偏移，或重新上传。
- 新导入生成独立项目，主项目不被覆盖。原文件复制到该项目的 sources 目录。
- `list_projects` 只列出此连接拥有的项目。每次调用显式携带不透明 `project_id`；不能通过项目参数访问任意本地目录。
- 导出 ZIP 经分块取回和浏览器 SHA256 检查后显示保存链接。浏览器下载预算 128 MiB。宿主禁止 Blob 下载时，可在本地工作台保存；未将这个回退冒充 ChatGPT 附件集成。
- 暂存配额 512 MiB，暂存句柄 24 小时后失效。过期文件在下次上传/导出时清理；可显式删除暂存副本。已导入的原件和审阅导出保留，不受暂存清理影响。

大型原生数据仍使用已有本地流式导入路径，不强行经过聊天或浏览器上传。

## 开通权限后的私有云端连接

需要 ChatGPT 开发模式、Platform Tunnels Read/Use 权限，以及创建隧道时的 Manage 权限；隧道应关联目标 ChatGPT 工作区。按官方 Platform 页面创建并获取 tunnel_id，安装官方 `tunnel-client`。运行密钥只放在安全的进程环境，不进入源码、日志、清单或聊天。

```powershell
# CONTROL_PLANE_API_KEY 通过本地安全方式设置，勿写入本文件。
python scripts/connect_private.py --tunnel-id <实际 tunnel_id> --run
```

此脚本使用官方客户端的 `init` 和 `doctor`，然后可前台保持隧道运行；不会申请账号权限或假造 ID。按官方指引，在 ChatGPT Plugins 的开发连接选择 Tunnel，使用相同 ID；随后在新对话调用 `open_project`。

官方依据：

- https://developers.openai.com/api/docs/guides/secure-mcp-tunnels
- https://developers.openai.com/plugins/deploy/connect-chatgpt
- https://developers.openai.com/plugins/build/chatgpt-ui

## 真实宿主验收

使用单独的合成/公开测试项目，记录实际账号界面、工具结果和失败，不填造研究者观察。

1. ChatGPT 发现工具，`open_project` 显示真正 iframe。
2. 用户框选后，对话中的 `get_selection` 取得相同精确 IDs、坐标框架和 revision。
3. ChatGPT 提出修改，确认前 head 不变；用户明确确认后只产生预览中的修改。
4. 历史提案因过期被拒绝，RNA-only 数据保留无蛋白、未测量 marker 保留未知，未知尺度拒绝微米邻域计算。
5. 导入文件产生新项目，切换后所有调用绑定相同 project_id，旧项目保持原样。
6. 下载包在用户端取得且 SHA256 一致；独立复算通过。
7. 关闭并重开界面、断线重连后，重新读取服务端状态，保留任务失败与取消记录。
8. 跨项目、未授权身份和只读写入请求被拒绝；不得将人为填入的 reviewer 当成认证身份。

本地 `test_host_bridge.py` 运行真实前端代码和实际 ToolService，模拟 postMessage 宿主。它验证工程通信契约，不验证 ChatGPT 的渲染、安全沙箱和下载策略。

## 受认证的远程入口

`python -m spatial_collab.gateway --config <运营者配置> --port 8790` 提供 OAuth 资源服务器。使用 `[gateway]` 可选依赖。仅监听 127.0.0.1，应位于实际 TLS 反向代理之后，原本的 local-only 工作台限制保持原样。

配置字段：`public_origin`、`issuer`、`jwks_url` 和 `projects`；每个项目包含可信本地 `path`、按 OAuth `sub` 明确声明的 `subjects`（read/write），可选 `max_jobs`。项目根目录不得重叠。每个项目端点为 `/projects/<alias>/mcp`，且它必须与访问令牌的 audience 精确匹配。

资源服务器验证 RS256 签名、issuer、audience、subject、iat、exp（最长一小时）、scope；写入需要项目 write 权限及 spatial:write。授权服务需自行支持 ChatGPT 的 OAuth/PKCE/客户端注册流程。项目 ACL 更改需重启加载，身份服务的即时令牌撤销需要运营者另行接入；本实现依靠过期时间限制既有 JWT 生命周期。

提供标准 protected-resource 元数据、有限请求体和并发请求、每项目任务数量限制，以及不含令牌/文件内容的访问日志。不同账号若被明确授予同一项目权限，会共享该项目的选区与修订。子项目继承所属连接项目权限，不是未声明的跨租户共享。

这些能力已经有本地负例检查，但尚未与真实身份提供方、公共域名或 ChatGPT OAuth 对话联调，不能称为生产服务已部署或安全审计完成。

## 数据处理

工具摘要、所选 IDs、坐标和检查结果可被宿主模型处理。文件选择器主动上传的文件进入本次连接的服务，私有部署时在本地保存；不能将“私有存储”宣传为“没有任何数据发送给模型”。不默认将基因矩阵塞入聊天。没有任何上架、自动公开或临床用途承诺。

## 发布边界

当前阶段仅私有/本地安装。图标与 Scientific Research 分类已补齐，但公开官网、支持、隐私与条款 URL、发布者认证和审核流程必须由实际发布者提供。没有虚构联系方式或把本地 Markdown 当成有效公共 URL。
