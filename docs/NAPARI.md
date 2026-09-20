# napari 与 Agent 共用一个研究对象

Spatial Collab 的 napari 适配器读取现有项目，使用 napari 自己的 `Points`、`Shapes`、图层和 dock。ChatGPT、Claude Code、本地浏览器和 napari 指向同一项目目录时，共用源数据、精确选区、提案、修订历史与比较结果。

## 安装与打开

在独立 napari 环境中安装本项目的可选依赖；基础 MCP 安装不会导入 Qt 或 napari。

```powershell
python -m pip install '.[napari]'
python -m napari
```

在 napari 的 Plugins 菜单中打开 **Spatial Collab → Shared spatial annotation review**。输入已通过 Spatial Collab 导入命令创建的项目目录，然后点击 **Open / refresh shared project**。也可以把项目目录或其中的 `project.sqlite3` 文件交给 napari 的 Open 功能，读取观测点和 ROI 图层；要进行持久化协作，请再打开 dock 并指定同一项目。Reader 读取已导入且声明坐标的项目，不会自行猜测原始文件的轴顺序、尺度或注释列。

已经打开 napari 时，也可从 napari 控制台添加 dock：

```python
from spatial_collab.napari_adapter import make_dock_widget
dock = make_dock_widget(viewer)
viewer.window.add_dock_widget(dock, name="Spatial Collab")
```

## 一次协作循环

1. **Explore**：在 observations 图层选择点，点击 **Share selected point IDs**。或在 ROI 图层画一个 polygon / rectangle、选中它，点击 **Share selected polygon / rectangle**。服务保存精确 ID、几何、切片、坐标系、来源散列和修订。Agent 可以读取相同选区。
2. 输入 marker 基因并检查选区。未测基因与已测但计数为零的基因有不同状态；marker 汇总包括被排除的观测。指标由共享项目计算，不从画面猜测细胞身份。
3. **Revise**：填写标签或纳入状态及理由，点击 **Preview annotation change**。审阅完整变化列表后，填写研究者姓名并勾选确认，再提交。姓名仅表示用户提供的归属信息。
4. **Compare**：选择 before / after 修订，固定 selection ID、半径、图范围、source / target 标签和描述性阈值。留空 selection ID 明确表示所有导入的观测。比较保存为共享分析记录，Agent 可读取完整结果。修订后保留先前选区 ID，可在同一选区上做前后比较。
5. Agent 或另一窗口提交了修订后，点击刷新。旧图层无法继续提交选区或应用提案。另一客户端仅改变选区时，marker 检查提示刷新；已选定的修订提案始终保留自己的精确选区。
6. 需要撤销时，在 **Revise** 选择 restore target，先预览，再具名确认。恢复生成新修订，保留整个历史。

**Share** 才会持久化选区；在 napari 中拖动点、修改 feature 表或移动图层不会直接修改科研项目。检测到这些变化时，适配器拒绝继续共享选区，并要求刷新。刷新会清空 dock 的临时 ROI 绘图，已持久化的选区记录仍保留在项目中。

如果提案由 Agent 或浏览器创建，在 **Revise → Shared proposal** 填入 Agent 返回的 proposal ID，点击 **Load shared proposal for review**；留空会尝试读取最新提案。dock 会显示原始提案与其完整精确选区，并高亮这个提案对应的点。即使其他客户端已换了活动选区，提交仍只使用提案保存的那些 ID；加载操作不会改写全局活动选区，也不会自动确认。过期或已应用的提案不能再次提交。若当前全局选区与载入提案不同，marker 检查会明确拒绝混用对象；应先与 Agent 协调活动选区，再重新检查并载入提案。请求生成或载入替代提案失败时，旧提案不会继续作为可提交对象。

## 架构与坐标边界

- `Project` 是共享的事务状态与修订模型，`NapariController` 是不依赖 Qt 的适配器，dock 只处理界面事件。`napari.yaml` 采用 npe2 reader 和 widget contributions。
- 项目坐标是 `(x, y)`；napari 点数组显式转换为 `(y, x)`，保留项目的 `micrometer`、`pixel`、`array_index` 或 `unknown` 单位。只有显式 µm 坐标可运行物理半径比较；其余仍可预览、选区、检查标记和提出工作标签。图层带有项目、source hash、切片、坐标系、修订和 snapshot ID。完整显示最多 250,000 个观测，超限会拒绝读取，不做未声明抽样；此限制不构成性能承诺。
- 观测点的顺序、ID、坐标、标签、纳入状态和区域必须与快照一致，points 层必须保持已校准的恒等变换。
- ROI 支持二维有限可逆仿射变换，读取公开的 `data_to_world` 映射并验证后转回项目 `(x, y)`。无效、奇异、非线性、三维或错配 frame 被拒绝；成员关系由后端按质心包含计算，包含边界。
- SQLite 读写和分析在 napari 的 worker 中运行，返回信号更新 Qt；不会从后台线程修改图层，也不会重建整个 viewer。只刷新本 dock 创建的层，保留其他插件的图层。切换来源项目时移除本 dock 的旧来源图像/分割层，防止混看。npe2 注入当前 dock 所属的 viewer，不依赖全局最后打开的窗口。
- 同步是显式刷新，尚无实时订阅或自动合并冲突。只读图像/分割接入见下；不支持分割编辑回写、转录本渲染、3D 配准或任意 napari 图层向科研项目回写。Visium spot 等观测的科学单位由项目导入元数据记录；画成点不使其变成单细胞。

以上结构参考 napari 官方的 [models and events](https://napari.org/stable/developers/architecture/napari_models.html)、[plugin contributions](https://napari.org/stable/plugins/building_a_plugin/guides.html) 和 [thread workers](https://napari.org/stable/guides/threading.html)，是对其公开机制的使用。它不意味着本插件获得 napari 官方验证或科学认可。

## 来源绑定的图像与分割

`assets.py` 接受本地 **2D NPY** 或 **明确 TIFF page**，必须提供项目 source SHA、slice ID、坐标系、单位、`pixel_to_world` 与来源说明；项目记录 sample ID 时，也必须精确提供它。矩阵是 3×3 有限可逆仿射，把像素中心 `(x,y)` 映射到项目 `(x,y)`，不会推断比例尺或把未知单位改写成 µm。每次登记保存原文件 SHA、像素内容 SHA、矩阵、映射和原始来源文件记录；JSON 侧车和渲染副本不会改动原数据或注释修订。

分割的非零正整数标签必须通过显式 `label_to_cell_id` 字典对应确切细胞 ID；0 仅是背景。缺失或部分映射可以加载观察，但未映射标签显示 unknown。映射到未知细胞、不存在的 mask 标签、重复或歧义整数键会被拒绝。质心附近最近像素的标签对应仅是可检查的描述性记录：缺图像、mask、映射或坐标依据时不能判定对应问题，也不生成“合并细胞”“doublet”“真共表达”等结论。

本地登记示例（绑定 JSON 必须先依据来源记录明确填写，而不是复制猜测的尺度）：

```python
import json
from spatial_collab.assets import register_asset
from spatial_collab.store import Project

project = Project("path/to/project")
# JSON keys: kind, source_sha256, slice_id, coordinate_system, units,
# pixel_to_world, registration_note; optional sample_id, tiff_page,
# label_to_cell_id, origin_sources=[{"path": "...", "sha256": "..."}], name.
with open("reviewed_asset_binding.json", encoding="utf-8") as stream:
    declaration = json.load(stream)
asset = register_asset(project, "path/to/image_or_mask.tif", **declaration)
```

在 dock 的 Explore 页依次点击 **Refresh registered asset list → Load registered image / mask**。Points 和图像/Labels 使用同一声明坐标系；选择最多 100 个点后，用 **Inspect selected IDs in visible asset** 查看确切 ID、图像像素位置和 mask 对应。检查前会重新验证来源文件、图层内容和变换。图层被移动、像素替换或来源改变后必须重新检查，不会沿用旧对应。

既有 napari 图层也可用 `capture_asset_layer(layer)` 冻结，再调用 `NapariController.register_layer_snapshot(capture, source_sha256=..., slice_id=..., coordinate_system=..., units=..., registration_note=..., label_to_cell_id=..., origin_sources=...)`。捕获仅接受已物化的二维 NumPy image/labels；lazy、multiscale、RGB、3D 被拒绝，不无界加载。已有图层的公开 `data_to_world` 仿射被明确记录，登记者仍须确认它对应所声明项目。

边界：每个 plane 最多 32,000,000 像素和 256 MiB，源文件最多 2 GiB，每项目最多 32 个资产，mask 最多 250,000 个标签。TIFF 使用可选 `tifffile`；多页必须选择 page，不做投影。大图须先生成有明确来源及位置记录的有界 crop。napari 的原生 Labels 渲染器要求可写缓冲区，因此使用独立渲染副本并设 `editable=False`；原文件和核心加载数组保持只读，控制台强行替换副本仍会被检查拒绝。该行为参考官方 [Labels non-editable mode](https://napari.org/dev/howtos/layers/labels.html)。

## 已验证范围与复现

2026-09-20，本地 headless controller 用例验证了精确 ID、轴转换、marker 状态、显式确认、历史恢复、前后比较、跨客户端 stale 状态、点和 frame 被改动时拒绝提交、ROI 仿射转换及错误变换拒绝。基础环境未安装 napari，真实宿主测试显示为 skip。

同日，在现有独立环境 **Python 3.12.13 / napari 0.8.0 / npe2 0.8.3** 中运行真实 Points / Shapes、npe2 manifest、dock 和后台 worker 集成测试，最终 **22 项通过**。其中真实 dock 用例覆盖双窗口正确注入、Agent 提案加载、原始精确选区高亮、加载后清除确认、点位置被改动时拒绝提交，以及刷新后仅提交提案保存的 ID；controller 用例验证未确认不得提交。使用 `viewer(show=False)` 的 Windows 原生后端；这属于隐藏窗口集成，不是可见界面人工验收。Qt 的 `offscreen` 后端尝试失败，因为该 Windows 后端无法创建 OpenGL 上下文，不能把该尝试计作通过。环境还存在一个与本插件无关的旧 `spatial-evidence-layer` entrypoint 缺失模块告警，以及 npe2 / Pydantic 弃用告警。

复现命令（本机示例，临时目录须选择未使用名称，且不能指向数据目录）：

```powershell
$env:PYTEST_DISABLE_PLUGIN_AUTOLOAD = '1'
$env:PYTHONDONTWRITEBYTECODE = '1'
$env:QT_QPA_PLATFORM = 'windows'
$env:PYTHONPATH = 'C:\Spatial\spatial-collab\src'
$env:NUMBA_CACHE_DIR = 'C:\Spatial\spatial-collab\projects\napari-numba-cache'
$env:NAPARI_CONFIG = 'C:\Spatial\spatial-collab\projects\napari-test-settings.yaml'
& 'C:\Users\13264\anaconda3\envs\napari-env\python.exe' -m pytest tests/test_napari_adapter.py -q -p no:cacheprovider --basetemp=C:\Spatial\spatial-collab\projects\napari-tests-UNUSED
```

测试采用合成数据；不提供真实组织准确率、研究者节省时间、ChatGPT / Claude 会话接入成功或生物学结论成立的证据。单切片描述性敏感性与独立生物学验证始终分开。

图像/分割功能另经现有 napari 环境的原生隐藏窗口验证。新增用例覆盖只读来源、明确 mask ID 映射、图像和点的 XY/YX 变换、错来源/切片/frame/单位拒绝、缺映射 unknown、来源 hash 变化、移动图层/替换像素拒绝、既有图层快照和非物理单位禁用物理分析。首次原生尝试曾报 `buffer source array is read-only`，已改成上述独立渲染副本机制；不能把失败尝试当作通过。

真实格式验收：`scripts/acceptance_assets.py` 读取本机 `D:\Spatial Evidence Layer\acceptance-data\xenium-kidney-v4-tiny\outs` 的匹配 bundle，登记 DAPI focus TIFF page 0 和 `cells.zarr.zip:masks/1` 的无损 NPY 快照。矩阵来自原生 `homogeneous_transform` 的显式逆变换，并与 OME/experiment 的 0.2125 µm 比例尺交叉核对；映射直接读取 `cell_boundaries.csv.gz` 的 `label_id` 与 `cell_id`，不按数组顺序猜测。这些字段的语义参考 [10x 官方 Zarr 输出说明](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/advanced/xoa-output-zarr)。

凭证保存在 `projects/assets-native-acceptance-20260920/native-assets-receipt.json`：358 个真实 tiny 对象，357 个质心最近像素对应相同声明 ID，1 个对应不同声明 ID（`mnlmjgpm-1` 的最近像素 label 221 声明为 `ipcgjdge-1`），差异完整保留。**这不是分割准确率，也不自动解释成分割错误或合并细胞。** 图像/分割/点已在原生隐藏 viewer 加载，移动 mask 被拒绝，共享状态前后相同，无修订或研究者批准。tiny 证明格式及宿主坐标行为；不证明淋巴结项目具备匹配图像、科学可靠性或可见界面人工验收。

复现真实格式验收（同样使用上述环境，输出须为新目录）：

```powershell
& 'C:\Users\13264\anaconda3\envs\napari-env\python.exe' scripts/acceptance_assets.py --outs 'D:\Spatial Evidence Layer\acceptance-data\xenium-kidney-v4-tiny\outs' --output 'C:\Spatial\spatial-collab\projects\assets-native-acceptance-UNUSED'
```
