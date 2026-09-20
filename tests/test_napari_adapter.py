"""Controller contracts run without Qt; optional host tests use real napari."""
from copy import deepcopy
import importlib.util
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from spatial_collab import napari_adapter as adapter
from spatial_collab.store import Project, SpatialError


@pytest.fixture
def controller(tmp_path):
    cells = [{"cell_id": "b", "x": 30, "y": 20, "label": "B", "counts": {"G": 0}},
             {"cell_id": "a", "x": 10, "y": 20, "label": "A", "counts": {"G": 4}},
             {"cell_id": "c", "x": 50, "y": 40, "label": "B", "counts": {}}]
    metadata = {"name": "synthetic test", "slice_id": "s", "coordinate_system": "s-xy",
                "units": "micrometer", "panel_genes": ["G"], "source_kind": "synthetic"}
    project = Project.create(tmp_path / "project", cells, metadata)
    return adapter.NapariController(project.root)


def fake_points(controller, selected=(0,)):
    data, kwargs, _ = controller.layer_data()[0]
    return SimpleNamespace(data=data, features=kwargs["features"], metadata=kwargs["metadata"],
                           selected_data=set(selected), ndim=2, data_to_world=lambda point: point)


def fake_shape(controller, scale=(1, 1), translate=(0, 0)):
    return SimpleNamespace(data=[np.array([[15, 5], [15, 35], [25, 35], [25, 5]])],
                           metadata=deepcopy(controller.layer_data()[1][1]["metadata"]),
                           selected_data={0}, shape_type=["polygon"], ndim=2,
                           data_to_world=lambda point: np.asarray(point) * scale + translate)


def test_import_is_independent_of_qt_and_napari():
    result = subprocess.run([sys.executable, "-c",
                             "import spatial_collab.napari_adapter; import sys; "
                             "assert 'napari' not in sys.modules; assert 'qtpy' not in sys.modules"],
                            env={**__import__("os").environ,
                                 "PYTHONPATH": str(Path(adapter.__file__).parents[1])},
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_reader_exact_axis_identity_and_metadata(controller):
    assert adapter.napari_get_reader(str(controller.project.root)) is adapter.read_project_layers
    assert adapter.napari_get_reader(str(controller.project.db_path)) is adapter.read_project_layers
    assert adapter.napari_get_reader(["one", "two"]) is None
    assert adapter.napari_get_reader("missing.h5ad") is None
    layers = adapter.read_project_layers(str(controller.project.db_path))
    points, props, kind = layers[0]
    assert kind == "points" and layers[1][2] == "shapes"
    assert points.tolist() == [[20, 10], [20, 30], [40, 50]]
    assert props["features"]["cell_id"] == ["a", "b", "c"]
    binding = props["metadata"][adapter.KEY]
    assert binding["revision_id"] == controller.revision_id
    assert binding["axis_order"] == ["y", "x"] and binding["units"] == "micrometer"


def test_exact_selection_marker_proposal_confirmation_and_restore(controller):
    base = controller.revision_id
    selection = controller.publish_points(adapter.capture_points(fake_points(controller)))
    assert selection["cell_ids"] == ["a"]
    markers = controller.inspect(["G", "ABSENT"])["expression"]
    assert markers["G"]["sum"] == 4 and markers["ABSENT"]["status"] == "unmeasured"
    proposal = controller.propose({"label": "B", "included": False}, "researcher marker review")
    assert proposal["deltas"][0]["cell_id"] == "a"
    assert controller.project.summary()["head_revision"] == base
    with pytest.raises(SpatialError, match="confirmation"):
        controller.commit("researcher", False)
    committed = controller.commit("researcher", True)
    assert committed["revision_id"] != base
    assert controller.project.cells()[0]["label"] == "B"
    result = controller.compare(base, committed["revision_id"], selection["selection_id"],
                                radius_um=25, source_label="A", target_label="B", graph_scope="whole_slice")
    assert result["comparison"]["status"] == "indeterminate"
    assert result["scientific_authorization"] == "NOT_ESTABLISHED"
    restored = controller.restore(base, "researcher", True)
    assert restored["revision_id"] not in {base, committed["revision_id"]}
    assert controller.project.cells()[0]["label"] == "A"


def test_polygon_applies_known_transform_then_yx_conversion(controller):
    # Displayed world y=15..25, x=5..35, selecting a and b.
    shape = fake_shape(controller, scale=(2, 3), translate=(1, -2))
    shape.data = [(shape.data[0] - (1, -2)) / (2, 3)]
    capture = adapter.capture_polygon(shape)
    assert capture["polygon"] == [[5, 15], [35, 15], [35, 25], [5, 25]]
    result = controller.publish_polygon(capture, adapter.capture_points(fake_points(controller)))
    assert result["cell_ids"] == ["a", "b"]
    assert result["method"] == "centroid_containment_boundary_included"


@pytest.mark.parametrize("mutation", ["move", "reorder", "add", "identity", "label", "metadata", "transform"])
def test_tampered_points_never_publish(controller, mutation):
    layer = fake_points(controller)
    if mutation == "move":
        layer.data[0, 0] += 1
    elif mutation == "reorder":
        layer.data = layer.data[::-1]
    elif mutation == "add":
        layer.data = np.vstack([layer.data, [0, 0]])
    elif mutation == "identity":
        layer.features["cell_id"][0] = "c"
    elif mutation == "label":
        layer.features["label"][0] = "new"
    elif mutation == "metadata":
        layer.metadata[adapter.KEY]["units"] = "pixel"
    else:
        layer.data_to_world = lambda point: np.asarray(point) + 100
    with pytest.raises(SpatialError):
        controller.publish_points(adapter.capture_points(layer))
    assert controller.project.get_selection() is None


def test_stale_revision_and_cross_client_selection_are_explicit(controller):
    points = adapter.capture_points(fake_points(controller))
    controller.publish_points(points)
    other = Project(controller.project.root)
    changed_selection = other.set_selection(controller.revision_id, cell_ids=["b"])
    with pytest.raises(SpatialError, match="changed the selection"):
        controller.inspect()
    # An already chosen exact ID set remains the proposal's selection, never the
    # other client's newly active selection.
    assert controller.propose({"label": "B"}, "review a")["deltas"][0]["cell_id"] == "a"
    proposal = other.propose_revision(controller.revision_id, changed_selection["selection_id"],
                                     {"label": "A"}, "external edit")
    other.apply_revision(proposal["proposal_id"], controller.revision_id, "other", True)
    with pytest.raises(SpatialError, match="Stale"):
        controller.commit("me", True)
    with pytest.raises(SpatialError, match="Stale"):
        controller.publish_points(points)
    controller.refresh()
    assert controller.preview is None
    assert controller.features["label"] == ["A", "A", "B"]


def test_load_agent_proposal_keeps_its_exact_selection_and_needs_confirmation(controller):
    agent = Project(controller.project.root)
    selection = agent.set_selection(controller.revision_id, cell_ids=["b", "c"])
    proposal = agent.propose_revision(controller.revision_id, selection["selection_id"],
                                      {"label": "A"}, "agent proposes two cells", actor="agent")
    # An unrelated active selection must not redirect the reviewed proposal.
    other_selection = agent.set_selection(controller.revision_id, cell_ids=["a"])
    loaded = controller.load_proposal()
    assert loaded["proposal"]["proposal_id"] == proposal["proposal_id"]
    assert loaded["selection"]["cell_ids"] == ["b", "c"]
    assert not loaded["shared_selection_changed"]
    assert agent.get_selection()["selection_id"] == other_selection["selection_id"]
    with pytest.raises(SpatialError, match="changed the selection"):
        controller.inspect(["G"])
    with pytest.raises(SpatialError, match="confirmation"):
        controller.commit("researcher", False)
    committed = controller.commit("researcher", True)
    assert set(committed["overlays"]) == {"b", "c"}
    with pytest.raises(SpatialError, match="stale|already applied"):
        controller.load_proposal(proposal["proposal_id"])
    assert controller.preview is None


def test_failed_replacement_preview_clears_old_intent(controller):
    controller.publish_points(adapter.capture_points(fake_points(controller)))
    controller.propose({"label": "B"}, "initial intent")
    with pytest.raises(SpatialError):
        controller.propose({}, "invalid replacement")
    with pytest.raises(SpatialError, match="preview first"):
        controller.commit("researcher", True)
    controller.propose({"label": "B"}, "initial intent")
    with pytest.raises(SpatialError):
        controller.load_proposal("unknown-proposal")
    with pytest.raises(SpatialError, match="preview first"):
        controller.commit("researcher", True)


@pytest.mark.parametrize("mutation", ["3d", "singular", "nonlinear", "ellipse", "unselected", "wrong_frame"])
def test_unsupported_shape_fails_without_selection_write(controller, mutation):
    shape = fake_shape(controller)
    if mutation == "3d":
        shape.ndim = 3
    elif mutation == "singular":
        shape.data_to_world = lambda point: (0, 0)
    elif mutation == "nonlinear":
        shape.data_to_world = lambda point: np.asarray(point) ** 2
    elif mutation == "ellipse":
        shape.shape_type = ["ellipse"]
    elif mutation == "unselected":
        shape.selected_data = set()
    else:
        shape.metadata[adapter.KEY]["coordinate_system"] = "unregistered image pixels"
    with pytest.raises(SpatialError):
        controller.publish_polygon(adapter.capture_polygon(shape), adapter.capture_points(fake_points(controller)))
    assert controller.project.get_selection() is None


def test_resource_limit_never_silently_samples(controller, monkeypatch):
    monkeypatch.setattr(adapter, "MAX_LAYER_CELLS", 2)
    with pytest.raises(SpatialError, match="no cells were sampled"):
        controller.refresh()


@pytest.mark.skipif(importlib.util.find_spec("napari") is None, reason="optional real napari runtime unavailable")
def test_real_napari_layers_and_dock(controller):
    # Use a native hidden viewer; Windows Qt offscreen lacks an OpenGL context.
    # This is host integration, not visible GUI usability acceptance.
    import napari
    from npe2 import PluginManifest, PluginManager
    from qtpy.QtWidgets import QApplication
    manifest = PluginManifest.from_file(Path(adapter.__file__).with_name("napari.yaml"))
    assert manifest.name == "spatial-collab"
    viewer = napari.Viewer(show=False)
    second_viewer = napari.Viewer(show=False)
    manager = PluginManager.instance()
    if manifest.name in manager:
        manager.unregister(manifest.name)
    manager.register(manifest)
    try:
        unrelated = viewer.add_points([[1, 2]], name="unrelated")
        _, dock = viewer.window.add_plugin_dock_widget(manifest.name)
        dock.opened(controller)
        assert unrelated in viewer.layers and len(viewer.layers) == 3
        assert len(second_viewer.layers) == 0  # actual viewer injection, not global current_viewer
        dock.points.selected_data = {0}
        assert controller.publish_points(adapter.capture_points(dock.points))["cell_ids"] == ["a"]
        dock.shapes.add_polygons([[[15, 5], [15, 35], [25, 35], [25, 5]]])
        dock.shapes.selected_data = {0}
        dock.shapes.scale, dock.shapes.translate = (2, 3), (1, -2)
        captured = adapter.capture_polygon(dock.shapes)
        assert captured["polygon"][0] == [13, 31]
        dock.points.translate = (10, 20)
        with pytest.raises(SpatialError, match="transformed"):
            controller.publish_points(adapter.capture_points(dock.points))
        dock.opened(adapter.NapariController(controller.project.root))
        assert len(viewer.layers) == 3
        assert np.array_equal(adapter._world_affine(dock.points), np.eye(3))
        # Worker signals marshal result back to the Qt event loop.
        completed = []
        dock.run_job(lambda: {"worker": "finished"}, completed.append)
        import time
        deadline = time.monotonic() + 15
        while dock.worker is not None and time.monotonic() < deadline:
            QApplication.processEvents()
            time.sleep(0.01)
        assert completed == [{"worker": "finished"}] and dock.worker is None
        agent = Project(controller.project.root)
        selection = agent.set_selection(dock.controller.revision_id, cell_ids=["b"])
        proposal = agent.propose_revision(dock.controller.revision_id, selection["selection_id"],
                                         {"label": "A"}, "shared proposal")
        agent.set_selection(dock.controller.revision_id, cell_ids=["a"])
        dock.shared_proposal.setText(proposal["proposal_id"])
        dock.confirm.setChecked(True)
        dock.load_proposal()
        deadline = time.monotonic() + 15
        while dock.worker is not None and time.monotonic() < deadline:
            QApplication.processEvents()
            time.sleep(0.01)
        assert not dock.confirm.isChecked()
        assert dock.points.selected_data == {1}
        assert dock.controller.preview["proposal_id"] == proposal["proposal_id"]
        assert '"cell_ids": [\n      "b"' in dock.output.toPlainText()
        before_commit = dock.controller.revision_id
        dock.points.data[1, 0] += 1
        dock.reviewer.setText("Human fixture")
        dock.confirm.setChecked(True)
        dock.commit()
        deadline = time.monotonic() + 15
        while dock.worker is not None and time.monotonic() < deadline:
            QApplication.processEvents()
            time.sleep(0.01)
        assert agent.summary()["head_revision"] == before_commit
        assert not dock.confirm.isChecked()
        assert "moved" in dock.status.text()
        dock.opened(adapter.NapariController(controller.project.root))
        dock.load_proposal()
        deadline = time.monotonic() + 15
        while dock.worker is not None and time.monotonic() < deadline:
            QApplication.processEvents()
            time.sleep(0.01)
        dock.confirm.setChecked(True)
        dock.commit()
        deadline = time.monotonic() + 15
        while dock.worker is not None and time.monotonic() < deadline:
            QApplication.processEvents()
            time.sleep(0.01)
        final_revision = agent.get_revision(agent.summary()["head_revision"])
        assert set(final_revision["overlays"]) == {"b"}
    finally:
        viewer.close()
        second_viewer.close()
        manager.unregister(manifest.name)
