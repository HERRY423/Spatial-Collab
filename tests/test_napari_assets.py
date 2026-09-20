"""Native hidden-window assets test; this is not visible GUI or biological acceptance."""
import importlib.util
import time

import numpy as np
import pytest

from spatial_collab import assets
from spatial_collab.napari_adapter import NapariController, capture_asset_layer, capture_points, make_dock_widget
from spatial_collab.store import Project, SpatialError


def create_project(path, units="micrometer"):
    return Project.create(path, [{"cell_id": "z", "x": 10, "y": 20, "label": "A", "counts": {"G": 1}},
                                 {"cell_id": "a", "x": 14, "y": 23, "label": "B", "counts": {}}],
                          {"name": "asset host fixture", "slice_id": "s", "coordinate_system": "s-xy",
                           "units": units, "panel_genes": ["G"], "source_kind": "synthetic"})


@pytest.mark.parametrize("units", ["pixel", "array_index", "unknown"])
def test_nonphysical_frame_retained_and_physical_analysis_refused(tmp_path, units):
    project = create_project(tmp_path / "project", units)
    controller = NapariController(project.root)
    assert controller.binding["units"] == units
    with pytest.raises(SpatialError, match="micrometer"):
        controller.compare(controller.revision_id, controller.revision_id,
                           radius_um=10, source_label="A", target_label="B")


@pytest.mark.skipif(importlib.util.find_spec("napari") is None, reason="optional native napari runtime unavailable")
def test_native_asset_load_registration_inspection_and_movement_rejection(tmp_path):
    import napari
    from qtpy.QtWidgets import QApplication
    project = create_project(tmp_path / "project")
    controller = NapariController(project.root)
    before = project.context()
    path = tmp_path / "mask.npy"
    np.save(path, np.array([[91, 0, 0], [0, 0, 7]], dtype=np.uint32))
    record = assets.register_asset(project, path, kind="segmentation", source_sha256=controller.binding["source_sha256"],
                                    slice_id="s", coordinate_system="s-xy", units="micrometer",
                                    pixel_to_world=[[2, 0, 10], [0, 3, 20], [0, 0, 1]],
                                    registration_note="Explicit synthetic format check", label_to_cell_id={"91": "z", "7": "a"})
    viewer = napari.Viewer(show=False)
    def wait(dock):
        deadline = time.monotonic() + 15
        while dock.worker is not None and time.monotonic() < deadline:
            QApplication.processEvents()
            time.sleep(.01)
        assert dock.worker is None
    try:
        unrelated = viewer.add_image(np.ones((3, 3)), name="unrelated image")
        dock = make_dock_widget(viewer)
        viewer.window.add_dock_widget(dock)
        dock.opened(controller)
        dock.refresh_assets()
        wait(dock)
        assert dock.asset_choice.currentData() == record["asset_id"]
        dock.load_asset()
        wait(dock)
        assert record["asset_id"] in dock.asset_layers, dock.status.text()
        layer = dock.asset_layers[record["asset_id"]]
        assert not layer.editable
        assert layer.data.flags.writeable  # native renderer needs its own writable buffer
        assert np.allclose(layer.data_to_world([1, 2]), [23, 14])
        dock.points.selected_data = {0, 1}
        result = controller.inspect_asset(record["asset_id"], capture_asset_layer(layer), capture_points(dock.points))
        assert [(i["cell_id"], i["label_id"]) for i in result["observations"]] == [("a", 7), ("z", 91)]
        assert all(i["correspondence"] == "same_declared_cell" for i in result["observations"])
        layer.translate = (3, 0)
        with pytest.raises(SpatialError, match="transformed"):
            controller.inspect_asset(record["asset_id"], capture_asset_layer(layer), capture_points(dock.points))
        layer.translate = (0, 0)
        layer.data = np.zeros_like(layer.data)
        with pytest.raises(SpatialError, match="pixels changed"):
            controller.inspect_asset(record["asset_id"], capture_asset_layer(layer), capture_points(dock.points))
        # Existing native image -> explicit snapshot, retaining world transform.
        existing = viewer.add_image(np.full((2, 3), 17, dtype=np.uint16), scale=(3, 2), translate=(20, 10))
        snapshot = controller.register_layer_snapshot(capture_asset_layer(existing),
                    source_sha256=controller.binding["source_sha256"], slice_id="s", coordinate_system="s-xy",
                    units="micrometer", registration_note="Synthetic existing image with explicit declared frame",
                    origin_sources=[assets.file_digest(path)], name="existing image snapshot")
        assert assets.inspect_asset(project, snapshot["asset_id"], ["a"])["observations"][0]["pixel_value"] == 17
        assert unrelated in viewer.layers
        assert project.context() == before
        # Uncalibrated project closes only this dock's old source-bound assets.
        pixel_project = create_project(tmp_path / "pixel-project", "pixel")
        dock.opened(NapariController(pixel_project.root))
        assert not dock.compare_button.isEnabled()
        assert layer not in viewer.layers and unrelated in viewer.layers and existing in viewer.layers
    finally:
        viewer.close()
