import base64

import numpy as np
import pytest

from spatial_collab.assets import _encode_png_rgba, get_asset_raster, register_asset
from spatial_collab.server import ToolService
from spatial_collab.store import Project


@pytest.fixture
def sample_project(tmp_path):
    cells = [
        {"cell_id": f"cell_{i}", "x": float(i * 10), "y": float(i * 10), "label": "Neuron", "counts": {"G1": 5.0}}
        for i in range(10)
    ]
    return Project.create(
        tmp_path / "asset_proj",
        cells,
        {
            "name": "asset raster test",
            "sample_id": "s_test",
            "slice_id": "slice_asset",
            "source_kind": "synthetic",
            "units": "micrometer",
            "coordinate_system": "cartesian",
            "panel_genes": ["G1"],
        },
    )


def test_encode_png_rgba():
    rgba = np.zeros((4, 4, 4), dtype=np.uint8)
    rgba[:, :, 0] = 255  # Red channel
    rgba[:, :, 3] = 255  # Alpha

    png_bytes = _encode_png_rgba(rgba)
    assert png_bytes.startswith(b"\x89PNG\r\n\x1a\n")
    assert b"IHDR" in png_bytes
    assert b"IDAT" in png_bytes
    assert png_bytes.endswith(b"IEND\xaeB`\x82")


def test_get_asset_raster_image(sample_project, tmp_path):
    img_data = np.linspace(0, 100, 100 * 100, dtype=np.float32).reshape(100, 100)
    npy_path = tmp_path / "img.npy"
    np.save(npy_path, img_data)

    affine = [[0.5, 0.0, 10.0], [0.0, 0.5, 20.0], [0.0, 0.0, 1.0]]
    meta = sample_project.summary()["metadata"]
    asset = register_asset(
        sample_project,
        npy_path,
        kind="image",
        source_sha256=sample_project.summary()["source_sha256"],
        slice_id=meta["slice_id"],
        coordinate_system=meta["coordinate_system"],
        units=meta["units"],
        pixel_to_world=affine,
        registration_note="Test continuous image raster",
        name="TestImage",
        sample_id=meta["sample_id"],
    )

    # Downsample with max_dimension = 50 -> step = 2
    raster = get_asset_raster(sample_project, asset["asset_id"], max_dimension=50)

    assert raster["asset_id"] == asset["asset_id"]
    assert raster["step"] == 2
    assert raster["width"] == 50
    assert raster["height"] == 50
    assert raster["data_url"].startswith("data:image/png;base64,")

    # Verify PNG bytes decode properly
    b64_payload = raster["data_url"].split(",", 1)[1]
    raw_png = base64.b64decode(b64_payload)
    assert raw_png.startswith(b"\x89PNG\r\n\x1a\n")

    # Verify scaled affine transform incorporates step: [0.5 * 2, 0, 10] = [1.0, 0, 10]
    p2w = raster["pixel_to_world"]
    assert p2w[0][0] == pytest.approx(1.0)
    assert p2w[1][1] == pytest.approx(1.0)
    assert p2w[0][2] == pytest.approx(10.0)
    assert p2w[1][2] == pytest.approx(20.0)
    assert len(raster["world_bounds"]) == 4


def test_get_asset_raster_segmentation(sample_project, tmp_path):
    mask_data = np.zeros((60, 60), dtype=np.uint32)
    mask_data[10:30, 10:30] = 1
    mask_data[35:55, 35:55] = 2

    npy_path = tmp_path / "mask.npy"
    np.save(npy_path, mask_data)

    affine = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    meta = sample_project.summary()["metadata"]
    asset = register_asset(
        sample_project,
        npy_path,
        kind="segmentation",
        source_sha256=sample_project.summary()["source_sha256"],
        slice_id=meta["slice_id"],
        coordinate_system=meta["coordinate_system"],
        units=meta["units"],
        pixel_to_world=affine,
        registration_note="Test segmentation mask raster",
        name="TestMask",
        sample_id=meta["sample_id"],
    )

    raster = get_asset_raster(sample_project, asset["asset_id"], max_dimension=100)
    assert raster["kind"] == "segmentation"
    assert raster["step"] == 1
    assert raster["width"] == 60
    assert raster["height"] == 60
    assert raster["data_url"].startswith("data:image/png;base64,")


def test_get_asset_raster_tool_service(sample_project, tmp_path):
    img_data = np.ones((20, 20), dtype=np.float32)
    npy_path = tmp_path / "img_tool.npy"
    np.save(npy_path, img_data)

    affine = [[1.0, 0.0, 5.0], [0.0, 1.0, 5.0], [0.0, 0.0, 1.0]]
    meta = sample_project.summary()["metadata"]
    asset = register_asset(
        sample_project,
        npy_path,
        kind="image",
        source_sha256=sample_project.summary()["source_sha256"],
        slice_id=meta["slice_id"],
        coordinate_system=meta["coordinate_system"],
        units=meta["units"],
        pixel_to_world=affine,
        registration_note="Test tool service dispatch",
        name="ToolImage",
        sample_id=meta["sample_id"],
    )

    service = ToolService(sample_project.root)
    res = service.get_asset_raster(asset_id=asset["asset_id"])
    assert res["asset_id"] == asset["asset_id"]
    assert res["data_url"].startswith("data:image/png;base64,")
