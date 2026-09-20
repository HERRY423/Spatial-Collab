"""Exact source/frame/mask identity and bounded read-only asset contracts."""
from copy import deepcopy
import importlib.util
import json

import numpy as np
import pytest

from spatial_collab import assets
from spatial_collab.store import Project, SpatialError


@pytest.fixture
def project(tmp_path):
    cells = [{"cell_id": "cell_z", "x": 10, "y": 20, "label": "working", "counts": {"G": 1}},
             {"cell_id": "cell_a", "x": 14, "y": 23, "label": "working", "counts": {}}]
    return Project.create(tmp_path / "project", cells,
                          {"name": "Asset fixture", "slice_id": "declared-tissue", "coordinate_system": "slice-xy",
                           "units": "micrometer", "panel_genes": ["G"], "source_kind": "synthetic"})


def declaration(project, **overrides):
    summary = project.summary()
    result = {"kind": "segmentation", "source_sha256": summary["source_sha256"],
              "slice_id": "declared-tissue", "coordinate_system": "slice-xy", "units": "micrometer",
              "pixel_to_world": [[2, 0, 10], [0, 3, 20], [0, 0, 1]],
              "registration_note": "Synthetic declared transform; no tissue identity inference",
              "label_to_cell_id": {"91": "cell_z", "7": "cell_a"}}
    result.update(overrides)
    return result


def source(tmp_path):
    path = tmp_path / "mask.npy"
    np.save(path, np.array([[91, 0, 0], [0, 0, 7]], dtype=np.uint32))
    return path


def test_explicit_labels_not_array_order_and_xy_yx_transform(project, tmp_path):
    before = project.context()
    record = assets.register_asset(project, source(tmp_path), **declaration(project))
    assert assets.register_asset(project, tmp_path / "mask.npy", **declaration(project)) == record
    result = assets.inspect_asset(project, record["asset_id"], ["cell_z", "cell_a"])
    assert [(c["label_id"], c["mapped_cell_id"], c["nearest_pixel_xy"]) for c in result["observations"]] == [
        (91, "cell_z", [0, 0]), (7, "cell_a", [2, 1])]
    assert all(c["correspondence"] == "same_declared_cell" for c in result["observations"])
    assert result["merge_or_true_coexpression"] == "NOT_DETERMINED"
    data, kwargs, kind = assets.load_asset_layer(project, record["asset_id"])
    assert not data.flags.writeable and kind == "labels"
    assert kwargs["affine"].tolist() == [[3, 0, 20], [0, 2, 10], [0, 0, 1]]
    assert project.context() == before
    listing = assets.list_assets(project)
    assert listing["segmentation_available"] and not listing["image_available"]
    assert listing["assets"][0]["mapping_status"] == "complete"


@pytest.mark.parametrize("mapping,status", [(None, "missing"), ({"91": "cell_z"}, "partial")])
def test_missing_mapping_stays_unknown(project, tmp_path, mapping, status):
    record = assets.register_asset(project, source(tmp_path), **declaration(project, label_to_cell_id=mapping))
    assert record["mapping_status"] == status
    result = assets.inspect_asset(project, record["asset_id"], ["cell_a"])
    assert result["observations"][0]["correspondence"] == "mapping_unknown"
    assert result["scientific_authorization"] == "NOT_ESTABLISHED"


@pytest.mark.parametrize("field,value", [("source_sha256", "0" * 64), ("slice_id", "other-tissue"),
                                        ("coordinate_system", "other-frame"), ("units", "pixel")])
def test_mismatched_source_frame_rejected_before_registration(project, tmp_path, field, value):
    with pytest.raises(SpatialError, match="exactly match"):
        assets.register_asset(project, source(tmp_path), **declaration(project, **{field: value}))
    assert assets.list_assets(project)["assets"] == []


@pytest.mark.parametrize("mapping", [{"0": "cell_z"}, {"01": "cell_z"}, {"91": "invented"},
                                     {"8": "cell_z"}, {91: "cell_z", "91": "cell_a"}])
def test_ambiguous_or_fabricated_mapping_rejected(project, tmp_path, mapping):
    with pytest.raises(SpatialError):
        assets.register_asset(project, source(tmp_path), **declaration(project, label_to_cell_id=mapping))


@pytest.mark.parametrize("matrix", [[[1, 0, 0], [0, 0, 0], [0, 0, 1]],
                                    [[1, 0, 0], [0, 1, 0], [1, 0, 1]],
                                    [[float("inf"), 0, 0], [0, 1, 0], [0, 0, 1]],
                                    [[1e308, 0, 0], [0, 1e308, 0], [0, 0, 1]]])
def test_unusable_transform_rejected(project, tmp_path, matrix):
    with pytest.raises(SpatialError, match="affine"):
        assets.register_asset(project, source(tmp_path), **declaration(project, pixel_to_world=matrix))


def test_resource_shape_dtype_and_lazy_budget_before_loading(project, tmp_path, monkeypatch):
    path = source(tmp_path)
    monkeypatch.setattr(assets, "MAX_PIXELS", 5)
    with pytest.raises(SpatialError, match="budget"):
        assets.register_asset(project, path, **declaration(project))
    class Lazy:
        def __array__(self):
            raise AssertionError("must never densify")
    with pytest.raises(SpatialError, match="materialized"):
        assets.array_digest(Lazy())
    with pytest.raises(SpatialError, match="two-dimensional"):
        assets.validate_array_budget((100, 100, 100), np.dtype("float64"))
    with pytest.raises(SpatialError, match="numeric"):
        assets.validate_array_budget((2, 2), np.dtype("object"))


def test_source_record_and_visible_layer_tamper_are_rejected(project, tmp_path):
    path = source(tmp_path)
    record = assets.register_asset(project, path, **declaration(project))
    data, kwargs, _ = assets.load_asset_layer(project, record["asset_id"])
    capture = {"data": data, "metadata": kwargs["metadata"], "transform": kwargs["affine"]}
    assets.validate_layer_capture(project, record["asset_id"], capture)
    moved = deepcopy(capture)
    moved["transform"][0, 2] += 1
    with pytest.raises(SpatialError, match="transformed"):
        assets.validate_layer_capture(project, record["asset_id"], moved)
    changed = deepcopy(capture)
    changed["data"][0, 0] = 7
    with pytest.raises(SpatialError, match="pixels changed"):
        assets.validate_layer_capture(project, record["asset_id"], changed)
    np.save(path, np.zeros_like(data))
    with pytest.raises(SpatialError, match="file changed"):
        assets.inspect_asset(project, record["asset_id"])
    record_path = project.root / "assets" / (record["asset_id"] + ".json")
    raw = json.loads(record_path.read_text(encoding="utf-8"))
    raw["units"] = "pixel"
    record_path.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(SpatialError, match="integrity"):
        assets.list_assets(project)


def test_origin_hashes_and_array_snapshots(project, tmp_path):
    origin = source(tmp_path)
    digest = assets.file_digest(origin)
    data = np.load(origin)
    record = assets.register_array_snapshot(project, data, origin_sources=[digest], **declaration(project))
    assert record["origin_sources"] == [digest]
    assert record["source"]["path"] != str(origin)
    data[:] = 0
    result = assets.inspect_asset(project, record["asset_id"], ["cell_z"])
    assert result["observations"][0]["label_id"] == 91
    origin.write_bytes(b"changed origin")
    with pytest.raises(SpatialError, match="origin file changed"):
        assets.inspect_asset(project, record["asset_id"])


def test_declared_sample_binding_is_required_and_exact(project, tmp_path):
    meta = deepcopy(project.summary()["metadata"])
    meta["sample_id"] = "donor-1"
    sample_project = Project.create(tmp_path / "sample-project", project.cells(), meta)
    kwargs = declaration(sample_project)
    path = source(tmp_path)
    with pytest.raises(SpatialError, match="exactly match"):
        assets.register_asset(sample_project, path, **kwargs)
    with pytest.raises(SpatialError, match="exactly match"):
        assets.register_asset(sample_project, path, sample_id="donor-2", **kwargs)
    record = assets.register_asset(sample_project, path, sample_id="donor-1", **kwargs)
    assert record["sample_id"] == "donor-1"
    assert assets.inspect_asset(sample_project, record["asset_id"])["content_integrity"] == "verified"


@pytest.mark.skipif(importlib.util.find_spec("tifffile") is None, reason="optional TIFF runtime unavailable")
def test_tiff_plane_is_explicit_and_no_rgb_or_implicit_z_projection(project, tmp_path):
    import tifffile
    path = tmp_path / "image.tif"
    with tifffile.TiffWriter(path) as writer:
        writer.write(np.ones((2, 3), dtype=np.uint16))
        writer.write(np.full((2, 3), 20, dtype=np.uint16))
    kwargs = declaration(project, kind="image", label_to_cell_id=None)
    with pytest.raises(SpatialError, match="explicit tiff_page"):
        assets.register_asset(project, path, **kwargs)
    record = assets.register_asset(project, path, tiff_page=1, **kwargs)
    assert assets.inspect_asset(project, record["asset_id"], ["cell_z"])["observations"][0]["pixel_value"] == 20
    assert record["mapping_status"] == "not_applicable"
    with pytest.raises(SpatialError, match="TIFF page"):
        assets.register_asset(project, path, tiff_page=2, **kwargs)
