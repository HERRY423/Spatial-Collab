"""Real SpatialData serialization: exact joins, physical transforms and exclusions."""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import pytest
from scipy import sparse

from spatial_collab.spatialdata_adapter import import_spatialdata
from spatial_collab.store import SpatialError

spatialdata = pytest.importorskip("spatialdata")
ad = pytest.importorskip("anndata")
gpd = pytest.importorskip("geopandas")
box = pytest.importorskip("shapely.geometry").box
models = pytest.importorskip("spatialdata.models")
transforms = pytest.importorskip("spatialdata.transformations")
PointsModel, ShapesModel, TableModel = models.PointsModel, models.ShapesModel, models.TableModel
Affine, Identity = transforms.Affine, transforms.Identity


def _files(root: Path):
    return {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in root.rglob("*") if path.is_file()}


def _write_store(tmp_path, *, category="points", instance_ids=(33, 11, 22), three_dimensional=False,
                 transform=None, sample_ids=None, duplicate_geometry=False):
    coordinate_frame = "calibrated_um"
    if transform is None:
        transform = Affine([[2, 0, 100], [0, 3, 200], [0, 0, 1]],
                           input_axes=("x", "y"), output_axes=("x", "y"))
    geometry_ids = [11, 22, 22 if duplicate_geometry else 33]
    if category == "points":
        coordinates = pd.DataFrame({"x": [1., 4., 8.], "y": [2., 6., 10.]}, index=geometry_ids)
        if three_dimensional:
            coordinates["z"] = [0., 1., 2.]
        element = PointsModel.parse(coordinates, transformations={coordinate_frame: transform})
    else:
        element = ShapesModel.parse(gpd.GeoDataFrame(
            {"geometry": [box(0, 1, 2, 3), box(3, 5, 5, 7), box(7, 9, 9, 11)]},
            index=geometry_ids), transformations={coordinate_frame: transform})
    obs = pd.DataFrame({"region": pd.Categorical(["geometry"] * 3), "instance_id": list(instance_ids),
                        "cell_type": ["Third", "First", "Second"]},
                       index=["observation-third", "observation-first", "observation-second"])
    if sample_ids is not None:
        obs["sample_id"] = sample_ids
    table = ad.AnnData(X=sparse.csr_matrix([[30., 3.], [10., 1.], [20., 2.]]),
                      obs=obs, var=pd.DataFrame(index=["Marker", "Control"]))
    table.layers["counts"] = table.X.copy()
    table = TableModel.parse(table, region="geometry", region_key="region", instance_key="instance_id")
    data = spatialdata.SpatialData(tables={"expression": table}, **{category: {"geometry": element}})
    path = tmp_path / "input.zarr"
    data.write(path)
    return path


def _import(path, destination, **overrides):
    options = {"table_name": "expression", "element_name": "geometry", "coordinate_system": "calibrated_um",
               "units": "micrometer", "slice_id": "section-one", "observation_unit": "cell"}
    options.update(overrides)
    return import_spatialdata(path, destination, **options)


@pytest.mark.parametrize("category", ["points", "shapes"])
def test_public_spatialdata_roundtrip_joins_shuffled_ids_after_nonidentity_transform(tmp_path, category):
    path = _write_store(tmp_path, category=category)
    original = _files(path)
    project = _import(path, tmp_path / "review")
    rows = {row["cell_id"]: row for row in project.cells()}
    for identifier, expected in {
        "observation-first": (102., 206., 10., "First"),
        "observation-second": (108., 218., 20., "Second"),
        "observation-third": (116., 230., 30., "Third"),
    }.items():
        x, y, marker, label = expected
        assert rows[identifier]["x"] == pytest.approx(x)
        assert rows[identifier]["y"] == pytest.approx(y)
        assert rows[identifier]["counts"]["Marker"] == marker
        assert rows[identifier]["label"] == label
    meta = project.summary()["metadata"]
    assert meta["spatialdata"]["instance_ids"] == ["33", "11", "22"]
    assert meta["spatialdata"]["element_type"] == category
    assert meta["units"] == "micrometer"
    assert meta["coordinate_system"] == "calibrated_um"
    assert meta["spatialdata"]["transform_xy"] == [[2., 0., 100.], [0., 3., 200.], [0., 0., 1.]]
    assert _files(path) == original


def test_wrong_geometry_join_rejects_instead_of_positionally_assigning(tmp_path):
    path = _write_store(tmp_path, instance_ids=(99, 11, 22))
    with pytest.raises(SpatialError, match="instance IDs and geometry IDs differ"):
        _import(path, tmp_path / "review")
    assert not (tmp_path / "review" / "project.sqlite3").exists()


def test_missing_coordinate_frame_and_uncalibrated_units_never_fall_back(tmp_path):
    path = _write_store(tmp_path)
    with pytest.raises((SpatialError, ValueError), match="(?i)coordinate|transformation|system"):
        _import(path, tmp_path / "missing", coordinate_system="absent_frame")
    with pytest.raises(SpatialError, match="calibration"):
        _import(path, tmp_path / "pixels", units="pixel")
    assert not (tmp_path / "missing" / "project.sqlite3").exists()
    assert not (tmp_path / "pixels" / "project.sqlite3").exists()


def test_three_dimensions_cannot_be_projected_implicitly(tmp_path):
    path = _write_store(tmp_path, three_dimensional=True, transform=Identity())
    with pytest.raises(SpatialError, match="3D"):
        _import(path, tmp_path / "review")


def test_singular_transform_cannot_create_false_distances(tmp_path):
    path = _write_store(tmp_path, transform=Affine([[0, 0, 100], [0, 3, 200], [0, 0, 1]],
                                                 input_axes=("x", "y"), output_axes=("x", "y")))
    with pytest.raises(SpatialError, match="invertible"):
        _import(path, tmp_path / "review")


def test_raw_layer_missing_and_mixed_samples_are_explicit_rejections(tmp_path):
    path = _write_store(tmp_path)
    with pytest.raises(SpatialError, match="raw-count layer missing"):
        _import(path, tmp_path / "wrong-layer", counts_layer="not_counts")
    mixed_dir = tmp_path / "mixed"
    mixed_dir.mkdir()
    mixed = _write_store(mixed_dir, sample_ids=["donor-a", "donor-a", "donor-b"])
    with pytest.raises(SpatialError, match="mixed or missing sample"):
        _import(mixed, tmp_path / "mixed-project")


def test_source_store_cannot_be_project_destination(tmp_path):
    path = _write_store(tmp_path)
    original = _files(path)
    with pytest.raises(SpatialError, match="outside the immutable source"):
        _import(path, path / "project")
    assert _files(path) == original


def test_duplicate_geometry_ids_are_rejected(tmp_path):
    path = _write_store(tmp_path, duplicate_geometry=True)
    with pytest.raises(SpatialError, match="duplicate instance IDs"):
        _import(path, tmp_path / "review")


def test_geometry_budget_rejects_before_materializing_public_reader(tmp_path, monkeypatch):
    from spatial_collab import spatialdata_adapter
    path = _write_store(tmp_path)
    called = []

    def must_not_read(*args, **kwargs):
        called.append(True)
        raise AssertionError("Full SpatialData reader must not run before resource checks pass")

    monkeypatch.setattr(spatialdata_adapter, "MAX_CELLS", 2)
    monkeypatch.setattr(spatialdata, "read_zarr", must_not_read)
    with pytest.raises(SpatialError, match="geometry category"):
        _import(path, tmp_path / "review")
    assert not called
