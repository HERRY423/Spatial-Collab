import csv
from io import BytesIO
import base64

import numpy as np
import pytest

from spatial_collab import atlas, pyramid
from spatial_collab.demo import create_demo
from spatial_collab.store import SpatialError
from spatial_collab.workflow_inputs import load_counts


def write_csv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def atlas_spec(tmp_path, n=120):
    write_csv(tmp_path / "cells.csv", ["id", "x", "y"], [{"id": str(i), "x": i % 20, "y": i // 20} for i in range(n)])
    write_csv(tmp_path / "features.csv", ["id", "symbol"], [{"id": "G1", "symbol": "One"}, {"id": "G2", "symbol": "Two"}, {"id": "zero", "symbol": "Zero"}])
    write_csv(tmp_path / "counts.csv", ["cell_id", "feature", "value"], [{"cell_id": str(i), "feature": "G1", "value": i+1} for i in range(n)])
    write_csv(tmp_path / "transcripts.csv", ["id", "feature", "x", "y", "qv"], [{"id": str(i), "feature": "G1", "x": i % 20 + .1, "y": i // 20 + .1, "qv": 5 if i % 2 else 30} for i in range(n)])
    write_csv(tmp_path / "boundaries.csv", ["id", "cell_id", "kind", "vertex", "x", "y"], [{"id": "b0", "cell_id": "0", "kind": "cell", "vertex": i, "x": x, "y": y} for i, (x, y) in enumerate([[-.2, -.2], [.2, -.2], [.2, .2], [-.2, .2]])])
    spec = {"metadata": {"name": "Synthetic index control", "sample_id": "test", "coordinate_system": "measured_xy", "units": "micrometer", "platform": "synthetic", "observation_unit": "cell", "species": "mouse"}}
    for kind in ("cells", "features", "counts", "transcripts", "boundaries"):
        path = tmp_path / f"{kind}.csv"
        with path.open() as stream:
            fields = next(csv.reader(stream))
        spec[kind] = {"path": str(path), "columns": {f: f for f in fields}}
    spec["transcripts"]["projection"] = "xy_projection"
    return spec


def test_empty_workspace_import_view_and_replay(tmp_path):
    from spatial_collab.store import Project
    from spatial_collab.server import ToolService
    from spatial_collab.replay import verify_bundle
    p = Project.create_workspace(tmp_path / "workspace")
    assert p.cells() == []
    service = ToolService(p.root)
    assert service.open_project()["cell_count"] == 0
    assert service.get_overview()["observation_count"] == 0
    record = atlas.build(p, atlas_spec(tmp_path, 6))
    assert atlas.view(p, record["object_id"])["total"] == 6
    assert atlas.extract(p, record["object_id"], [-1, -1, 10, 1])["shape"] == [6, 3]
    verify_bundle(p.export_bundle(compact=True)["export_path"])
    with pytest.raises(SpatialError, match="already exists"):
        Project.create_workspace(p.root)
    metadata = p.summary()["metadata"]
    metadata["source_kind"] = "ordinary_import"
    with pytest.raises(SpatialError, match="nonempty"):
        Project.create(tmp_path / "bad", [], metadata)


@pytest.fixture
def data(tmp_path):
    p = create_demo(tmp_path / "project")
    spec = atlas_spec(tmp_path)
    record = atlas.build(p, spec)
    return p, record, spec


def test_density_counts_every_object_without_sampling(data):
    p, record, _ = data
    result = atlas.view(p, record["object_id"], limit=10, grid=16)
    assert result["mode"] == "aggregate_density"
    assert result["aggregate_total"] == result["total"] == 120
    cropped = atlas.view(p, record["object_id"], [.5, .5, 4.5, 4.5], limit=1)
    assert cropped["aggregate_total"] == 16
    assert all("cell_id" not in r for r in result["records"])
    assert atlas.view(p, record["object_id"], layer="transcripts", limit=1)["aggregate_total"] == 120


def test_exact_roi_molecules_boundaries_and_sparse_counts(data):
    p, record, _ = data
    oid = record["object_id"]
    result = atlas.view(p, oid, [-.5, -.5, 1.5, .5])
    assert [r["source_id"] for r in result["records"]] == ["0", "1"]
    molecules = atlas.view(p, oid, [-.5, -.5, 1.5, .5], layer="transcripts", feature="G1")
    assert [r["qv"] for r in molecules["records"]] == [30, 5]
    boundary = atlas.view(p, oid, [-.1, -.1, .1, .1], layer="boundaries")
    assert boundary["records"][0]["cell_id"] == "0"
    frozen = atlas.extract(p, oid, [-.5, -.5, 1.5, .5])
    r, x = load_counts(p, frozen["object_id"])
    assert r["features"] == ["G1", "G2", "zero"]
    np.testing.assert_array_equal(x.toarray(), [[1, 0, 0], [2, 0, 0]])
    with pytest.raises(SpatialError, match="No sampling"):
        atlas.extract(p, oid, record["bounds"], max_cells=3)


def test_atlas_detects_modified_storage(data):
    p, r, _ = data
    path = p.root / r["database"]
    with path.open("ab") as stream:
        stream.write(b"tamper")
    with pytest.raises(SpatialError, match="changed"):
        atlas.view(p, r["object_id"])


def test_shared_atlas_view_is_distinct_from_primary_selection(data):
    p, r, _ = data
    before = p.context()
    shared = atlas.share_view(p, r["object_id"], [-.5, -.5, 1.5, .5], layer="transcripts", feature="G1")
    after = p.context()
    assert after["atlas_view"] == shared
    assert before["cursor"] != after["cursor"]
    assert before["head_revision"] == after["head_revision"]
    assert before["selection_id"] == after["selection_id"]


@pytest.mark.parametrize("value", [-1, .2, "nan"])
def test_stream_import_rejects_nonraw_counts(tmp_path, value):
    p = create_demo(tmp_path / "project")
    spec = atlas_spec(tmp_path)
    write_csv(tmp_path / "counts.csv", ["cell_id", "feature", "value"], [{"cell_id": "0", "feature": "G1", "value": value}])
    with pytest.raises(SpatialError):
        atlas.build(p, spec)
    assert not list((p.root / "atlases").glob("*.sqlite3"))


def test_native_image_pyramid_reads_exact_pixels_and_affine(data, tmp_path):
    tifffile = pytest.importorskip("tifffile")
    pytest.importorskip("zarr")
    from PIL import Image
    p, atlas_record, _ = data
    pixels = np.arange(512*512, dtype=np.uint8).reshape(512, 512)
    path = tmp_path / "pyramid.ome.tif"
    with tifffile.TiffWriter(path, ome=True) as writer:
        writer.write(pixels, tile=(128, 128), subifds=1, metadata={"axes": "YX"})
        writer.write(pixels[::2, ::2], tile=(128, 128), subfiletype=1)
    r = pyramid.register(p, path, atlas_id=atlas_record["object_id"], pixel_to_world=[[2, 0, 10], [0, 3, 20], [0, 0, 1]])
    assert len(r["levels"]) == 2
    full_tile = pyramid.tile(p, r["object_id"], 0, 1, 1)
    actual_full = np.asarray(Image.open(BytesIO(base64.b64decode(full_tile["data_url"].split(",")[1]))))
    np.testing.assert_array_equal(actual_full, pixels[256:, 256:])
    tile = pyramid.tile(p, r["object_id"], 1, 0, 0)
    actual = np.asarray(Image.open(BytesIO(base64.b64decode(tile["data_url"].split(",")[1]))))
    np.testing.assert_array_equal(actual, pixels[::2, ::2])
    assert tile["world_corners"] == [[9.0, 18.5], [1033.0, 18.5], [9.0, 1554.5]]
    with pytest.raises(SpatialError):
        pyramid.tile(p, r["object_id"], 1, 1, 0)


def test_h5ad_streaming_count_axis_is_producer_identity(tmp_path):
    ad = pytest.importorskip("anndata")
    from scipy import sparse
    p = create_demo(tmp_path / "project")
    spec = atlas_spec(tmp_path, 6)
    a = ad.AnnData(sparse.csr_matrix(np.arange(18).reshape(6, 3)))
    a.obs_names = [str(i) for i in range(6)]
    a.var_names = ["wrong1", "wrong2", "wrong3"]
    a.var["gene_id"] = ["G1", "G2", "zero"]
    path = tmp_path / "raw.h5ad"
    a.write(path)
    spec["counts"] = {"path": str(path), "format": "h5ad_csr", "counts_layer": "X", "feature_field": "gene_id"}
    r = atlas.build(p, spec)
    frozen = atlas.extract(p, r["object_id"], [-1, -1, 10, 1])
    _, actual = load_counts(p, frozen["object_id"])
    np.testing.assert_array_equal(actual.toarray(), np.arange(18).reshape(6, 3))
