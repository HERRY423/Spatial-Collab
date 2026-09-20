"""Format-profile fixtures exercise identity, geometry, units and failure paths."""
from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from spatial_collab import format_readers as readers
from spatial_collab.import_registry import import_source, list_formats, probe_source, register_reader
from spatial_collab.importers import import_h5ad, import_xenium, source_digest


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return path


def matrix(root, *, ids=("a", "b"), compressed=False, data=None):
    suffix = ".gz" if compressed else ""
    write(root / ("barcodes.tsv" + suffix), "\n".join(ids) + "\n")
    write(root / ("features.tsv" + suffix), "ENSG1\tG1\tGene Expression\nENSG2\tG2\tGene Expression\ncontrol\tNeg1\tNegative Control Probe\n")
    write(root / ("matrix.mtx" + suffix), data or f"%%MatrixMarket matrix coordinate integer general\n% test\n3 {len(ids)} 3\n1 1 7\n2 2 9\n3 1 99\n")
    return root


def sidecars(root, ids=("a", "b")):
    positions = write(root / "positions.csv", "cell_id,x,y\n" + "".join(f"{identifier},{i * 2},{i * 3}\n" for i, identifier in enumerate(ids)))
    labels = write(root / "annotations.csv", "cell_id,label\n" + "".join(f"{identifier},Type {i}\n" for i, identifier in enumerate(ids)))
    return positions, labels


def generic(root, destination, positions, annotations, **options):
    return import_source("tenx_mtx", root, destination, positions=positions, annotations=annotations,
                         slice_id="slice-one", coordinate_system="native", units="micrometer", **options)


@pytest.mark.parametrize("compressed", [False, True])
def test_generic_mtx_counts_ids_and_source_hashes(tmp_path, compressed):
    root = matrix(tmp_path / "mtx", ids=("b", "a"), compressed=compressed)
    positions, labels = sidecars(tmp_path)
    before = source_digest(positions)
    project = generic(root, tmp_path / "project", positions, labels, observation_unit="bin", platform="declared assay")
    cells = {cell["cell_id"]: cell for cell in project.cells()}
    assert cells["b"]["counts"] == {"G1": 7}
    assert cells["a"]["counts"] == {"G2": 9}
    assert cells["b"]["x"] == 2 and cells["b"]["y"] == 3
    assert cells["b"]["label"] == "Type 1"
    metadata = project.summary()["metadata"]
    assert metadata["observation_unit"] == "bin" and metadata["label_semantics"] == "bin_annotation"
    assert metadata["panel_genes"] == ["G1", "G2"]
    assert len(metadata["source_files"]) == 5
    assert source_digest(positions) == before
    for record in metadata["source_files"]:
        assert source_digest(record["path"])["sha256"] == record["sha256"]


def test_annotation_free_import_requires_explicit_option(tmp_path):
    root = matrix(tmp_path / "mtx")
    positions, _ = sidecars(tmp_path)
    with pytest.raises(ValueError, match="Annotations are required"):
        generic(root, tmp_path / "strict", positions, None)
    project = generic(root, tmp_path / "exploration", positions, None, allow_unannotated=True)
    assert {cell["label"] for cell in project.cells()} == {"Unannotated"}
    assert project.summary()["metadata"]["annotation_status"] == "unannotated"


@pytest.mark.parametrize("file,content,message", [
    ("barcodes.tsv", "a\na\n", "Duplicate"),
    ("barcodes.tsv", "a\n", "dimensions"),
    ("features.tsv", "g1\tG1\tGene Expression\ng2\tG1\tGene Expression\nc\tC\tNegative\n", "Duplicate"),
    ("matrix.mtx", "%%MatrixMarket matrix coordinate integer general\n3 2 1\n1 1 1.5\n", "integer"),
    ("matrix.mtx", "%%MatrixMarket matrix coordinate integer general\n3 2 1\n1 1 -2\n", "integer"),
    ("matrix.mtx", "%%MatrixMarket matrix coordinate real general\n3 2 1\n1 1 nan\n", "integer"),
    ("matrix.mtx", "%%MatrixMarket matrix coordinate integer general\n3 2 1\n4 1 2\n", "outside"),
    ("matrix.mtx", "%%MatrixMarket matrix coordinate integer general\n3 2 2\n1 1 2\n", "truncated"),
    ("matrix.mtx", "%%MatrixMarket matrix coordinate integer symmetric\n3 2 1\n1 1 2\n", "general coordinate"),
])
def test_matrix_market_rejects_malformed_inputs(tmp_path, file, content, message):
    root = matrix(tmp_path / "mtx")
    write(root / file, content)
    positions, labels = sidecars(tmp_path)
    with pytest.raises(ValueError, match=message):
        generic(root, tmp_path / "project", positions, labels)
    assert not (tmp_path / "project").exists()


def test_matrix_duplicate_values_aggregate_without_unsigned_wrap(tmp_path):
    root = matrix(tmp_path / "mtx", data="%%MatrixMarket matrix coordinate integer general\n3 2 2\n1 1 3000000000\n1 1 3000000000\n")
    positions, labels = sidecars(tmp_path)
    project = generic(root, tmp_path / "project", positions, labels)
    assert project.cells()[0]["counts"]["G1"] == 6000000000


def test_matrix_resource_budget_before_array_allocation(tmp_path, monkeypatch):
    root = matrix(tmp_path / "mtx", data="%%MatrixMarket matrix coordinate integer general\n3 2 2000001\n")
    positions, labels = sidecars(tmp_path)
    def forbidden(*_args, **_kwargs):
        raise AssertionError("Allocation occurred before resource check")
    monkeypatch.setattr(readers.np, "empty", forbidden)
    with pytest.raises(ValueError, match="resource budget"):
        generic(root, tmp_path / "project", positions, labels)


def test_positions_and_labels_require_exact_alignment(tmp_path):
    root = matrix(tmp_path / "mtx")
    positions, labels = sidecars(tmp_path, ids=("a", "other"))
    with pytest.raises(ValueError, match="Position IDs"):
        generic(root, tmp_path / "project", positions, labels)
    positions, _ = sidecars(tmp_path / "correct")
    with pytest.raises(ValueError, match="Annotation IDs"):
        generic(root, tmp_path / "project", positions, labels, allow_unannotated=True)


def test_source_mutation_during_import_is_rejected(tmp_path, monkeypatch):
    root = matrix(tmp_path / "mtx")
    positions, labels = sidecars(tmp_path)
    original = readers._annotations
    def altered(*args, **kwargs):
        result = original(*args, **kwargs)
        labels.write_text("cell_id,label\na,Changed\nb,Changed\n")
        return result
    monkeypatch.setattr(readers, "_annotations", altered)
    with pytest.raises(ValueError, match="Source changed"):
        generic(root, tmp_path / "project", positions, labels)
    assert not (tmp_path / "project").exists()


def test_xenium_mtx_fallback_and_unannotated(tmp_path):
    outs = tmp_path / "outs"
    matrix(outs / "cell_feature_matrix")
    write(outs / "cells.csv", "cell_id,x_centroid,y_centroid\nb,20,30\na,1,2\n")
    project = import_xenium(outs, None, tmp_path / "project", slice_id="s", allow_unannotated=True)
    assert project.cells()[0]["counts"] == {"G1": 7}
    assert project.cells()[1]["x"] == 20
    assert project.summary()["metadata"]["platform"] == "Xenium"
    assert project.summary()["metadata"]["import_parameters"]["matrix_format"] == "10x_mtx"


def test_xenium_requires_feature_types_even_in_mtx_fallback(tmp_path):
    outs = tmp_path / "outs"
    matrix(outs / "cell_feature_matrix")
    write(outs / "cell_feature_matrix" / "features.tsv", "g1\tG1\ng2\tG2\ncontrol\tNeg1\n")
    write(outs / "cells.csv", "cell_id,x_centroid,y_centroid\na,1,2\nb,20,30\n")
    with pytest.raises(ValueError, match="typed three-column"):
        import_xenium(outs, None, tmp_path / "project", slice_id="s", allow_unannotated=True)


def xenium_h5_window_fixture(root):
    import h5py
    import numpy as np
    from scipy import sparse
    write(root / "cells.csv", "cell_id,x_centroid,y_centroid,transcript_counts,genomic_control_counts,cell_area,nucleus_area,segmentation_method\nb,2,0,9,0,15,3,nuclear_expansion\na,0,0,7,2,12,,cell_boundary\n")
    mat = sparse.csc_matrix([[7, 0], [0, 9], [99, 88]])
    with h5py.File(root / "cell_feature_matrix.h5", "w") as handle:
        group = handle.create_group("matrix")
        for key in ("data", "indices", "indptr"):
            group[key] = getattr(mat, key)
        group["shape"] = mat.shape
        group["barcodes"] = np.array([b"a", b"b"])
        features = group.create_group("features")
        features["name"] = np.array([b"G1", b"G2", b"Control"])
        features["feature_type"] = np.array([b"Gene Expression", b"Gene Expression", b"Negative Control Probe"])
    return root


def test_xenium_window_streams_source_without_loading_full_count_arrays(tmp_path, monkeypatch):
    import h5py
    import spatial_collab.importers as importers
    root = xenium_h5_window_fixture(tmp_path / "xenium")
    # Source has two cells but selected budget is one: full import must still fail.
    monkeypatch.setattr(importers, "MAX_CELLS", 1)
    with pytest.raises(ValueError, match="CSV exceeds"):
        import_xenium(root, None, tmp_path / "full", slice_id="s", allow_unannotated=True)
    original = h5py.Dataset.__getitem__
    reads = []
    def guarded(dataset, key):
        if dataset.name in {"/matrix/data", "/matrix/indices"}:
            reads.append((dataset.name, key))
            assert isinstance(key, slice) and key.start is not None and key.stop is not None, "Full count-array read"
            assert key.start == 0 and key.stop == 2, "Read counts outside selected CSC column"
        return original(dataset, key)
    monkeypatch.setattr(h5py.Dataset, "__getitem__", guarded)
    project = import_xenium(root, None, tmp_path / "window", slice_id="s", allow_unannotated=True,
                            bounds=[-1, -1, 1, 1])
    assert len(reads) == 2
    assert [cell["cell_id"] for cell in project.cells()] == ["a"]
    assert project.cells()[0]["counts"] == {"G1": 7}
    metadata = project.summary()["metadata"]
    assert metadata["panel_genes"] == ["G1", "G2"]
    scope = metadata["import_scope"]
    assert scope["source_observation_count"] == 2 and scope["selected_observation_count"] == 1
    assert scope["kind"] == "spatial_window" and scope["sampling"] == "none"
    attributes = project.cells()[0]["attributes"]
    assert attributes["genomic_control_counts"] == 2 and attributes["transcript_counts"] == 7
    assert attributes["nucleus_area"] is None and "total_counts" not in attributes
    assert attributes["segmentation_method"] == "cell_boundary"
    project.set_selection(project.summary()["head_revision"], cell_ids=["a"])
    assert project.inspect_selection(["G2"])["expression"]["G2"]["sum"] == 0


def test_xenium_window_closed_boundary_includes_all_centroids(tmp_path):
    root = xenium_h5_window_fixture(tmp_path / "xenium")
    project = import_xenium(root, None, tmp_path / "window", slice_id="s", allow_unannotated=True,
                            bounds=[0, -1, 2, 1])
    assert {cell["cell_id"] for cell in project.cells()} == {"a", "b"}


@pytest.mark.parametrize("bounds", [[], [1, 2, 3], [1, 1, 1, 2], [2, 1, 1, 2], [False, 0, 1, 1],
                                    [0, 0, "1", 1], [0, 0, float("inf"), 1]])
def test_xenium_window_bounds_strict(tmp_path, bounds):
    root = xenium_h5_window_fixture(tmp_path / "xenium")
    with pytest.raises(ValueError, match="bounds"):
        import_xenium(root, None, tmp_path / "window", slice_id="s", allow_unannotated=True, bounds=bounds)


def test_xenium_window_empty_and_entry_budget(tmp_path, monkeypatch):
    import spatial_collab.importers as importers
    root = xenium_h5_window_fixture(tmp_path / "xenium")
    with pytest.raises(ValueError, match="contains no cell"):
        import_xenium(root, None, tmp_path / "empty", slice_id="s", allow_unannotated=True, bounds=[10, 10, 11, 11])
    monkeypatch.setattr(importers, "MAX_NNZ", 1)
    with pytest.raises(ValueError, match="stored matrix entries"):
        import_xenium(root, None, tmp_path / "budget", slice_id="s", allow_unannotated=True, bounds=[-1, -1, 1, 1])
    assert not (tmp_path / "budget").exists()


def test_xenium_window_validates_ids_outside_selected_area(tmp_path):
    import h5py
    root = xenium_h5_window_fixture(tmp_path / "xenium")
    with h5py.File(root / "cell_feature_matrix.h5", "r+") as handle:
        handle["matrix/barcodes"][1] = b"z"
    with pytest.raises(ValueError, match="unmatched source barcodes"):
        import_xenium(root, None, tmp_path / "window", slice_id="s", allow_unannotated=True, bounds=[-1, -1, 1, 1])


def test_xenium_window_annotation_scope_must_match_imported_cells(tmp_path):
    root = xenium_h5_window_fixture(tmp_path / "xenium")
    labels = write(tmp_path / "labels.csv", "cell_id,label\na,T cell\nb,Myeloid\n")
    with pytest.raises(ValueError, match="Annotation IDs"):
        import_xenium(root, labels, tmp_path / "window", slice_id="s", bounds=[-1, -1, 1, 1])


def test_xenium_window_qc_nonfinite_rejected(tmp_path):
    root = xenium_h5_window_fixture(tmp_path / "xenium")
    write(root / "cells.csv", "cell_id,x_centroid,y_centroid,cell_area\na,0,0,inf\nb,2,0,3\n")
    with pytest.raises(ValueError, match="Nonfinite cell_area"):
        import_xenium(root, None, tmp_path / "window", slice_id="s", allow_unannotated=True, bounds=[-1, -1, 1, 1])


def test_xenium_window_optional_qc_nan_is_missing_not_zero(tmp_path):
    root = xenium_h5_window_fixture(tmp_path / "xenium")
    write(root / "cells.csv", "cell_id,x_centroid,y_centroid,nucleus_area,nucleus_count\na,0,0,NaN,0\nb,2,0,3,1\n")
    project = import_xenium(root, None, tmp_path / "window", slice_id="s", allow_unannotated=True, bounds=[-1, -1, 1, 1])
    attributes = project.cells()[0]["attributes"]
    assert attributes["nucleus_area"] is None and attributes["nucleus_count"] == 0


def visium_fixture(root, legacy=False):
    matrix(root / "raw_feature_bc_matrix")
    header = "" if legacy else "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n"
    name = "tissue_positions_list.csv" if legacy else "tissue_positions.csv"
    write(root / "spatial" / name, header + "b,0,1,2,200,300\na,1,0,0,10,30\n")
    return root


@pytest.mark.parametrize("legacy", [False, True])
def test_visium_axes_calibration_spot_semantics_and_subset(tmp_path, legacy):
    root = visium_fixture(tmp_path / "visium", legacy)
    project = import_source("visium", root, tmp_path / "project", slice_id="s", microns_per_pixel=.5,
                            in_tissue_only=True, allow_unannotated=True)
    assert len(project.cells()) == 1
    spot = project.cells()[0]
    assert spot["cell_id"] == "a" and spot["x"] == 15 and spot["y"] == 5
    assert spot["counts"] == {"G1": 7}
    metadata = project.summary()["metadata"]
    assert metadata["observation_unit"] == "spot" and metadata["label_semantics"] == "spot_annotation"
    assert metadata["import_parameters"]["excluded_off_tissue_count"] == 1
    assert "not pure cell identities" in " ".join(metadata["limitations"])


def test_visium_all_spots_and_exact_annotation_scope(tmp_path):
    root = visium_fixture(tmp_path / "visium")
    _, labels = sidecars(tmp_path)
    project = import_source("visium", root, tmp_path / "all", slice_id="s", microns_per_pixel=2,
                            in_tissue_only=False, annotations=labels)
    assert len(project.cells()) == 2
    assert project.cells()[1]["region"] == "off_tissue"
    with pytest.raises(ValueError, match="Annotation IDs"):
        import_source("visium", root, tmp_path / "subset", slice_id="s", microns_per_pixel=2,
                      in_tissue_only=True, annotations=labels)


@pytest.mark.parametrize("factor", [0, -1, True, float("inf"), float("nan"), "0.5", 10 ** 1000])
def test_visium_rejects_invalid_calibration(tmp_path, factor):
    with pytest.raises(ValueError):
        readers.import_visium(tmp_path, tmp_path / "project", slice_id="s", microns_per_pixel=factor, in_tissue_only=True)


def test_visium_rejects_hd(tmp_path):
    root = visium_fixture(tmp_path / "visium")
    (root / "binned_outputs").mkdir()
    with pytest.raises(ValueError, match="HD/bin"):
        readers.import_visium(root, tmp_path / "project", slice_id="s", microns_per_pixel=.5, in_tissue_only=True)


def cosmx_fixture(root, *, local=False, mixed=False, unassigned=False):
    coords = "CenterX_local_px,CenterY_local_px" if local else "CenterX_global_px,CenterY_global_px"
    write(root / "metadata.csv", f"fov,cell_ID,{coords},slide_ID,assay_type\n1,1,10,20,s1,RNA\n2,1,100,200,{'s2' if mixed else 's1'},RNA\n")
    write(root / "exprMat.csv", "fov,cell_ID,G1,G2,Negative1,SystemControl1\n2,1,7,0,10,20\n1,1,0,9,30,40\n" + ("1,0,3,4,5,6\n" if unassigned else ""))
    return root


def test_cosmx_global_fov_identity_explicit_gene_scope(tmp_path):
    root = cosmx_fixture(tmp_path / "cosmx")
    labels = write(tmp_path / "annotations.csv", "fov,cell_ID,label\n1,1,T cell\n2,1,Myeloid\n")
    project = readers.import_cosmx(root, tmp_path / "project", slice_id="s", microns_per_pixel=.12,
                                    gene_columns=["G1", "G2"], annotations=labels)
    cells = {cell["cell_id"]: cell for cell in project.cells()}
    assert set(cells) == {"fov=1;cell=1", "fov=2;cell=1"}
    assert cells["fov=1;cell=1"]["counts"] == {"G2": 9}
    assert cells["fov=2;cell=1"]["x"] == 12
    assert cells["fov=2;cell=1"]["y"] == 24
    assert project.summary()["metadata"]["import_parameters"]["excluded_count_columns"] == ["Negative1", "SystemControl1"]


@pytest.mark.parametrize("options,message", [({"local": True}, "FOV-local"), ({"mixed": True}, "Mixed")])
def test_cosmx_rejects_local_or_mixed_slide_coordinates(tmp_path, options, message):
    root = cosmx_fixture(tmp_path / "cosmx", **options)
    with pytest.raises(ValueError, match=message):
        readers.import_cosmx(root, tmp_path / "project", slice_id="s", microns_per_pixel=.12,
                              gene_columns=["G1", "G2"], allow_unannotated=True)


def test_cosmx_unassigned_requires_explicit_exclusion(tmp_path):
    root = cosmx_fixture(tmp_path / "cosmx", unassigned=True)
    with pytest.raises(ValueError, match="unassigned transcripts"):
        readers.import_cosmx(root, tmp_path / "strict", slice_id="s", microns_per_pixel=.12,
                              gene_columns=["G1", "G2"], allow_unannotated=True)
    project = readers.import_cosmx(root, tmp_path / "exclude", slice_id="s", microns_per_pixel=.12,
                                    gene_columns=["G1", "G2"], allow_unannotated=True, exclude_unassigned=True)
    assert len(project.cells()) == 2
    assert project.summary()["metadata"]["import_parameters"]["unassigned_count_rows"] == 1


def test_cosmx_count_context_cannot_contradict_metadata(tmp_path):
    root = cosmx_fixture(tmp_path / "cosmx")
    write(root / "exprMat.csv", "fov,cell_ID,G1,G2,slide_ID\n1,1,2,0,wrong-slide\n2,1,0,4,wrong-slide\n")
    with pytest.raises(ValueError, match="slide_ID does not match metadata"):
        readers.import_cosmx(root, tmp_path / "project", slice_id="s", microns_per_pixel=.12,
                              gene_columns=["G1", "G2"], allow_unannotated=True)


def merscope_fixture(root):
    write(root / "cell_metadata.csv", "EntityID,fov,center_x,center_y\na,1,0.5,-1\nb,2,100,200\n")
    write(root / "cell_by_gene.csv", "cell,G1,G2\nb,7,0\na,0,9\n")
    return root


def test_merscope_global_coordinates_and_count_alignment(tmp_path):
    root = merscope_fixture(tmp_path / "merscope")
    project = import_source("merscope_csv", root, tmp_path / "project", slice_id="s", allow_unannotated=True)
    assert project.cells()[0]["counts"] == {"G2": 9}
    assert project.cells()[0]["x"] == .5 and project.cells()[0]["y"] == -1
    assert project.summary()["metadata"]["coordinate_system"] == "merscope_global_xy"


@pytest.mark.parametrize("content,message", [
    ("cell,G1,G2\nb,7,0\nc,0,9\n", "match exactly"),
    ("cell,G1,G2\na,7,0\na,0,9\n", "Duplicate"),
    ("cell,G1,G2\nb,7.5,0\na,0,9\n", "integers"),
    ("cell,G1,G2\nb,7,0,extra\na,0,9\n", "row width"),
    ("cell,G1,G1\nb,7,0\na,0,9\n", "unique headers"),
])
def test_csv_counts_malformed_or_misaligned(tmp_path, content, message):
    root = merscope_fixture(tmp_path / "merscope")
    write(root / "cell_by_gene.csv", content)
    with pytest.raises(ValueError, match=message):
        readers.import_merscope(root, tmp_path / "project", slice_id="s", allow_unannotated=True)


def test_csv_streaming_count_budget(tmp_path, monkeypatch):
    root = merscope_fixture(tmp_path / "merscope")
    monkeypatch.setattr(readers, "MAX_CSV_COUNT_VALUES", 3)
    with pytest.raises(ValueError, match="streaming alpha budget"):
        readers.import_merscope(root, tmp_path / "project", slice_id="s", allow_unannotated=True)


def test_gzip_decompression_budget_is_checked(tmp_path, monkeypatch):
    root = matrix(tmp_path / "mtx", compressed=True)
    positions, labels = sidecars(tmp_path)
    monkeypatch.setattr(readers, "MAX_TEXT_BYTES", 20)
    with pytest.raises(ValueError, match="Decompressed text"):
        generic(root, tmp_path / "project", positions, labels)


def test_merscope_marker_zero_remains_measured_and_bundle_verifies(tmp_path):
    from spatial_collab.replay import verify_bundle
    root = merscope_fixture(tmp_path / "merscope")
    project = readers.import_merscope(root, tmp_path / "project", slice_id="s", allow_unannotated=True)
    project.set_selection(project.summary()["head_revision"], cell_ids=["a"])
    inspection = project.inspect_selection(["G1", "G2", "NOT_MEASURED"])
    assert inspection["expression"]["G1"]["measured"] is True
    assert inspection["expression"]["G1"]["sum"] == 0
    assert inspection["expression"]["NOT_MEASURED"]["measured"] is False
    bundle = project.export_bundle()
    verification = verify_bundle(bundle["export_path"])
    assert verification["integrity_verified"] is True


def test_h5ad_can_explicitly_start_without_labels(tmp_path):
    ad = pytest.importorskip("anndata")
    import numpy as np
    import pandas as pd
    data = ad.AnnData(np.array([[0, 3], [1, 0]]), obs=pd.DataFrame(index=["a", "b"]), var=pd.DataFrame(index=["G1", "G2"]))
    data.obsm["spatial"] = np.array([[0., 1.], [2., 3.]])
    path = tmp_path / "input.h5ad"
    data.write_h5ad(path)
    project = import_h5ad(path, tmp_path / "project", slice_id="s", coordinate_system="native", units="micrometer",
                          counts_layer="X", allow_unannotated=True, observation_unit="spot", platform="declared")
    assert project.summary()["metadata"]["annotation_status"] == "unannotated"
    assert project.summary()["metadata"]["label_semantics"] == "spot_annotation"


def test_registry_is_lazy_and_declarative():
    environment = os.environ.copy()
    source_path = str(Path(__file__).resolve().parents[1] / "src")
    environment["PYTHONPATH"] = source_path + (os.pathsep + environment["PYTHONPATH"] if environment.get("PYTHONPATH") else "")
    result = subprocess.run([sys.executable, "-c", "import sys,json; from spatial_collab.import_registry import list_formats; print(json.dumps({'formats':len(list_formats()),'heavy':[x for x in ('anndata','spatialdata','h5py','numpy','scipy') if x in sys.modules]}))"],
                            capture_output=True, text=True, env=environment)
    assert result.returncode == 0, f"Registry subprocess failed: {result.stderr}\n{result.stdout}"
    value = json.loads(result.stdout)
    assert value["formats"] == 7 and value["heavy"] == []
    capabilities = list_formats()
    assert {item["id"] for item in capabilities} == {"anndata_h5ad", "xenium", "tenx_mtx", "visium", "cosmx_csv", "merscope_csv", "spatialdata_zarr"}
    capabilities[0]["id"] = "mutated"
    assert list_formats()[0]["id"] != "mutated"


def test_register_reader_explicit_dispatch_and_duplicate_guard(monkeypatch, tmp_path):
    import spatial_collab.import_registry as registry
    monkeypatch.setattr(registry, "_READERS", registry._READERS.copy())
    declaration = {"id": "local_test", "description": "test", "platform": "test", "observation_units": ["cell"],
                   "required_inputs": [], "required_options": [], "dependencies": [], "ceilings": ["test only"]}
    def reader(path, destination, **options):
        return {"path": str(path), "destination": str(destination), "options": options}
    register_reader(declaration, reader)
    assert import_source("local_test", tmp_path, tmp_path / "out", explicit=True)["options"] == {"explicit": True}
    with pytest.raises(ValueError, match="already registered"):
        register_reader(declaration, reader)
    with pytest.raises(ValueError, match="Unknown format"):
        import_source("unknown", tmp_path, tmp_path / "out")


def test_probe_has_no_automatic_choice_and_reports_ambiguity(tmp_path):
    root = matrix(tmp_path / "source")
    write(root / "cells.csv", "cell_id,x_centroid,y_centroid\na,0,0\n")
    probe = probe_source(root)
    assert {item["format_id"] for item in probe["candidates"]} == {"tenx_mtx", "xenium"}
    assert probe["ambiguous"] is True and probe["automatic_selection"] is False
    assert not (root / "project.sqlite3").exists()


def test_probe_cannot_claim_zarr_is_spatialdata(tmp_path):
    root = tmp_path / "generic.zarr"
    root.mkdir()
    probe = probe_source(root)
    assert probe["candidates"][0]["format_id"] == "spatialdata_zarr"
    assert "not yet validated" in probe["candidates"][0]["evidence"][0]
