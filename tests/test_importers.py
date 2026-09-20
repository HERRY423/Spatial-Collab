import numpy as np
import pytest
from scipy import sparse

from spatial_collab.importers import _matrix_counts, import_h5ad, import_xenium, source_digest


def h5ad_fixture(tmp_path, *, fractional=False, mixed=False, duplicate=False):
    ad = pytest.importorskip("anndata")
    pd = pytest.importorskip("pandas")
    matrix = sparse.csr_matrix([[2.5 if fractional else 2, 0], [0, 4]])
    data = ad.AnnData(X=matrix.copy(), obs=pd.DataFrame({"cell_type": ["T cell", "Myeloid"],
        "slice_id": ["s1", "s2" if mixed else "s1"]}, index=["a", "b"]),
        var=pd.DataFrame(index=["CD3D", "CD3D" if duplicate else "LYZ"]))
    data.layers["counts"] = matrix
    data.obsm["spatial"] = np.array([[0., 0.], [1., 0.]])
    path = tmp_path / "cells.h5ad"
    data.write_h5ad(path)
    return path


def do_h5ad(path, target, **kwargs):
    return import_h5ad(path, target, slice_id="s1", coordinate_system="native", units="micrometer", **kwargs)


def test_h5ad_raw_counts_and_source_unchanged(tmp_path):
    path = h5ad_fixture(tmp_path)
    digest = source_digest(path)
    project = do_h5ad(path, tmp_path / "project")
    assert project.cells()[0]["counts"].get("LYZ", 0) == 0
    assert project.cells()[1]["counts"]["LYZ"] == 4
    assert project.summary()["metadata"]["biological_replicates"] == 0
    assert source_digest(path) == digest


@pytest.mark.parametrize("options,match", [({"fractional": True}, "integer"), ({"mixed": True}, "slice"),
                                          ({"duplicate": True}, "Duplicate")])
def test_h5ad_rejects_ambiguous_inputs(tmp_path, options, match):
    path = h5ad_fixture(tmp_path, **options)
    with pytest.raises(ValueError, match=match):
        do_h5ad(path, tmp_path / "project")
    assert not (tmp_path / "project").exists()


def test_h5ad_missing_layer_no_silent_fallback(tmp_path):
    path = h5ad_fixture(tmp_path)
    with pytest.raises(ValueError, match="missing"):
        do_h5ad(path, tmp_path / "project", counts_layer="absent")


def test_h5ad_units_not_inferred(tmp_path):
    path = h5ad_fixture(tmp_path)
    project = import_h5ad(path, tmp_path / "project", slice_id="s", coordinate_system="x", units="pixel")
    assert project.summary()["metadata"]["units"] == "pixel"
    assert project.cells()[1]["x"] == 1


def xenium_fixture(tmp_path):
    h5py = pytest.importorskip("h5py")
    outs = tmp_path / "outs"
    outs.mkdir()
    (outs / "cells.csv").write_text("cell_id,x_centroid,y_centroid\nb,2,0\na,0,0\n", encoding="utf-8")
    annotations = tmp_path / "labels.csv"
    annotations.write_text("cell_id,label\na,T cell\nb,Myeloid\n", encoding="utf-8")
    # Matrix barcode order intentionally differs from cells.csv.
    mat = sparse.csc_matrix([[7, 0], [0, 9], [99, 88]])
    with h5py.File(outs / "cell_feature_matrix.h5", "w") as handle:
        group = handle.create_group("matrix")
        for key in ("data", "indices", "indptr"):
            group[key] = getattr(mat, key)
        group["shape"] = mat.shape
        group["barcodes"] = np.array([b"a", b"b"])
        features = group.create_group("features")
        features["name"] = np.array([b"CD3D", b"LYZ", b"control"])
        features["feature_type"] = np.array([b"Gene Expression", b"Gene Expression", b"Negative Control Probe"])
    return outs, annotations


def test_xenium_barcode_alignment_controls_excluded(tmp_path):
    outs, annotations = xenium_fixture(tmp_path)
    project = import_xenium(outs, annotations, tmp_path / "project", slice_id="s1")
    cells = {c["cell_id"]: c for c in project.cells()}
    assert cells["a"]["counts"] == {"CD3D": 7}
    assert cells["b"]["counts"] == {"LYZ": 9}
    assert set(project.summary()["metadata"]["panel_genes"]) == {"CD3D", "LYZ"}


def test_xenium_missing_annotation_rejected(tmp_path):
    outs, annotations = xenium_fixture(tmp_path)
    annotations.write_text("cell_id,label\na,T cell\n", encoding="utf-8")
    with pytest.raises(ValueError, match="match Xenium"):
        import_xenium(outs, annotations, tmp_path / "project", slice_id="s1")


def test_xenium_duplicate_annotation_rejected(tmp_path):
    outs, annotations = xenium_fixture(tmp_path)
    annotations.write_text("cell_id,label\na,T cell\na,Myeloid\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate"):
        import_xenium(outs, annotations, tmp_path / "project", slice_id="s1")


@pytest.mark.parametrize("case", ["short_feature_types", "short_barcodes"])
def test_xenium_metadata_dimensions_no_silent_dropping(tmp_path, case):
    h5py = pytest.importorskip("h5py")
    outs, annotations = xenium_fixture(tmp_path)
    with h5py.File(outs / "cell_feature_matrix.h5", "a") as handle:
        group = handle["matrix"]
        if case == "short_feature_types":
            del group["features"]["feature_type"]
            group["features"]["feature_type"] = np.array([b"Gene Expression", b"Gene Expression"])
        else:
            group["shape"][:] = [3, 3]
    with pytest.raises(ValueError, match="dimensions"):
        import_xenium(outs, annotations, tmp_path / "project", slice_id="s1")
    assert not (tmp_path / "project").exists()


def test_unsigned_duplicate_counts_do_not_wrap():
    matrix = sparse.csc_matrix((np.array([3_000_000_000, 3_000_000_000], dtype=np.uint32),
                                np.array([0, 0]), np.array([0, 2])), shape=(1, 1))
    assert _matrix_counts(matrix, ["a"], ["gene"]) == [{"gene": 6_000_000_000.0}]


def test_counts_exceeding_exact_integer_range_rejected():
    with pytest.raises(ValueError, match="exact integer"):
        _matrix_counts(sparse.csr_matrix([[2**54]]), ["a"], ["gene"])
