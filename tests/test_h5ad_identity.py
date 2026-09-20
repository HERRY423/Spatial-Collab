"""Raw H5AD imports preserve sample namespaces and the exact feature axis."""
import numpy as np
import pytest
from scipy import sparse

from spatial_collab.identity import describe_feature_query, qualify_observation_id
from spatial_collab.importers import import_h5ad, source_digest


def source(tmp_path, *, ids=("E1", "E2", "E3"), symbols=("Dup", "Dup", "Other"), sample_column=None):
    ad = pytest.importorskip("anndata")
    pd = pytest.importorskip("pandas")
    obs = pd.DataFrame(index=["same-barcode", "b"])
    if sample_column is not None:
        obs["sample_id"] = sample_column
    data = ad.AnnData(X=sparse.csr_matrix([[2, 3, 0], [0, 4, 7]]), obs=obs,
                     var=pd.DataFrame({"gene_ids": list(ids), "symbols": list(symbols)}, index=list(symbols)))
    data.obsm["spatial"] = np.array([[1., 2.], [3., 4.]])
    path = tmp_path / "source.h5ad"
    data.write_h5ad(path)
    return path


def load(path, destination, **kwargs):
    options = dict(slice_id="slice", sample_id="sample-one", coordinate_system="source_array", units="array_index",
                   counts_layer="X", allow_unannotated=True, feature_id_key="gene_ids", observation_unit="spot")
    options.update(kwargs)
    return import_h5ad(path, destination, **options)


def test_two_slices_same_barcode_never_merge_and_duplicate_symbols_preserved(tmp_path):
    path = source(tmp_path)
    digest = source_digest(path)
    first = load(path, tmp_path / "first")
    second = load(path, tmp_path / "second", sample_id="sample-two")
    ids1, ids2 = [{c["cell_id"] for c in p.cells()} for p in (first, second)]
    assert ids1.isdisjoint(ids2)
    cell = next(c for c in first.cells() if c["source_cell_id"] == "same-barcode")
    assert cell["sample_id"] == "sample-one"
    assert cell["cell_id"] == qualify_observation_id("sample-one", "same-barcode")
    assert cell["counts"] == {"E1": 2., "E2": 3.}
    meta = first.summary()["metadata"]
    assert meta["identity_scope"] == "sample_qualified" and meta["observation_unit"] == "spot"
    assert meta["features"] == [{"feature_id": "E1", "symbol": "Dup"}, {"feature_id": "E2", "symbol": "Dup"},
                                {"feature_id": "E3", "symbol": "Other"}]
    assert describe_feature_query(meta, "Dup")["status"] == "ambiguous"
    assert meta["units"] == "array_index" and cell["x"] == 1
    assert source_digest(path) == digest


@pytest.mark.parametrize("units", ["micrometer", "array_index", "pixel", "unknown"])
def test_declared_units_preserve_coordinates_without_scaling(tmp_path, units):
    project = load(source(tmp_path), tmp_path / "project", units=units)
    assert project.summary()["metadata"]["units"] == units
    assert {(c["x"], c["y"]) for c in project.cells()} == {(1., 2.), (3., 4.)}


@pytest.mark.parametrize("options,match", [({"feature_id_key": None}, "Duplicate"),
    ({"feature_id_key": "missing"}, "missing"), ({"feature_symbol_key": "missing"}, "missing"),
    ({"units": "guess"}, "units"), ({"max_features": 2}, "features"), ({"max_nnz": 3}, "nonzero"),
    ({"max_features": True}, "max_features"), ({"max_nnz": 10_000_001}, "max_nnz"),
    ({"max_features": 100_001}, "max_features"), ({"max_nnz": "100"}, "max_nnz")])
def test_import_preflight_rejects_ambiguity_or_outside_declared_budget(tmp_path, options, match):
    with pytest.raises(ValueError, match=match):
        load(source(tmp_path), tmp_path / "project", **options)
    assert not (tmp_path / "project").exists()


def test_duplicate_stable_ids_rejected_even_with_different_symbols(tmp_path):
    with pytest.raises(ValueError, match="Duplicate stable feature ID"):
        load(source(tmp_path, ids=("E1", "E1", "E3")), tmp_path / "project")


def test_sample_declared_does_not_overwrite_source_sample(tmp_path):
    with pytest.raises(ValueError, match="differs from the source"):
        load(source(tmp_path, sample_column=["real", "real"]), tmp_path / "project")


def test_mixed_source_sample_rejected_even_with_explicit_namespace(tmp_path):
    with pytest.raises(ValueError, match="Mixed"):
        load(source(tmp_path, sample_column=["sample-one", "sample-two"]), tmp_path / "project")


def test_custom_slice_column_cannot_hide_mixed_sample_ids(tmp_path):
    ad = pytest.importorskip("anndata")
    path = source(tmp_path, sample_column=["sample-one", "sample-two"])
    data = ad.read_h5ad(path)
    data.obs["declared_slice"] = "same"
    data.write_h5ad(path)
    with pytest.raises(ValueError, match="Mixed"):
        load(path, tmp_path / "project", slice_key="declared_slice")


def test_explicit_feature_budget_preserves_large_complete_axis(tmp_path):
    ad = pytest.importorskip("anndata")
    pd = pytest.importorskip("pandas")
    genes = [f"stable-{i}" for i in range(20001)]
    matrix = sparse.csr_matrix(([3, 9], ([0, 0], [0, 20000])), shape=(1, len(genes)))
    data = ad.AnnData(X=matrix, obs=pd.DataFrame(index=["a"]),
                     var=pd.DataFrame({"gene_ids": genes}, index=[f"G{i}" for i in range(len(genes))]))
    data.obsm["spatial"] = np.array([[0., 0.]])
    path = tmp_path / "full.h5ad"
    data.write_h5ad(path)
    with pytest.raises(ValueError, match="features"):
        load(path, tmp_path / "default")
    project = load(path, tmp_path / "expanded", max_features=30000)
    meta = project.summary()["metadata"]
    assert len(meta["features"]) == len(meta["panel_genes"]) == 20001
    assert project.cells()[0]["counts"] == {"stable-0": 3., "stable-20000": 9.}


def test_feature_symbol_column_explicit_and_budget_recorded(tmp_path):
    project = load(source(tmp_path), tmp_path / "project", feature_symbol_key="symbols", max_nnz=5, max_features=3)
    params = project.summary()["metadata"]["import_parameters"]
    assert params["feature_symbol_key"] == "symbols" and params["max_nnz"] == 5 and params["max_features"] == 3


def test_serialized_source_budget_rejects_before_creating_project(tmp_path, monkeypatch):
    from spatial_collab import importers
    path = source(tmp_path)
    digest = source_digest(path)
    monkeypatch.setattr(importers, "MAX_H5AD_SNAPSHOT_BYTES", 500)
    with pytest.raises(ValueError, match="replayable-source budget"):
        load(path, tmp_path / "project")
    assert not (tmp_path / "project").exists()
    assert source_digest(path) == digest


def test_snapshot_budget_matches_actual_canonical_export_bytes(tmp_path):
    import json
    from spatial_collab.importers import _h5ad_snapshot_budget
    project = load(source(tmp_path), tmp_path / "project")
    cells, metadata = project.cells(), project.summary()["metadata"]
    actual = len(json.dumps({"metadata": metadata, "cells": cells}, ensure_ascii=False,
                            sort_keys=True, separators=(",", ":")).encode("utf-8")) + 1
    assert _h5ad_snapshot_budget(cells, metadata) == actual
