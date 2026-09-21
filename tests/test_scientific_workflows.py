from copy import deepcopy

import numpy as np
import pytest
from scipy import sparse

from spatial_collab import objects
from spatial_collab.demo import create_demo
from spatial_collab.store import SpatialError
from spatial_collab.workflow_inputs import register_counts, load_counts, snapshot
from spatial_collab.workflow_methods import run
from spatial_collab.workflow_jobs import submit, execute, cancel, inspect
from spatial_collab.spatial_statistics import graph, moran_scores, bh


@pytest.fixture
def inputs(tmp_path):
    p = create_demo(tmp_path / "p")
    # Pure and mixed observations with an exactly recoverable RNA contribution.
    signature = np.array([[80, 5, 10, 5], [5, 80, 5, 10]])
    ref = np.repeat(signature, 8, axis=0)
    fractions = np.array([1, 1, 0.8, 0.8, 0.6, 0.6, 0.4, 0.4, 0.2, 0.2, 0, 0])
    x = np.round(fractions[:, None] * signature[0] + (1 - fractions[:, None]) * signature[1]).astype(int)
    xy = np.array([[i % 4, i // 4] for i in range(len(x))])
    genes = ["A", "B", "C", "D"]
    r = register_counts(
        p,
        x,
        [f"p{i}" for i in range(len(x))],
        genes,
        sample_id="s",
        species="mouse",
        kind="spatial",
        coordinates=xy,
        provenance={"synthetic": True, "units": "array_index"},
    )
    ref_record = register_counts(
        p,
        ref,
        [f"r{i}" for i in range(len(ref))],
        genes,
        sample_id="ref",
        species="mouse",
        kind="reference",
        labels=["typeA"] * 8 + ["typeB"] * 8,
        provenance={"synthetic": True},
    )
    return p, r, ref_record, fractions


def test_sparse_snapshot_preserves_exact_counts_and_rejects_normalized_input(inputs):
    p, r, _, _ = inputs
    record, x = load_counts(p, r["object_id"])
    assert sparse.isspmatrix_csr(x) and x.shape == (12, 4)
    assert record["observation_ids"][0] == "p0"
    with pytest.raises(SpatialError, match="raw counts"):
        register_counts(
            p,
            x / 3,
            record["observation_ids"],
            record["features"],
            sample_id="s",
            species="mouse",
            kind="spatial",
            coordinates=record["coordinates"],
            provenance={"test": True},
        )


def test_nnls_recovers_known_mixtures_and_is_not_named_cell_counts(inputs):
    p, r, ref, truth = inputs
    result = run(
        p,
        {
            "method": "nnls",
            "input_ids": [r["object_id"], ref["object_id"]],
            "parameters": {"features": r["features"]},
        },
    )
    np.testing.assert_allclose(
        np.asarray(result["output"]["rna_contribution_fraction"])[:, 0], truth, atol=0.015
    )
    assert "absolute cell" in result["output"]["interpretation"]
    page = inspect(p, result["object_id"], limit=3)
    assert [r["cell_id"] for r in page["rows"]] == ["p0", "p1", "p2"]
    assert page["next_offset"] == 3 and not page["same_current_project_revision"]


def test_reference_species_and_unknown_features_are_not_silently_intersected(inputs):
    p, r, ref, _ = inputs
    spec = {
        "method": "nnls",
        "input_ids": [r["object_id"], ref["object_id"]],
        "parameters": {"features": ["A", "absent"]},
    }
    with pytest.raises(SpatialError, match="measured"):
        run(p, spec)
    foreign = deepcopy(ref)
    foreign.pop("object_id")
    foreign.pop("object_sha256")
    foreign["species"] = "human"
    foreign = objects.put(p, "analysisinput", foreign)
    with pytest.raises(SpatialError, match="Cross-species"):
        run(p, {**spec, "input_ids": [r["object_id"], foreign["object_id"]]})


def test_moran_matches_direct_quadratic_form_and_bh():
    xy = np.array([[i % 5, i // 5] for i in range(25)])
    w = graph(xy, 4)
    x = np.array([float(i // 5 > 2) for i in range(25)])
    z = x - x.mean()
    assert moran_scores(x[:, None], w)[0] == pytest.approx(25 / w.sum() * (z @ w.toarray() @ z) / (z @ z))
    assert np.isnan(moran_scores(np.ones((25, 1)), w)[0])
    np.testing.assert_allclose(bh([0.01, 0.04, 0.03, 0.2]), [0.04, 0.053333333333, 0.053333333333, 0.2])


def test_sparse_graph_20000_positions_without_dense_quadratic_storage():
    xy = np.array([[i % 200, i // 200] for i in range(20000)])
    w = graph(xy)
    assert sparse.issparse(w) and w.nnz <= 240000
    assert w.data.nbytes + w.indices.nbytes + w.indptr.nbytes < 4_000_000


def test_svg_reproducibility_and_family_wide_fdr(inputs):
    p, r, _, _ = inputs
    spec = {
        "method": "moran_svg",
        "input_ids": [r["object_id"]],
        "parameters": {"neighbors": 2, "permutations": 19, "min_detected": 2},
        "seed": 23,
    }
    a, b = run(p, spec), run(p, spec)
    assert a["output"] == b["output"]
    assert a["output"]["tested_features"] == 4
    assert all(row["q_value"] >= row["p_value"] for row in a["output"]["rows"])


def test_job_failure_cache_and_cancellation_are_durable(inputs):
    p, r, ref, _ = inputs
    spec = {"method": "nnls", "input_ids": [r["object_id"], ref["object_id"]]}
    first = submit(p, spec, launch=False)
    assert execute(p, first["id"])["status"] == "succeeded"
    assert submit(p, spec, launch=False)["cache_hit"]
    bad = submit(p, {**spec, "parameters": {"unknown": True}}, launch=False)
    failure = execute(p, bad["id"])
    assert failure["status"] == "failed" and "Unknown method parameters" in failure["error"]
    job = submit(p, {**spec, "seed": 33}, launch=False)
    assert cancel(p, job["id"])["status"] == "cancelled"
    assert execute(p, job["id"])["status"] == "cancelled"


def test_immutable_workflow_inputs_export_and_offline_recompute(inputs):
    from spatial_collab.replay import verify_bundle

    p, r, ref, _ = inputs
    result = run(p, {"method": "nnls", "input_ids": [r["object_id"], ref["object_id"]]})
    bundle = p.export_bundle(compact=True)
    verified = verify_bundle(bundle["export_path"], recompute_workflows=True)
    assert verified["recompute_matches"] is True
    assert verified["workflow_verification"] == [
        {"result_id": result["object_id"], "method": "nnls", "status": "recomputed_matched"}
    ]


def test_snapshot_does_not_drop_unknown_roi_ids(inputs):
    p, _, _, _ = inputs
    with pytest.raises(SpatialError, match="ROI"):
        snapshot(p, p.context()["head_revision"], "human", ["not-a-real-cell"])


def test_database_mapping_never_guesses_or_merges_symbols():
    from spatial_collab.workflow_methods import _symbol_lookup

    with pytest.raises(SpatialError, match="declared gene symbols"):
        _symbol_lookup({"features": ["ENSG1", "ENSG2"]})
    assert _symbol_lookup({"feature_symbols": ["CD3D", "DUP", "DUP", None]}) == {"CD3D": 0}


def test_count_import_uses_explicit_producer_gene_id_mapping(inputs, tmp_path):
    ad = pytest.importorskip("anndata")
    import pandas as pd
    from spatial_collab.workflow_inputs import import_h5ad

    p, _, _, _ = inputs
    data = ad.AnnData(
        sparse.csr_matrix([[2, 3], [4, 1], [5, 7], [1, 1]]),
        obs=pd.DataFrame({"type": ["A", "A", "B", "B"]}, index=["a", "b", "c", "d"]),
        var=pd.DataFrame({"gene_ids": ["ENSMUSG1", "ENSMUSG2"]}, index=["Gene1", "Gene2"]),
    )
    path = tmp_path / "producer.h5ad"
    data.write_h5ad(path)
    r = import_h5ad(
        p,
        path,
        sample_id="reference",
        species="mouse",
        kind="reference",
        counts_layer="X",
        label_key="type",
        feature_id_key="gene_ids",
        feature_symbol_key="_index",
    )
    assert r["features"] == ["ENSMUSG1", "ENSMUSG2"]
    assert r["feature_symbols"] == ["Gene1", "Gene2"]


def test_reconnect_lists_failed_and_cancelled_jobs(inputs):
    from spatial_collab.workflow_jobs import list_jobs

    p, r, ref, _ = inputs
    spec = {"method": "nnls", "input_ids": [r["object_id"], ref["object_id"]]}
    first = submit(p, spec, launch=False)
    cancel(p, first["id"])
    bad = submit(p, {**spec, "parameters": {"unknown": True}}, launch=False)
    execute(p, bad["id"])
    assert {r["status"] for r in list_jobs(p)["jobs"]} == {"failed", "cancelled"}


def test_harmony_rejects_complete_condition_batch_confounding(inputs):
    pytest.importorskip("harmonypy")
    p, r, _, _ = inputs
    other = {k: v for k, v in r.items() if k not in {"object_id", "object_sha256"}}
    other["sample_id"] = "second-sample"
    other = objects.put(p, "analysisinput", other)
    with pytest.raises(SpatialError, match="confounded"):
        run(
            p,
            {
                "method": "harmony",
                "input_ids": [r["object_id"], other["object_id"]],
                "parameters": {"sample_conditions": ["control", "treated"], "sample_batches": ["A", "B"]},
            },
        )


def test_paste_signature_bridge_preserves_optimizer_and_restores_library():
    ot = pytest.importorskip("ot")
    pytest.importorskip("paste")
    from spatial_collab.paste_compat import line_search_compatibility

    original = ot.optim.cg
    a = np.array([0.5, 0.5])
    cost = np.array([[0.0, 1.0], [1.0, 0.0]])

    def legacy_line_search(fun, coupling, direction, gradient, cost_value, **kwargs):
        return 1.0, 1, fun(coupling + direction)

    with line_search_compatibility() as details:
        if not details["active"]:
            pytest.skip("Signature bridge is intentionally limited to validated upstream versions")
        result = ot.optim.cg(
            a,
            a,
            cost,
            0.01,
            lambda g: (g * g).sum(),
            lambda g: 2 * g,
            line_search=legacy_line_search,
            numItermax=5,
        )
        np.testing.assert_allclose(result, np.diag(a), atol=1e-12)
    assert ot.optim.cg is original
    with pytest.raises(RuntimeError), line_search_compatibility():
        raise RuntimeError("restore even after cancellation/failure")
    assert ot.optim.cg is original


def test_training_tables_have_portable_finite_serialization():
    pd = pytest.importorskip("pandas")
    import json
    from spatial_collab.workflow_methods import _clean

    table = pd.DataFrame({"ELBO": [np.nan, -100.0, -80.0]}, index=[0, 1, 2])
    value = json.loads(json.dumps(_clean(table), allow_nan=False))
    assert value["data"] == [[None], [-100.0], [-80.0]]
