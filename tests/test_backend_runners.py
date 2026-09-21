import json

import numpy as np
import pytest

from spatial_collab.external_runner import export_input_contract, ingest_external_result
from spatial_collab.integration import run_baseline
from spatial_collab.integration_jobs import execute, submit
from spatial_collab.integration_runners import (
    compute_morans_i,
    compute_spatial_weights,
    fit_spatial_graph,
    select_spatial_features,
    smooth_features,
)
from spatial_collab.proteomics import register_protein
from spatial_collab.store import Project


@pytest.fixture
def paired(tmp_path):
    rng = np.random.default_rng(42)
    n = 36
    # 6x6 spatial grid
    cells = [
        {
            "cell_id": f"c{i}",
            "x": float(i % 6),
            "y": float(i // 6),
            "label": "T cell" if i < 18 else "B cell",
            "counts": {
                "G_clustered": 50.0 if (i // 6) < 3 else 2.0,  # Spatially patterned
                "G_noisy": float(rng.integers(1, 40)),  # Spatially random
                "G_low": 1.0,
            },
        }
        for i in range(n)
    ]
    project = Project.create(
        tmp_path / "proj",
        cells,
        {
            "name": "backend runner test",
            "sample_id": "s1",
            "slice_id": "slice_01",
            "source_kind": "synthetic",
            "units": "micrometer",
            "coordinate_system": "cartesian",
            "panel_genes": ["G_clustered", "G_noisy", "G_low"],
        },
    )
    prot_path = tmp_path / "protein.csv"
    with prot_path.open("w", encoding="utf-8") as f:
        f.write("cell_id,x,y,CD3,CD19\n")
        for c in cells:
            f.write(f"{c['cell_id']},{c['x']},{c['y']},{rng.integers(1, 30)},{rng.integers(1, 30)}\n")

    assay = register_protein(
        project,
        prot_path,
        sample_id="s1",
        source_sha256=project.context()["source_sha256"],
        name="TestProtein",
        measurement_type="antibody_count",
        coordinate_system="cartesian",
        units="micrometer",
        registration_note="Synthetic paired protein",
        format_id="csv",
    )
    return project, assay["assay_id"], prot_path


def test_spatial_weights_and_smoothing():
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [5.0, 5.0]])
    W = compute_spatial_weights(coords, k=2)
    assert W.shape == (5, 5)
    # Row normalization check
    np.testing.assert_allclose(W.sum(axis=1), np.ones(5), rtol=1e-5)
    # Diagonal should be zero
    assert np.all(np.diag(W) == 0.0)

    # Feature smoothing test
    X = np.array([[10.0], [10.0], [10.0], [10.0], [0.0]])
    # Alpha = 1 gives identical matrix
    np.testing.assert_allclose(smooth_features(X, W, alpha=1.0), X)
    # Alpha < 1 smooths isolated outlier towards its nearest neighbors
    smoothed = smooth_features(X, W, alpha=0.5)
    assert smoothed[4, 0] > 0.0


def test_morans_i_spatial_autocorrelation():
    coords = np.array([[float(i % 5), float(i // 5)] for i in range(25)])
    W = compute_spatial_weights(coords, k=4)

    # Spatially segregated feature (top 2 rows high, bottom 3 rows low)
    f_clustered = np.array([20.0 if (i // 5) < 2 else 1.0 for i in range(25)]).reshape(-1, 1)
    rng = np.random.default_rng(99)
    f_random = rng.uniform(1.0, 20.0, size=(25, 1))

    I_clustered = compute_morans_i(f_clustered, W)[0]
    I_random = compute_morans_i(f_random, W)[0]

    assert I_clustered > 0.3
    assert I_clustered > I_random


def test_select_spatial_features():
    coords = np.array([[float(i % 5), float(i // 5)] for i in range(25)])

    # Feature 0: high variance but spatially random noise
    rng = np.random.default_rng(123)
    f_noise = rng.normal(0, 10.0, size=(25, 1))
    # Feature 1: moderate variance but strongly spatially clustered
    f_clustered = np.array([10.0 if (i // 5) < 2 else -10.0 for i in range(25)]).reshape(-1, 1)
    cells = [
        {
            "cell_id": f"c_{i}",
            "counts": {
                "G_noise": float(f_noise[i, 0]),
                "G_clustered": float(f_clustered[i, 0]),
            },
        }
        for i in range(25)
    ]

    # Moran's I selection should prioritize G_clustered
    selected = select_spatial_features(cells, ["G_noise", "G_clustered"], coords, max_features=1)
    assert selected == ["G_clustered"]


def test_fit_spatial_graph():
    coords = np.array([[float(i % 6), float(i // 6)] for i in range(36)])
    rng = np.random.default_rng(7)
    rna = rng.uniform(0, 10, size=(36, 5))
    protein = rng.uniform(0, 10, size=(36, 4))

    reps, domains, models = fit_spatial_graph(
        rna,
        protein,
        coords,
        components=3,
        clusters=4,
        seed=42,
    )
    assert len(reps["rna_only"]) == 36
    assert len(reps["protein_only"]) == 36
    assert len(reps["joint"]) == 36
    assert len(reps["joint"][0]) == 3
    assert len(domains["rna_only"]) == 36
    assert len(domains["protein_only"]) == 36
    assert len(domains["joint"]) == 36
    assert set(domains["joint"]).issubset({"0", "1", "2", "3"})
    assert "loadings" in models["joint"]
    assert models["joint"]["spatial_neighbors"] == 6


def test_run_baseline_spatial_graph(paired):
    project, assay_id, _ = paired
    head = project.context()["head_revision"]

    # Run with spatial_graph backend and morans_i feature selection
    result = run_baseline(
        project,
        head,
        assay_id,
        backend="spatial_graph",
        feature_selection="morans_i",
        components=2,
        clusters=3,
        seed=11,
    )
    assert result["method"]["name"] == "spatial_graph_kmeans"
    assert result["method"]["parameters"]["backend"] == "spatial_graph"
    assert result["method"]["parameters"]["feature_selection"] == "morans_i"
    assert len(result["domains"]["joint"]) == 36
    assert len(result["representations"]["joint"]) == 36


def test_integration_job_submission(paired):
    project, assay_id, _ = paired
    head = project.context()["head_revision"]

    spec = {
        "revision_id": head,
        "assay_id": assay_id,
        "backend": "spatial_graph",
        "feature_selection": "morans_i",
        "clusters": 3,
        "components": 2,
        "seed": 2026,
    }
    job = submit(project, spec)
    assert job["status"] in ("queued", "pending", "succeeded")
    res = execute(project, job["id"])
    assert res["status"] == "succeeded"
    assert res["result_id"] is not None


def test_external_runner_export_and_ingest(paired, tmp_path):
    project, assay_id, _ = paired
    head = project.context()["head_revision"]
    export_dir = tmp_path / "export_pkg"

    contract = export_input_contract(project, export_dir, revision_id=head, assay_id=assay_id)
    assert (export_dir / "input_contract.json").exists()
    assert (export_dir / "observations.csv").exists()
    assert (export_dir / "rna_counts.csv").exists()
    assert (export_dir / "protein_counts.csv").exists()
    assert contract["observation_count"] == 36

    # Simulate an external pipeline generating predictions conforming to SCHEMA
    n = contract["observation_count"]
    obs_ids = contract["observation_ids"]
    external_pred = {
        "schema": "spatial-collab.integration.v1",
        "output_origin": "model_prediction",
        "input": contract["input"],
        "preprocessing": {
            "rna": "log1p_cpm",
            "protein": "log1p",
        },
        "method": {
            "name": "external_seurat_v5",
            "version": "5.1.0",
            "code_reference": "https://github.com/satijalab/seurat",
            "environment": {"R": "4.3.0", "Seurat": "5.1.0"},
            "seed": 42,
            "uses_annotations": False,
            "parameters": {
                "backend": "external",
                "notes": "External R Seurat bridge",
            },
        },
        "observation_ids": obs_ids,
        "representations": {
            "rna_only": [[0.1, 0.2] for _ in range(n)],
            "protein_only": [[0.3, 0.4] for _ in range(n)],
            "joint": [[0.5, 0.6] for _ in range(n)],
        },
        "domains": {
            "rna_only": [f"ext_rna_{i % 2}" for i in range(n)],
            "protein_only": [f"ext_prot_{i % 2}" for i in range(n)],
            "joint": [f"ext_joint_{i % 3}" for i in range(n)],
        },
        "fit_scope": "whole_slice",
    }
    pred_path = tmp_path / "external_pred.json"
    pred_path.write_text(json.dumps(external_pred), encoding="utf-8")

    res = ingest_external_result(project, pred_path, execution_evidence={"runtime_seconds": 1.5})
    assert res["status"] == "ingested"
    assert res["result_id"].startswith("integration_")
    assert res["registered"]["method"]["parameters"]["backend"] == "external"
    assert res["registered"]["provenance"]["execution_evidence"]["runtime_seconds"] == 1.5


def test_external_bridge_rejects_unknown_subset_and_implicit_fit_metadata(paired, tmp_path):
    from spatial_collab.store import SpatialError

    project, assay_id, _ = paired
    with pytest.raises(SpatialError, match="Explicit export IDs"):
        export_input_contract(
            project, tmp_path / "bad", assay_id=assay_id, observation_ids=["c0", "not_measured"]
        )
    with pytest.raises(SpatialError, match="environment"):
        ingest_external_result(project, {"method": {"uses_annotations": False}})
    with pytest.raises(SpatialError, match="annotations"):
        ingest_external_result(project, {"method": {"environment": {"python": "test"}}})
