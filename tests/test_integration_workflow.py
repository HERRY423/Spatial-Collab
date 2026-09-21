from copy import deepcopy
import csv

import numpy as np
import pytest

from spatial_collab import objects
from spatial_collab.integration import run_baseline, register_result, compare_results, fixed_filter
from spatial_collab.integration_jobs import submit, execute, cancel
from spatial_collab.multimodal import register_correspondence, register_layer, register_molecular_relation
from spatial_collab.proteomics import register_protein, inspect_protein, get_assay
from spatial_collab.server import ToolService
from spatial_collab.store import Project, SpatialError


@pytest.fixture
def paired(tmp_path):
    rng = np.random.default_rng(12)
    cells = [{"cell_id": f"c{i}", "x": i % 6, "y": i // 6,
              "label": "Working", "counts": {"R": int(rng.integers(1, 30)), "S": int(rng.integers(1, 30))}} for i in range(30)]
    project = Project.create(tmp_path / "p", cells, {"name": "integration fixture", "sample_id": "s", "slice_id": "s",
        "source_kind": "synthetic", "units": "array_index", "coordinate_system": "array", "panel_genes": ["R", "S"]})
    path = tmp_path / "protein.csv"
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cell_id", "x", "y", "P", "Q"])
        writer.writerows([[c["cell_id"], c["x"], c["y"], int(rng.integers(1, 30)), int(rng.integers(1, 30))] for c in cells])
    assay = register_protein(project, path, sample_id="s", source_sha256=project.context()["source_sha256"], name="Protein",
        measurement_type="antibody_count", coordinate_system="array", units="array_index", registration_note="Synthetic paired observations", format_id="csv")
    return project, assay["assay_id"], path


def fit(paired):
    p, assay, _ = paired
    return run_baseline(p, p.context()["head_revision"], assay, components=2, clusters=3, seed=8)


def test_exact_import_and_permutation_invariance(paired):
    p, _, _ = paired
    result = fit(paired)
    imported = {k: v for k, v in result.items() if not k.startswith("object_")}
    permuted = deepcopy(imported)
    permuted["domains"]["joint"] = ["permuted_" + x for x in result["domains"]["joint"]]
    new = register_result(p, permuted)
    comparison = compare_results(p, result["object_id"], new["object_id"], "joint", "joint")
    assert comparison["adjusted_rand_index"] == pytest.approx(1)
    assert comparison["mean_boundary_disagreement"] == 0
    # External row order may differ if all arrays are reordered together.
    reordered = deepcopy(imported)
    reordered["observation_ids"].reverse()
    for key in ("representations", "domains"):
        for rows in reordered[key].values():
            rows.reverse()
    rr = register_result(p, reordered)
    assert compare_results(p, result["object_id"], rr["object_id"], "joint", "joint")["adjusted_rand_index"] == 1


@pytest.mark.parametrize("fault", ["missing", "duplicate", "foreign", "version", "nan", "feature", "source"])
def test_reject_misleading_external_results(paired, fault):
    p, _, _ = paired
    spec = fit(paired)
    if fault == "missing":
        spec["observation_ids"].pop()
    elif fault == "duplicate":
        spec["observation_ids"][0] = spec["observation_ids"][1]
    elif fault == "foreign":
        spec["observation_ids"][0] = "another_sample:c0"
    elif fault == "version":
        spec["input"]["fit_revision_sha256"] = "0" * 64
    elif fault == "nan":
        spec["representations"]["joint"][0][0] = float("nan")
    elif fault == "feature":
        spec["input"]["rna_features"] = ["unmeasured"]
    else:
        spec["input"]["source_sha256"] = "0" * 64
    with pytest.raises(SpatialError):
        register_result(p, spec)


def test_revision_filter_refit_cache_and_cancellation(paired):
    p, assay, _ = paired
    base = p.context()["head_revision"]
    spec = dict(revision_id=base, assay_id=assay, components=2, clusters=3, seed=8)
    job = submit(p, spec, launch=False)
    done = execute(p, job["id"])
    assert done["status"] == "succeeded"
    assert submit(p, spec, launch=False)["cache_hit"]
    selection = p.set_selection(base, cell_ids=["c0", "c1"])
    proposal = p.propose_revision(base, selection["selection_id"], {"included": False}, "TEST_ONLY synthetic exclusion")
    newrev = p.apply_revision(proposal["proposal_id"], base, "TEST_ONLY", True)["revision_id"]
    filtered = fixed_filter(p, done["result_id"], newrev)
    old = objects.get(p, done["result_id"], "integration")
    assert filtered["input"] == old["input"]
    assert filtered["representations"]["joint"] == old["representations"]["joint"][2:]
    job2 = submit(p, {**spec, "revision_id": newrev}, launch=False)
    assert job2["cache_key"] != job["cache_key"]
    new = execute(p, job2["id"])
    assert new["status"] == "succeeded"
    assert objects.get(p, new["result_id"], "integration")["input"]["fit_revision"] == newrev
    cancelled = submit(p, {**spec, "seed": 9}, launch=False)
    cancel(p, cancelled["id"])
    assert execute(p, cancelled["id"])["status"] == "cancelled"


def test_chunked_protein_and_full_selection_statistics(paired):
    p, _, path = paired
    record = register_protein(p, path, sample_id="s", source_sha256=p.context()["source_sha256"], name="Chunked",
        measurement_type="antibody_count", coordinate_system="array", units="array_index", registration_note="Array test",
        format_id="csv", storage="chunked")
    assay = get_assay(p, record["assay_id"])
    assert len(assay["values"]) == 30
    view = inspect_protein(p, record["assay_id"], p.context()["head_revision"], "P", limit=5, bounds=[0, 0, 2, 2])
    assert len(view["points"]) == 5 and view["total_in_view"] == 9
    assert view["summary"]["measured_count"] == 30
    assert view["next_offset"] == 5
    from spatial_collab.replay import verify_bundle
    bundle = p.export_bundle(compact=True)
    assert verify_bundle(bundle["export_path"], recompute=False)["integrity"] == "verified"
    chunk = p.root / "assays" / record["chunks"][0]["file"]
    with chunk.open("r+b") as stream:
        stream.seek(-8, 2)
        stream.write(b"xxxxxxxx")
    with pytest.raises(SpatialError, match="integrity"):
        get_assay(p, record["assay_id"])


def test_semantics_and_mcp_surfaces(paired):
    p, assay, _ = paired
    relation = register_correspondence(p, {"relation": "same_object", "left_assay": "rna", "right_assay": assay,
        "evidence": "Shared synthetic IDs", "method_version": "fixture1",
        "edges": [{"left_id": "c0", "right_id": "c0", "weight": 1}]})
    assert relation["analysis_permissions"] == ["paired_descriptive"]
    layer = register_layer(p, {"assay_id": assay, "layer_kind": "normalized_corrected", "features": ["P"],
        "observation_ids": ["c0", "c1"], "values": [[-1.2], [None]], "output_scale": "centered", "method": "centering", "method_version": "1"})
    assert layer["values"][0][0] < 0
    molecular = register_molecular_relation(p, {"left_assay": "rna", "left_feature": "R", "right_assay": assay,
        "right_feature": "P", "relation": "marker_association", "evidence": "Test relation", "declared_by": "TEST_ONLY"})
    assert molecular["relation"] == "marker_association"
    assert ToolService(p.root).call("get_multimodal_context")["measurementlayer"]
    with pytest.raises(SpatialError, match="controls"):
        register_layer(p, {**layer, "correction_claim": True, "required_controls": ["isotype"], "controls": {}})


def test_export_recomputes_baselines_and_retains_external_outputs(paired):
    from spatial_collab.replay import verify_bundle
    p, _, _ = paired
    result = fit(paired)
    bundle = p.export_bundle(compact=True)
    verified = verify_bundle(bundle["export_path"])
    assert verified["recompute_matches"]
    assert verified["integration_verification"] == [{"result_id": result["object_id"], "status": "recomputed_matched"}]


def test_large_viewport_pages_do_not_change_statistics(tmp_path):
    n = 10023
    p = Project.create(tmp_path / "large", [{"cell_id": str(i), "x": i, "y": 0, "label": "Working", "counts": {"R": 1}} for i in range(n)],
        {"name": "Large synthetic fixture", "sample_id": "s", "slice_id": "s", "source_kind": "synthetic", "panel_genes": ["R"], "units": "array_index", "coordinate_system": "array"})
    path = tmp_path / "large.csv"
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cell_id", "x", "y", "P"])
        w.writerows([[str(i), i, 0, i % 9] for i in range(n)])
    r = register_protein(p, path, sample_id="s", source_sha256=p.context()["source_sha256"], name="large",
        measurement_type="antibody_count", coordinate_system="array", units="array_index", registration_note="Synthetic scale test", format_id="csv", storage="chunked")
    first = inspect_protein(p, r["assay_id"], p.context()["head_revision"], "P")
    second = inspect_protein(p, r["assay_id"], p.context()["head_revision"], "P", offset=first["next_offset"])
    assert len(first["points"]) == 10000 and len(second["points"]) == 23
    assert first["summary"] == second["summary"]
    assert len({p["cell_id"] for p in first["points"] + second["points"]}) == n


def test_study_counts_subjects_not_sections_and_blocks_confounding(paired):
    from spatial_collab.study import register_study, compare_study
    p, assay, _ = paired
    rows = [{"subject_id": subject, "sample_id": "s" + str(i), "section_id": "section" + str(i),
             "condition": condition, "batch": condition, "sample_source": "synthetic", "roi_definition": "fixed rectangle",
             "roi_class": "tissue_region", "roi_origin": "prespecified", "source_sha256": "a" * 64,
             "assay_ids": [assay], "run_ids": ["run" + str(i)], "pair_id": None}
            for i, (subject, condition) in enumerate([("patient1", "case"), ("patient1", "case"), ("patient2", "control")])]
    study = register_study(p, {"study_id": "synthetic-study", "sections": rows})
    assert study["subjects_per_condition"] == {"case": 1, "control": 1}
    results = [{"section_id": r["section_id"], "source_sha256": r["source_sha256"], "run_id": r["run_ids"][0],
                "value": i, "metric": "mean_signal", "units": "counts", "method_version": "1"} for i, r in enumerate(rows)]
    compared = compare_study(p, study["object_id"], results)
    assert compared["inference_status"] == "blocked_condition_batch_confounding"
    assert len(compared["subject_summaries"]) == 2
    assert compared["subject_summaries"][0]["mean_of_section_summaries"] == .5


def test_worker_subprocess_durable_result(paired):
    import time
    from spatial_collab.integration_jobs import get_job
    p, assay, _ = paired
    job = submit(p, {"revision_id": p.context()["head_revision"], "assay_id": assay, "components": 2, "clusters": 3})
    deadline = time.monotonic() + 30
    while time.monotonic() < deadline:
        current = get_job(p, job["id"])
        if current["status"] in {"succeeded", "failed"}:
            break
        time.sleep(.1)
    assert current["status"] == "succeeded", current
    assert current["worker_pid"] != __import__("os").getpid()


def test_frozen_hypotheses_extend_to_all_integration_alternatives(paired, monkeypatch):
    from spatial_collab.hypotheses import create_hypothesis_plan, freeze_hypothesis_plan
    from spatial_collab import integration_sensitivity
    p, assay, _ = paired
    revision = p.context()["head_revision"]
    selection = p.set_selection(revision, cell_ids=["c0", "c1"])
    plan = create_hypothesis_plan(p, {"name": "Integration hypotheses", "rationale": "Synthetic alternative test",
        "revision_id": revision, "selection_id": selection["selection_id"], "source_label": "Working", "target_label": "Other",
        "radii_um": [10], "graph_scopes": ["whole_slice"], "min_effect": .1, "background": {"kind": "imported_universe"},
        "variants": [{"name": "exclude", "rationale": "Test population sensitivity", "selector": {"cell_ids": ["c0"]}, "changes": {"included": False}}]})
    frozen = freeze_hypothesis_plan(p, plan["plan_id"], plan["version"])
    monkeypatch.setattr(integration_sensitivity, "submit", lambda project, spec: submit(project, spec, launch=False))
    run = integration_sensitivity.run_plan(p, frozen["plan_id"], frozen["version"], assay, [
        {"backend": "balanced_pca", "protein_transform": scale, "components": 2, "clusters": 3, "seed": 0} for scale in ("log1p", "asinh")])
    assert len(run["rows"]) == 4 and run["all_alternatives_retained"]
    assert p.context()["head_revision"] == revision
    for row in run["rows"]:
        assert execute(p, row["job_id"])["status"] == "succeeded"
    assert {r["status"] for r in integration_sensitivity.inspect_plan_run(p, run["object_id"])["rows"]} == {"succeeded"}


def test_explicit_coordinate_conversion_and_no_implicit_tolerance(paired):
    p, _, path = paired
    rows = list(csv.DictReader(path.open()))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["cell_id", "x", "y", "P", "Q"])
        writer.writeheader()
        writer.writerows([{**r, "x": float(r["x"]) + 1} for r in rows])
    options = dict(sample_id="s", source_sha256=p.context()["source_sha256"], name="converted",
        measurement_type="antibody_count", coordinate_system="array", units="array_index", registration_note="Known format translation",
        format_id="csv")
    with pytest.raises(SpatialError, match="coordinates"):
        register_protein(p, path, **options)
    result = register_protein(p, path, **options, coordinate_transform=[[1, 0, -1], [0, 1, 0], [0, 0, 1]])
    assert result["coordinate_conversion"]["maximum_absolute_residual"] == 0


def test_import_cannot_bypass_fixed_model_uncertainty_or_annotation_contract(paired):
    p, _, _ = paired
    original = fit(paired)
    original["uncertainty"] = {"score": [0.1] * len(original["observation_ids"])}
    original["uncertainty_meaning"] = {"score": "Synthetic non-probabilistic score"}
    original = register_result(p, original)
    spec = fixed_filter(p, original["object_id"], p.context()["head_revision"])
    spec["uncertainty"]["score"][0] = .8
    with pytest.raises(SpatialError, match="Fixed-model"):
        register_result(p, spec)
    base = p.context()["head_revision"]
    selected = p.set_selection(base, cell_ids=["c0"])
    proposal = p.propose_revision(base, selected["selection_id"], {"label": "Other"}, "TEST_ONLY")
    revision = p.apply_revision(proposal["proposal_id"], base, "TEST_ONLY", True)["revision_id"]
    original["method"]["uses_annotations"] = True
    supervised = register_result(p, original)
    forged = {**supervised, "mode": "fixed_model_filter", "view_revision": revision, "parent_result_id": supervised["object_id"]}
    with pytest.raises(SpatialError, match="Annotation-dependent"):
        register_result(p, forged)
