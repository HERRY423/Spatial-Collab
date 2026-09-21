import numpy as np
import pytest
from scipy import stats

from spatial_collab.demo import create_demo
from spatial_collab.study import register_study
from spatial_collab.study_inference import infer
from spatial_collab.store import SpatialError


def design_data(p, paired=False, batches=False, n=12, sections=3):
    rng = np.random.default_rng(93)
    design, records = [], []
    for i in range(n):
        u = rng.normal(0, 2)
        conditions = ["control", "case"] if paired else ["control" if i < n//2 else "case"]
        for condition in conditions:
            for j in range(sections):
                sid = f"s{i}-{condition}-{j}"
                batch = "b" + str(j % 2) if batches else "b0"
                design.append({"subject_id": str(i), "sample_id": f"sample{i}-{condition}", "section_id": sid, "condition": condition, "batch": batch, "sample_source": "synthetic_known_model", "roi_definition": "fixed tissue ROI", "roi_class": "tissue", "roi_origin": "prespecified", "source_sha256": "a"*64, "assay_ids": ["a"], "run_ids": [sid]})
                records.append({"section_id": sid, "source_sha256": "a"*64, "run_id": sid, "metric": "signal", "units": "score", "method_version": "control-v1", "value": float(u+(condition == "case")*1.5+(batch == "b1")*.5+rng.normal(0, .3))})
    d = register_study(p, {"study_id": "known_control", "sections": design})
    return d, records


def plan(method):
    return {"method": method, "control": "control", "case": "case", "roi_class": "tissue", "metrics": ["signal"], "outcome_scale": "continuous_section_summary", "rationale": "Synthetic numerical correctness control"}


def test_welch_uses_subject_not_section_replication(tmp_path):
    p = create_demo(tmp_path / "project")
    d, rows = design_data(p)
    r = infer(p, d["object_id"], rows, plan("welch"))["results"][0]
    a = [np.mean([x["value"] for x in rows if x["section_id"].startswith(f"s{i}-")]) for i in range(6)]
    b = [np.mean([x["value"] for x in rows if x["section_id"].startswith(f"s{i}-")]) for i in range(6, 12)]
    v = np.var(a, ddof=1)/6+np.var(b, ddof=1)/6
    df = v*v/((np.var(a, ddof=1)/6)**2/5+(np.var(b, ddof=1)/6)**2/5)
    expected = 2*stats.t.sf(abs((np.mean(b)-np.mean(a))/np.sqrt(v)), df)
    assert r["p_value"] == pytest.approx(expected)
    assert r["subjects_per_condition"] == {"control": 6, "case": 6}


def test_paired_test_preserves_subject_pairs(tmp_path):
    p = create_demo(tmp_path / "project")
    d, rows = design_data(p, paired=True)
    r = infer(p, d["object_id"], rows, plan("paired_t"))["results"][0]
    assert r["status"] == "tested"
    assert r["degrees_of_freedom"] == 11
    assert r["effect"] == pytest.approx(1.5, abs=.2)
    bad = infer(p, d["object_id"], rows, plan("welch"))["results"][0]
    assert bad["p_value"] is None and "both conditions" in bad["reason"]


def test_mixedlm_estimates_known_effect_and_batch(tmp_path):
    p = create_demo(tmp_path / "project")
    d, rows = design_data(p, paired=True, batches=True, n=24)
    r = infer(p, d["object_id"], rows, plan("mixedlm"))["results"][0]
    assert r["status"] == "tested", r
    assert r["converged"] and r["effect"] == pytest.approx(1.5, abs=.15)
    assert r["fixed_effect_estimates"][2] == pytest.approx(.5, abs=.15)
    assert r["subject_count"] == 24 and r["section_count"] == 144


def test_missing_family_members_are_retained_and_counted(tmp_path):
    p = create_demo(tmp_path / "project")
    d, rows = design_data(p)
    r = infer(p, d["object_id"], rows, {**plan("welch"), "metrics": ["signal", "missing"]})
    assert r["results"][1]["status"] == "not_tested"
    assert r["family_size"] == 2
    assert r["results"][0]["q_value"] == min(1, r["results"][0]["p_value"]*2)
    excluded = infer(p, d["object_id"], rows[1:], {**plan("welch"), "missing_policy": "complete_subjects"})
    assert excluded["results"][0]["excluded_subjects"] == ["0"]
    assert excluded["results"][0]["subjects_per_condition"]["control"] == 5


def test_confounding_and_source_mismatch_cannot_produce_pvalues(tmp_path):
    p = create_demo(tmp_path / "project")
    d, rows = design_data(p)
    for s in d["sections"]:
        s["batch"] = s["condition"]
    dd = register_study(p, {"study_id": "confounded", "sections": d["sections"]})
    r = infer(p, dd["object_id"], rows, plan("mixedlm"))["results"][0]
    assert r["p_value"] is None and "rank deficient" in r["reason"]
    rows[0]["source_sha256"] = "b"*64
    with pytest.raises(SpatialError, match="identity"):
        infer(p, d["object_id"], rows, plan("welch"))


def test_study_result_offline_replay(tmp_path):
    from spatial_collab.replay import verify_bundle
    p = create_demo(tmp_path / "project")
    d, rows = design_data(p)
    infer(p, d["object_id"], rows, plan("welch"))
    bundle = p.export_bundle(compact=True)
    result = verify_bundle(bundle["export_path"])
    assert result["study_verification"][0]["recompute_matches"]
