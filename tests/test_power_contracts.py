from copy import deepcopy

import numpy as np
import pytest

from spatial_collab import objects, power
from spatial_collab.method_contracts import contract
from spatial_collab.spatial_statistics import graph, bh
from spatial_collab.store import Project, SpatialError
from spatial_collab.workflow_inputs import register_counts
from spatial_collab.workflow_methods import run, catalog


def assumptions(method="moran_svg", **kwargs):
    return {"model": "gaussian_sar" if method == "moran_svg" else "lognormal_neighbor_coupling",
            "rationale": "Synthetic numerical control, not a biological effect estimate.",
            "simulations": 80, "effect_grid": [0, .5, .95], **kwargs}


@pytest.fixture
def workspace(tmp_path):
    p = Project.create_workspace(tmp_path / "p")
    rng = np.random.default_rng(12)
    x = rng.poisson(5, (36, 3)) + 1
    r = register_counts(p, x, [str(i) for i in range(36)], ["a", "b", "c"],
                        sample_id="control", species="mouse", kind="spatial",
                        coordinates=[[i % 6, i // 6] for i in range(36)], provenance={"synthetic": True})
    return p, {"method": "moran_svg", "input_ids": [r["object_id"]], "parameters": {"permutations": 99}}


def test_bh_resolution_bound_matches_actual_bh_not_bonferroni_claim():
    d = power.resolution(300, 12426, 99)
    assert d["minimum_possible_bh_rejection_rank"] == 2486
    assert d["minimum_permutations_at_assumed_rank"] == 248520
    assert not d["rank_resolution_attainable"]
    for signals, expected in [(1, 0), (2485, 0), (2486, 2486)]:
        q = bh([.01] * signals + [1.] * (12426 - signals))
        assert sum(q < .05) == expected
    assert not power.resolution(300, 100, 19, rank=100)["rank_resolution_attainable"]  # equality rejected
    assert power.resolution(300, 100, 20, rank=100)["rank_resolution_attainable"]


def test_blocked_rank_has_no_finite_mde_without_inventing_effect():
    r = power.simulate({"method": "moran_svg", "observations": 300, "family_size": 12426,
                        "permutations": 99}, None, assumptions())
    assert r["power"] == 0 and r["mde_rho"] is None and r["curve"] == []
    assert r["status"] == "resolution_blocked_at_declared_rank"


def test_moran_power_matches_independent_dense_reference():
    xy = [[i % 3, i // 3] for i in range(9)]
    w = graph(xy, 2)
    d = {"method": "moran_svg", "observations": 9, "family_size": 1, "permutations": 99}
    a = assumptions(simulations=40, seed=19, effect_grid=[0, .8])
    result = power.simulate(d, w, a)
    wd = w.toarray()
    scale = np.diag(1 / np.sqrt(wd.sum(axis=1)))
    def manual(matrix):
        z = matrix - matrix.mean(axis=0)
        return np.array([9 / wd.sum() * sum(wd[i, j]*col[i]*col[j] for i in range(9) for j in range(9)) / sum(col**2) for col in z.T])
    for row in result["curve"]:
        rng = np.random.default_rng(19)
        x = np.linalg.solve(np.eye(9)-row["rho"]*scale@wd@scale, rng.normal(size=(9, 40)))
        observed = manual(x)
        counts = np.zeros(40, dtype=int)
        for _ in range(99):
            counts += manual(rng.permuted(x, axis=0)) >= observed
        assert row["rejections"] == np.sum((counts+1)/100 < .05)
        assert row["generated_effect_median"] == pytest.approx(np.median(observed))


@pytest.mark.parametrize("method", ["moran_svg", "spatial_lr"])
def test_model_power_positive_control_null_and_reproducibility(method):
    xy = [[i % 8, i // 8] for i in range(64)]
    d = {"method": method, "observations": 64, "family_size": 1, "permutations": 99}
    a = assumptions(method, simulations=200, seed=7)
    r = power.simulate(d, graph(xy), a)
    assert r == power.simulate(d, graph(xy), a)
    assert r["curve"][0]["power"] < .12
    assert r["curve"][-1]["power"] > .8
    assert r["mde_rho"] is not None
    assert all(0 <= row["power_ci95"][0] <= row["power"] <= row["power_ci95"][1] <= 1 for row in r["curve"])


def test_freeze_bind_exact_family_and_retrospective_label(workspace):
    p, spec = workspace
    plan = power.create(p, spec, assumptions())
    assert plan["timing"] == "before_target_run_in_this_workspace"
    linked = {**spec, "power_plan_id": plan["object_id"]}
    r = run(p, linked)
    assert r["output"]["power_evidence"]["plan"]["plan_id"] == plan["object_id"]
    assert r["method_contract"]["tier"] == "internal_numerical"
    bad = deepcopy(linked)
    bad["parameters"]["permutations"] = 199
    with pytest.raises(SpatialError, match="does not match"):
        run(p, bad)
    later = power.create(p, spec, assumptions())
    assert later["timing"] == "retrospective_sensitivity"
    from spatial_collab.replay import verify_bundle
    result = verify_bundle(p.export_bundle(compact=True)["export_path"], recompute_workflows=True)
    assert all(row["recompute_matches"] for row in result["power_verification"])
    assert result["workflow_verification"][0]["status"] == "recomputed_matched"


def test_family_is_measured_filtered_and_normalization_constant_aware(workspace):
    p, spec = workspace
    spec["parameters"]["features"] = ["a", "c"]
    d, _ = power.design(p, spec)
    assert d["family_size"] == 2
    r = run(p, spec)
    assert r["output"]["tested_features"] == d["family_size"]
    assert r["output"]["power_evidence"]["plan"] is None


def test_method_partition_and_forged_receipts(workspace):
    c = catalog()
    assert len(c["methods"]) == 11
    groups = {r["method"]: r["method_contract"]["tier"] for r in c["methods"]}
    assert {m for m, tier in groups.items() if tier == "internal_numerical"} == {"nnls", "moran_svg", "spatial_lr"}
    assert all(groups[m] == "external_adapter" for m in ["mofa", "mefisto", "cell2location", "paste", "spagcn", "spatial_spectral", "harmony", "progeny"])
    assert contract("region_comparison")["tier"] == "internal_numerical"
    assert contract("external_bridge:moran_svg")["tier"] == "external_adapter"
    p, spec = workspace
    r = run(p, spec)
    from spatial_collab.workflow_validation import validate_result
    r["environment"]["forged"] = True
    with pytest.raises(SpatialError, match="receipt"):
        validate_result(p, r)


@pytest.mark.parametrize("change", [{"effect_grid": [0, float("nan")]}, {"simulations": True}, {"bh_rank": 0}, {"rationale": ""}, {"model": "observed_power"}])
def test_invalid_assumptions_fail_explicitly(workspace, change):
    p, spec = workspace
    with pytest.raises(SpatialError):
        power.create(p, spec, assumptions(**change))


def test_declaration_precedes_simulation_and_failures_are_retained(workspace, monkeypatch):
    p, spec = workspace
    def fail(*args):
        assert len(objects.catalog(p, "powerdesign")) == 1
        raise SpatialError("numerical failure")
    monkeypatch.setattr(power, "simulate", fail)
    with pytest.raises(SpatialError, match="numerical failure"):
        power.create(p, spec, assumptions())
    assert not objects.catalog(p, "powerplan")


def test_already_submitted_target_cannot_be_called_prospective(workspace):
    from spatial_collab.workflow_jobs import submit
    p, spec = workspace
    job = submit(p, spec, launch=False)
    plan = power.create(p, spec, assumptions())
    assert plan["timing"] == "retrospective_sensitivity"
    assert objects.get(p, plan["design_id"])["prior_job_ids"] == [job["id"]]


def test_custom_lr_matches_independent_edge_sum_and_permutations(monkeypatch):
    from scipy import sparse
    from spatial_collab import workflow_methods as methods
    monkeypatch.setattr(methods, "_lr_pairs", lambda record, params: ([("L", "R")], "synthetic", "0"*64))
    xy = [[i % 3, i // 3] for i in range(9)]
    raw = np.random.default_rng(3).poisson(5, (9, 3)) + 1
    record = {"coordinates": xy, "feature_symbols": ["L", "R", "background"]}
    actual = methods._lr(record, sparse.csr_matrix(raw), {"neighbors": 2, "permutations": 99}, 11, lambda: None)
    y = np.log1p(raw * (10000/raw.sum(axis=1))[:, None])
    w = graph(xy, 2).toarray()
    w /= w.sum(axis=1)[:, None]
    def statistic(b):
        return sum(y[i, 0]*w[i, j]*b[j] for i in range(9) for j in range(9))/9
    observed = statistic(y[:, 1])
    rng = np.random.default_rng(11)
    pvalue = (1+sum(statistic(y[rng.permutation(9), 1]) >= observed for _ in range(99)))/100
    assert actual["rows"][0]["score"] == pytest.approx(observed)
    assert actual["rows"][0]["p_value"] == pvalue
    assert actual["rows"][0]["q_value"] == pvalue
