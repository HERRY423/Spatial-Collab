"""Design-conditional permutation resolution and simulation power, never observed power.

BH rank is a declared sufficient rejection threshold, NOT a predicted number of
discoveries. A grid MDE is model/graph conditional, with Monte Carlo uncertainty.
"""
from fractions import Fraction
import hashlib
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from scipy.stats import norm
from threadpoolctl import threadpool_limits

from . import objects
from .spatial_statistics import graph, moran_scores, normalized_counts
from .store import SpatialError, _hash, _now


def resolution(n, family_size, permutations, alpha=0.05, rank=1):
    for name, value, low in (("observations", n, 3), ("family_size", family_size, 1),
                             ("permutations", permutations, 1), ("rank", rank, 1)):
        if type(value) is not int or value < low:
            raise SpatialError(f"{name} must be an integer >= {low}.")
    if rank > family_size or type(alpha) not in (int, float) or not 0 < alpha < 1:
        raise SpatialError("Require rank <= family size and 0 < alpha < 1.")
    a = Fraction(str(alpha))
    floor = Fraction(1, permutations + 1)
    threshold = a * rank / family_size
    # Strict q < alpha throughout. Fractions avoid exact-boundary rounding.
    required_rank = int(Fraction(family_size, permutations + 1) / a) + 1
    required_permutations = int(Fraction(family_size, rank) / a)
    return {"observations": n, "family_size": family_size, "permutations": permutations,
            "alpha": alpha, "decision_rule": "q < alpha", "minimum_p": float(floor),
            "assumed_bh_rank": rank, "rank_threshold": float(threshold),
            "minimum_possible_bh_rejection_rank": required_rank,
            "minimum_permutations_at_assumed_rank": required_permutations,
            "current_analysis_permutation_cap": 9999,
            "required_budget_within_current_runner": required_permutations <= 9999,
            "rank_resolution_attainable": floor < threshold,
            "rank_scope": "Sufficient BH cutoff conditional on declared rejection rank; actual BH rank is unknown. Other signals can rescue a lower-ranked p-value.",
            "effect_size_from_counts_alone": None,
            "interpretation": "Sample size and p-value resolution do not uniquely determine a minimum detectable effect. Zero discoveries do not establish biological absence."}


def _int(d, key, default, low, high):
    v = d.get(key, default)
    if type(v) is not int or not low <= v <= high:
        raise SpatialError(f"{key} must be an integer in {low}..{high}.")
    return v


def design(project, spec):
    """Resolve the exact test family without evaluating observed test statistics."""
    from .workflow_inputs import load_counts
    from .workflow_methods import METHOD_PARAMETERS, _lr_pairs, _symbol_lookup
    if not isinstance(spec, dict) or set(spec) - {"method", "input_ids", "parameters", "seed", "power_plan_id"}:
        raise SpatialError("Power design requires a workflow recipe.")
    method = spec.get("method")
    if method not in {"moran_svg", "spatial_lr"} or len(spec.get("input_ids", [])) != 1:
        raise SpatialError("Power planning needs one spatial input and Moran or spatial_lr.")
    params = spec.get("parameters", {})
    if not isinstance(params, dict) or set(params) - set(METHOD_PARAMETERS[method]):
        raise SpatialError("Unknown power design method parameters.")
    record, x = load_counts(project, spec["input_ids"][0])
    if record["kind"] != "spatial":
        raise SpatialError("Power planning requires measured spatial coordinates.")
    n = x.shape[0]
    k = _int(params, "neighbors", 6, 1, min(64, n - 1))
    b = _int(params, "permutations", 99, 19, 9999)
    w = graph(record["coordinates"], k)
    resource = None
    if method == "moran_svg":
        features = params.get("features", record["features"])
        if not isinstance(features, list) or not features or not all(isinstance(g, str) for g in features) or len(set(features)) != len(features) or not set(features) <= set(record["features"]):
            raise SpatialError("Power family must contain unique measured feature IDs.")
        minimum = _int(params, "min_detected", 5, 2, n)
        detected = np.asarray((x > 0).sum(axis=0)).ravel()
        lookup = {g: i for i, g in enumerate(record["features"])}
        ix = [lookup[g] for g in features if detected[lookup[g]] >= minimum]
        y = normalized_counts(x)
        family = []
        for start in range(0, len(ix), 64):
            cols = ix[start:start + 64]
            z = y[:, cols].toarray()
            valid = ((z - z.mean(axis=0)) ** 2).sum(axis=0) > 1e-12
            family.extend(record["features"][i] for i, keep in zip(cols, valid) if keep)
    else:
        pairs, _, resource = _lr_pairs(record, params)
        lookup = _symbol_lookup(record)
        family = [list(p) for p in pairs if set(p[0].split("_")) | set(p[1].split("_")) <= set(lookup)]
    if not family:
        raise SpatialError("No eligible tests in the declared family; no MDE can be estimated.")
    resolved = {"method": method, "input_id": record["object_id"], "input_sha256": record["object_sha256"],
                "observations": n, "neighbors": k, "permutations": b, "family_size": len(family),
                "family_sha256": _hash(sorted(family)), "resource_sha256": resource,
                "coordinates_sha256": _hash(record["coordinates"]),
                "eligibility": "Exact measured/filter/constant eligibility; no observed Moran/LR statistic or p-value is used."}
    return resolved, w


def _wilson(successes, n):
    p = successes / n
    z = norm.ppf(.975)
    center = (p + z*z/(2*n)) / (1 + z*z/n)
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/(1+z*z/n)
    return [float(max(0, center-half)), float(min(1, center+half))]


@threadpool_limits.wrap(limits=1, user_api="blas")
def simulate(resolved, weights, assumptions, checkpoint=lambda: None):
    if not isinstance(assumptions, dict) or set(assumptions) - {"alpha", "bh_rank", "target_power", "effect_grid", "simulations", "seed", "model", "rationale"}:
        raise SpatialError("Unknown power assumptions.")
    expected_model = "gaussian_sar" if resolved["method"] == "moran_svg" else "lognormal_neighbor_coupling"
    if assumptions.get("model") != expected_model:
        raise SpatialError(f"Declare model={expected_model}; there is no universal MDE.")
    if not isinstance(assumptions.get("rationale"), str) or not assumptions["rationale"].strip():
        raise SpatialError("Declare why this effect model and BH rank are useful for the proposed design.")
    rank = _int(assumptions, "bh_rank", 1, 1, resolved["family_size"])
    nsim = _int(assumptions, "simulations", 200, 40, 2000)
    seed = _int(assumptions, "seed", 0, 0, 2**32-1)
    target = assumptions.get("target_power", .8)
    if type(target) not in (float, int) or not .5 <= target < 1:
        raise SpatialError("Choose target_power in [0.5,1).")
    grid = assumptions.get("effect_grid", [0, .2, .4, .6, .8, .95])
    if not isinstance(grid, list) or not 2 <= len(grid) <= 20 or any(type(v) not in (int, float) or not 0 <= v < 1 for v in grid) or grid != sorted(set(grid)) or grid[0] != 0:
        raise SpatialError("effect_grid must be 2..20 increasing unique rho values in [0,1), including zero.")
    n, b = resolved["observations"], resolved["permutations"]
    limits = resolution(n, resolved["family_size"], b, assumptions.get("alpha", .05), rank)
    common = {"resolution": limits, "model": expected_model, "target_power": target,
              "effect_definition": "rho in the declared generative model; not log fold change, a universal Moran I, or biochemical interaction strength",
              "scope": "Power for p < alpha * declared_rank / family_size, a sufficient BH rejection condition; not full-family BH power or study-level power.",
              "uncertainty": "Pointwise 95% Wilson intervals, not simultaneous confidence bounds or a validated sample-size guarantee.",
              "biological_calibration": "NOT_ESTABLISHED; design simulation uses graph and declared model, not observed effect estimates, count dropout, complex subunit noise or tissue-specific nuisances."}
    if not limits["rank_resolution_attainable"]:
        return {**common, "status": "resolution_blocked_at_declared_rank", "power": 0,
                "mde_rho": None, "mde_reason": "No finite effect can meet this declared rank cutoff with the chosen permutation p-value grid.", "curve": []}
    if n > 5000 or n * nsim * b * len(grid) > 300_000_000:
        raise SpatialError("Power simulation budget exceeded; explicitly reduce simulations/effect grid/ROI. No subsampling executed.")
    roww = sparse.diags(1 / np.asarray(weights.sum(axis=1)).ravel()) @ weights
    sym = sparse.diags(1 / np.sqrt(np.asarray(weights.sum(axis=1)).ravel()))
    sym = sym @ weights @ sym
    curve = []
    for rho in grid:
        checkpoint()
        rng = np.random.default_rng(seed)  # common random numbers across effects, independent replicate columns
        e = rng.normal(size=(n, nsim))
        if resolved["method"] == "moran_svg":
            z = splu((sparse.eye(n) - rho * sym).tocsc()).solve(e)
            observed = moran_scores(z, weights)
            def score(shuffled):
                return moran_scores(shuffled, weights)
            effect = observed
        else:
            ligand = np.exp(e - .5)
            neighbor = roww.T @ ligand
            neighbor = (neighbor - neighbor.mean(axis=0)) / neighbor.std(axis=0)
            z = np.exp(rho * neighbor + np.sqrt(1-rho*rho) * rng.normal(size=(n, nsim)) - .5)
            def score(shuffled):
                return np.sum(ligand * (roww @ shuffled), axis=0) / n
            observed = score(z)
            effect = observed / (ligand.mean(axis=0) * z.mean(axis=0)) - 1
        exceed = np.zeros(nsim, dtype=int)
        for perm in range(b):
            if perm % 100 == 0:
                checkpoint()
            exceed += score(rng.permuted(z, axis=0)) >= observed
        p = (exceed + 1) / (b + 1)
        successes = int(np.sum(p < limits["rank_threshold"]))
        curve.append({"rho": rho, "rejections": successes, "simulations": nsim,
                      "power": successes/nsim, "power_ci95": _wilson(successes, nsim),
                      "generated_effect_median": float(np.median(effect)),
                      "generated_effect_range90": np.quantile(effect, [.05, .95]).tolist()})
    candidates = [i for i, row in enumerate(curve) if row["rho"] > 0 and all(x["power_ci95"][0] >= target for x in curve[i:])]
    mde = curve[candidates[0]] if candidates else None
    return {**common, "status": "simulated", "curve": curve, "mde_rho": mde["rho"] if mde else None,
            "mde_generated_effect_median": mde["generated_effect_median"] if mde else None,
            "generated_effect_name": "Moran I" if resolved["method"] == "moran_svg" else "LR score / shuffled-placement expectation - 1",
            "mde_reason": "Smallest tested nonzero rho whose pointwise lower confidence bound and all larger grid points meet target; no interpolation or monotonicity assumption. Not an exact continuous minimum." if mde else "Target not demonstrated within the declared model/grid; not proof of no biological signal."}


def create(project, spec, assumptions):
    resolved, w = design(project, spec)
    prior = [oid for oid in objects.catalog(project, "analysisresult") if
             (r := objects.get(project, oid))["method"] == resolved["method"] and
             resolved["input_id"] in r["recipe"]["input_ids"]]
    import json
    with project._db() as db:
        has_jobs = db.execute("SELECT 1 FROM sqlite_master WHERE name='analysis_jobs'").fetchone()
        prior_jobs = [row[0] for row in db.execute("SELECT id,recipe FROM analysis_jobs") if
                      (recipe := json.loads(row[1])).get("method") == resolved["method"] and
                      resolved["input_id"] in recipe.get("input_ids", [])] if has_jobs else []
    from .integration import numerical_environment
    env = numerical_environment()
    # Immutable declaration is saved BEFORE Monte Carlo or the target analysis.
    declaration = objects.put(project, "powerdesign", {"resolved": resolved, "assumptions": assumptions,
                               "created_at": _now(), "prior_result_ids": prior, "prior_job_ids": prior_jobs,
                               "timing": "retrospective_sensitivity" if prior or prior_jobs else "before_target_run_in_this_workspace",
                               "timing_scope": "Local record ordering only; no claim that the researcher has not seen these data elsewhere.",
                               "environment": env, "environment_sha256": _hash(env),
                               "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    result = simulate(resolved, w, assumptions)
    return objects.put(project, "powerplan", {"design_id": declaration["object_id"], "resolved": resolved,
                       "assumptions": assumptions, "timing": declaration["timing"], "result": result})


def bind(project, spec):
    if not spec.get("power_plan_id"):
        return None
    plan = objects.get(project, spec["power_plan_id"], "powerplan")
    resolved, _ = design(project, spec)
    if resolved != plan["resolved"]:
        raise SpatialError("Power plan does not match the exact input, graph, test family or permutation budget. Freeze a new plan.")
    return {"plan_id": plan["object_id"], "plan_sha256": plan["object_sha256"],
            "timing": plan["timing"], "result": plan["result"]}


def create_configured(project, spec, assumptions):
    """Use the same configured scientific runtime as the eventual analysis."""
    import json
    import os
    import subprocess
    from .workflow_jobs import runtime
    from .store import _json
    env = {**os.environ, "PYTHONPATH": str(Path(__file__).resolve().parents[1]),
           "PYTHONIOENCODING": "utf-8", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    try:
        child = subprocess.run([runtime(project)["python"], "-m", "spatial_collab.power"],
                               input=_json({"project": str(project.root), "spec": spec, "assumptions": assumptions}),
                               text=True, encoding="utf-8", capture_output=True, env=env, timeout=60,
                               creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
    except subprocess.TimeoutExpired as exc:
        raise SpatialError("Power runtime exceeded 60 seconds; any frozen declaration is retained. Reduce the declared simulation budget.") from exc
    if child.returncode:
        raise SpatialError("Power runtime failed; frozen declarations remain: " + child.stderr[-2000:])
    return objects.get(project, json.loads(child.stdout)["object_id"], "powerplan")


if __name__ == "__main__":
    import contextlib
    import json
    import sys
    from .store import Project
    request = json.load(sys.stdin)
    with contextlib.redirect_stdout(sys.stderr):
        plan = create(Project(request["project"]), request["spec"], request["assumptions"])
    print(json.dumps({"object_id": plan["object_id"]}))
