"""Retrospective public-input design checks; never relabel prior data as preregistered."""
import json
from pathlib import Path

from spatial_collab import objects, power
from spatial_collab.store import Project
from spatial_collab.workflow_methods import run


def main():
    root = Path(__file__).resolve().parents[1]
    project = Project(root / "projects/v05-real-data/spots-spleen-1")
    previous = json.loads((root / "output/v07/public-validation.json").read_text())
    output = root / "output/power"
    output.mkdir(exist_ok=True)
    report = {"scope": "Retrospective graph-conditional planning on previously analyzed public SPOTS ROI; no new biological validation", "methods": {}}
    for method in ("moran_svg", "spatial_lr"):
        old = objects.get(project, previous["methods"][method]["result_id"], "analysisresult")
        spec = old["recipe"]
        base = {"model": "gaussian_sar" if method == "moran_svg" else "lognormal_neighbor_coupling",
                "simulations": 200, "target_power": .8, "effect_grid": [0, .2, .4, .6, .8, .95], "seed": 731,
                "rationale": "Retrospective engineering comparison of sparse versus dense BH rank assumptions. Dense scenario is not a prediction of biological signal abundance."}
        sparse_plan = power.create(project, spec, {**base, "bh_rank": 1})
        rank = sparse_plan["result"]["resolution"]["minimum_possible_bh_rejection_rank"]
        dense_plan = power.create(project, spec, {**base, "bh_rank": rank})
        assert sparse_plan["timing"] == dense_plan["timing"] == "retrospective_sensitivity"
        result = run(project, {**spec, "power_plan_id": sparse_plan["object_id"]})
        record = {"source_input": old["inputs"][0], "previous_result_id": old["object_id"],
                  "sparse_plan": sparse_plan, "dense_rank_scenario": dense_plan,
                  "new_result_id": result["object_id"], "discovery_summary": result["output"]["discovery_summary"],
                  "method_contract": result["method_contract"], "execution_receipt": result["execution_receipt"]}
        report["methods"][method] = record
        (output / "public-power-validation.json").write_text(json.dumps(report, indent=2))
        print(method, sparse_plan["resolved"]["family_size"], "rank", rank,
              "conditional_mde", dense_plan["result"]["mde_rho"], flush=True)


if __name__ == "__main__":
    main()
