"""Descriptive paired human observations; refuses incomplete or synthetic records."""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
import statistics


def summarize(path):
    groups = defaultdict(dict)
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return {"status": "no_observations", "participants": 0, "efficiency_claim": "NOT_ESTABLISHED"}
    for r in rows:
        if r["reviewer_independence"] != "independent" or not r["result_bundle"]:
            raise ValueError("Require independent grading and a retained result bundle.")
        if r["workflow"] not in {"notebook_viewer", "generic_agent_tools", "agent_spatial_collab"}:
            raise ValueError("Unknown workflow.")
        if r["workflow"] in groups[r["participant_code"]]:
            raise ValueError("Repeated participant/workflow; inspect original observations.")
        for key in ("active_seconds", "waiting_seconds", "total_seconds"):
            value = float(r[key])
            if not 0 <= value < float("inf"):
                raise ValueError("Timings must be finite nonnegative observations.")
            r[key] = value
        if r["total_seconds"] < max(r["active_seconds"], r["waiting_seconds"]):
            raise ValueError("Total time is smaller than a component.")
        for key in ("correctly_completed", "replay_success", "critical_identity_error"):
            if r[key] not in {"0", "1"}:
                raise ValueError("Outcome fields require observed 0/1, never blank defaults.")
        groups[r["participant_code"]][r["workflow"]] = r
    complete = [g for g in groups.values() if len(g) == 3]
    comparisons = {}
    for arm in ("notebook_viewer", "generic_agent_tools"):
        ratios = [g["agent_spatial_collab"]["active_seconds"] / g[arm]["active_seconds"] for g in complete if g[arm]["active_seconds"] > 0]
        comparisons[arm] = {"paired_n": len(ratios), "median_active_time_ratio": statistics.median(ratios) if ratios else None,
                           "raw_ratios": ratios}
    return {"status": "descriptive_pilot", "participants": len(groups), "complete_crossover": len(complete),
            "comparisons": comparisons, "critical_errors": sum(int(r["critical_identity_error"]) for r in rows),
            "correct_completion_by_arm": {a: [int(r["correctly_completed"]) for r in rows if r["workflow"] == a] for a in ("notebook_viewer", "generic_agent_tools", "agent_spatial_collab")},
            "efficiency_claim": "Descriptive observed pilot only; review correctness, order/learning effects and missingness before inference."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("observations")
    print(json.dumps(summarize(parser.parse_args().observations), indent=2))
