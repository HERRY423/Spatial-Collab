"""Summarize actual saved results; failures and scientific evidence limits stay explicit."""

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from spatial_collab import objects
from spatial_collab.store import Project
from spatial_collab.workflow_validation import validate_result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    p = Project(args.project)
    attempts = [
        json.loads(line)
        for line in (args.output / "executions.jsonl").read_text(encoding="utf-8").splitlines()
    ]
    latest = {job["recipe"]["method"]: job for job in attempts}
    summary = {
        "scope": "public data execution and numerical checks, not independent biological or human-efficiency validation",
        "attempts": len(attempts),
        "failed_attempts": [
            {"id": j["id"], "method": j["recipe"]["method"], "error": j["error"]}
            for j in attempts
            if j["status"] != "succeeded"
        ],
        "methods": {},
    }
    for method, job in latest.items():
        row = {k: job[k] for k in ("status", "result_id", "elapsed_seconds", "error")}
        summary["methods"][method] = row
        if job["status"] != "succeeded":
            continue
        result = objects.get(p, job["result_id"], "analysisresult")
        validate_result(p, result)
        out = result["output"]
        row.update(
            contract_verified=True,
            source_observations=[len(i["observation_ids"]) for i in result["inputs"]],
            implementation_sha256=result["implementation_sha256"],
            package_version=result["environment"].get(
                {
                    "nnls": "scipy",
                    "spatial_spectral": "scikit-learn",
                    "moran_svg": "scipy",
                    "mofa": "mofapy2",
                    "mefisto": "mofapy2",
                    "harmony": "harmonypy",
                    "paste": "paste-bio",
                    "spagcn": "SpaGCN",
                    "spatial_lr": "liana",
                    "progeny": "decoupler",
                    "cell2location": "cell2location",
                }[method]
            ),
        )
        for key in ("embedding", "rna_contribution_fraction", "pathway_activity"):
            if key in out:
                row[key + "_shape"] = list(np.asarray(out[key]).shape)
        if "domains" in out:
            row["domain_sizes"] = dict(Counter(out["domains"]))
        if "rows" in out:
            row["row_statuses"] = dict(Counter(r["status"] for r in out["rows"]))
            row["fdr_below_005"] = sum(
                r.get("q_value") is not None and r["q_value"] < 0.05 for r in out["rows"]
            )
        if method == "nnls":
            row.update(
                cell_types=out["cell_types"],
                signature_condition_number=out["signature_condition_number"],
                median_residual_l2=float(np.median(out["residual_l2"])),
                max_fraction_sum_error=float(
                    np.max(np.abs(np.asarray(out["rna_contribution_fraction"]).sum(axis=1) - 1))
                ),
            )
        if method == "cell2location":
            a = {q: np.asarray(v) for q, v in out["abundance"].items()}
            row.update(
                abundance_shape=list(a["means"].shape),
                median_posterior_interval_width=float(np.median(a["q95"] - a["q05"])),
                training_epochs=out["max_epochs"],
                cells_per_location_prior=result["recipe"]["parameters"]["cells_per_location"],
                convergence_established=False,
            )
        if method == "harmony":
            row["cross_batch_neighbor_fraction"] = out["cross_batch_neighbor_fraction"]
            row["control_only_not_independent_slices"] = True
        if method == "paste":
            target = objects.get(p, result["inputs"][1]["input_id"], "analysisinput")
            if target["provenance"].get("control") == "DERIVED_RIGID_TRANSFORM_NOT_INDEPENDENT_SLICE":
                truth = np.asarray(target["coordinates"])
                row["derived_control_coordinate_rmse"] = float(
                    np.sqrt(
                        np.mean(
                            np.sum((np.asarray(out["target_barycentric_coordinates"]) - truth) ** 2, axis=1)
                        )
                    )
                )
                row["control_only_not_independent_slices"] = True
            row["marginal_max_error"] = out["marginal_max_error"]
        if method == "progeny":
            row.update(
                pathways=out["pathways"],
                coverage=out["coverage"],
                resource_sha256=out["resource_sha256"],
                resource_species_mapping=out["resource_species_mapping"],
            )
        if method in {"mofa", "mefisto"}:
            row.update(
                rna_features=len(out["features"]),
                protein_features=len(out["protein_features"]),
                spatial_covariates_used=out["spatial_covariates_used"],
                iterations_requested=out["iterations_requested"],
            )
    (args.output / "public-validation.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8"
    )
    print(json.dumps({m: r["status"] for m, r in summary["methods"].items()}))


if __name__ == "__main__":
    main()
