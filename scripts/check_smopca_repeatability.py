"""Report, never hide, cross-build and within-build SMOPCA output differences."""
import argparse
import json
from pathlib import Path

import numpy as np

from spatial_collab import objects
from spatial_collab.integration import run_baseline, compare_results
from spatial_collab.store import Project


def delta(a, b):
    results = {}
    for name in a["representations"]:
        x, y = np.asarray(a["representations"][name]), np.asarray(b["representations"][name])
        signs = np.sign(np.sum(x*y, axis=0))
        signs[signs == 0] = 1
        results[name] = {"sign_aligned_max_abs_difference": float(np.abs(x-y*signs).max()),
                         "embedding_allclose": bool(np.allclose(x, y*signs, rtol=1e-7, atol=1e-9)),
                         "domain_ids_equal": a["domains"][name] == b["domains"][name]}
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("project", type=Path)
    parser.add_argument("reference_result")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    p = Project(args.project)
    original = objects.get(p, args.reference_result, "integration")
    identity = original["input"]
    parameters = original["method"]["parameters"]
    fitted = []
    for i in range(2):
        print(f"RUN {i+1}", flush=True)
        fitted.append(run_baseline(p, identity["fit_revision"], identity["protein_assay_id"],
            rna_features=identity["rna_features"], protein_features=identity["protein_features"],
            components=parameters["components"], clusters=parameters["clusters"], seed=original["method"]["seed"],
            backend="smopca", protein_transform=original["preprocessing"]["protein"]))
    report = {"reference_id": original["object_id"], "new_result_ids": [r["object_id"] for r in fitted],
              "reference_environment": original["method"]["environment"], "environment": fitted[0]["method"]["environment"],
              "cross_build": delta(original, fitted[0]), "same_build": delta(*fitted),
              "cross_build_joint_comparison": compare_results(p, original["object_id"], fitted[0]["object_id"], "joint", "joint", limit=5)}
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"cross_build": report["cross_build"], "same_build": report["same_build"]}), flush=True)
    if not all(r["embedding_allclose"] and r["domain_ids_equal"] for r in report["same_build"].values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
