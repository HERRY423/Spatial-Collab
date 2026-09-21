"""Explicit synthetic scale check beyond the old 3000-position GP boundary."""
import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np

from spatial_collab.store import Project
from spatial_collab.proteomics import register_protein
from spatial_collab.workflow_inputs import snapshot
from spatial_collab.workflow_methods import run


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(71)
    n = 4000
    xy = np.column_stack([np.arange(n) % 80, np.arange(n) // 80])
    latent = 1 + np.sin(xy[:, 0] / 13) + np.cos(xy[:, 1] / 11)
    raw = rng.poisson(np.exp(latent[:, None] * rng.normal(0, .4, (1, 32))) * 8)
    protein = rng.poisson(np.exp(latent[:, None] * rng.normal(0, .4, (1, 6))) * 12)
    genes = [f"synthetic_gene_{i}" for i in range(32)]
    ids = [f"p{i}" for i in range(n)]
    cells = [{"cell_id": ids[i], "x": float(xy[i, 0]), "y": float(xy[i, 1]), "label": "Unannotated", "counts": dict(zip(genes, raw[i].astype(float).tolist()))} for i in range(n)]
    p = Project.create(args.workspace, cells, {"name": "Synthetic scale control, not biological data", "sample_id": "synthetic-scale", "slice_id": "synthetic-scale", "source_kind": "synthetic", "units": "array_index", "coordinate_system": "synthetic_grid", "panel_genes": genes})
    path = p.root / "synthetic-protein.csv"
    with path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["cell_id", "x", "y", *[f"P{i}" for i in range(6)]])
        for i in range(n):
            writer.writerow([ids[i], *xy[i], *protein[i]])
    assay = register_protein(p, path, sample_id="synthetic-scale", source_sha256=p.context()["source_sha256"], name="Synthetic counts", measurement_type="antibody_count", coordinate_system="synthetic_grid", units="array_index", registration_note="Known synthetic same-object pairing", format_id="csv")
    frozen = snapshot(p, p.context()["head_revision"], "mouse")
    start = time.perf_counter()
    result = run(p, {"method": "mefisto", "input_ids": [frozen["object_id"]], "parameters": {"assay_id": assay["assay_id"], "n_features": 32, "components": 3, "iterations": 20, "frac_inducing": .02}, "seed": 71})
    receipt = {"scope": "SYNTHETIC_SCALE_ONLY_NOT_BIOLOGICAL_VALIDATION", "observations": n, "rna_features": 32, "protein_features": 6, "embedding_shape": list(np.asarray(result["output"]["embedding"]).shape), "spatial_covariates_used": result["output"]["spatial_covariates_used"], "elapsed_seconds": time.perf_counter()-start, "result_id": result["object_id"], "implementation_sha256": result["implementation_sha256"], "environment": result["environment"], "recipe": result["recipe"]}
    args.output.write_text(json.dumps(receipt, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps({k: v for k, v in receipt.items() if k not in {"environment", "recipe"}}))


if __name__ == "__main__":
    main()
