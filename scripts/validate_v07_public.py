"""Real public data execution receipts, with explicit derived regression controls."""

import argparse
import hashlib
import json
from pathlib import Path
import time

import anndata as ad
import numpy as np

from spatial_collab import objects
from spatial_collab.store import Project
from spatial_collab.workflow_inputs import snapshot, register_counts, load_counts
from spatial_collab.workflow_jobs import configure, submit, execute
from spatial_collab.proteomics import list_assays

REFERENCE_SHA = "2098ca88739a9d2d789cb092481c326d7bcbfa73ca815ce3e80b7668cffbe452"


def prepare(p, reference):
    head = p.context()["head_revision"]
    cells = [c for c in p.cells(head) if c["included"]]
    # Explicit geometric ROI chosen before the outcomes, not selected for good performance.
    xy = np.array([[c["x"], c["y"]] for c in cells])
    ranked = np.argsort(np.sum((xy - np.median(xy, axis=0)) ** 2, axis=1), kind="stable")[:300]
    roi = snapshot(p, head, "mouse", [cells[i]["cell_id"] for i in ranked])
    with reference.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    if digest != REFERENCE_SHA:
        raise ValueError("Reference bytes do not match the recorded public source.")
    data = ad.read_h5ad(reference)
    labels = data.obs["cell_types"].astype(str)
    good = (data.obs["tissue"].astype(str) == "Spleen") & ~labels.str.contains(
        "doublet|Low quality|Cycling B/T", case=False, regex=True
    )
    sizes = labels[good].value_counts()
    good &= labels.isin(sizes[sizes >= 5].index)
    # Keep at most 200 actual cells per producer type in original source order.
    selected = np.concatenate(
        [np.flatnonzero(good & (labels == label))[:200] for label in sorted(labels[good].unique())]
    )
    selected.sort()
    ref = register_counts(
        p,
        data.X[selected],
        data.obs_names[selected].tolist(),
        data.var["gene_ids"].astype(str).tolist(),
        feature_symbols=data.var_names.tolist(),
        sample_id="public-sln-spleen-reference",
        species="mouse",
        kind="reference",
        labels=labels.iloc[selected].tolist(),
        provenance={
            "source_url": "https://exampledata.scverse.org/scvi-tools/sln_111.h5ad",
            "file_sha256": digest,
            "scope": "Spleen; exclude producer doublet/low-quality/mixed cycling labels; >=5 cells/type; first <=200 original rows/type",
            "labels": "Original producer cell_types, not independent biological ground truth",
            "feature_mapping": "Original producer var.gene_ids to var index symbols, no inferred mapping",
        },
    )
    record, x = load_counts(p, roi["object_id"])
    # A known transformed copy tests registration mechanics, not cross-slice biological accuracy.
    coordinates = np.asarray(record["coordinates"])
    rotation = np.array([[0.0, -1.0], [1.0, 0.0]])
    moved = register_counts(
        p,
        x,
        ["transformed:" + c for c in record["observation_ids"]],
        record["features"],
        feature_symbols=record.get("feature_symbols"),
        sample_id="spatial-derived-registration-control",
        species="mouse",
        kind="spatial",
        coordinates=coordinates @ rotation + [50, 100],
        provenance={
            "source_input": roi["object_id"],
            "control": "DERIVED_RIGID_TRANSFORM_NOT_INDEPENDENT_SLICE",
            "units": record["provenance"]["units"],
        },
    )
    return roi, ref, moved


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--methods",
        nargs="+",
        default=[
            "mofa",
            "mefisto",
            "nnls",
            "cell2location",
            "harmony",
            "paste",
            "spagcn",
            "spatial_spectral",
            "moran_svg",
            "spatial_lr",
            "progeny",
        ],
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    p = Project(args.project)
    import sys

    configure(p, sys.executable)
    manifest = args.output / "inputs-v2.json"
    if manifest.exists():
        entries = json.loads(manifest.read_text())
        roi, ref, moved = [
            objects.get(p, entries[key], "analysisinput")
            for key in ("roi", "reference", "transformed_control")
        ]
    else:
        roi, ref, moved = prepare(p, args.reference)
        manifest.write_text(
            json.dumps(
                {
                    "roi": roi["object_id"],
                    "reference": ref["object_id"],
                    "transformed_control": moved["object_id"],
                    "reference_shape": ref["shape"],
                    "roi_shape": roi["shape"],
                },
                indent=2,
            )
        )
    assay = list_assays(p)["assays"][0]["assay_id"]
    params = {
        "mofa": {"assay_id": assay, "n_features": 256, "components": 4, "iterations": 100},
        "mefisto": {
            "assay_id": assay,
            "n_features": 128,
            "components": 4,
            "iterations": 60,
            "frac_inducing": 0.1,
        },
        "nnls": {"n_features": 512},
        "cell2location": {
            "n_features": 512,
            "max_epochs": 200,
            "cells_per_location": 8,
            "posterior_samples": 50,
        },
        "harmony": {
            "n_features": 256,
            "components": 8,
            "sample_batches": ["original", "derived_control"],
            "sample_conditions": ["same_source", "same_source"],
        },
        "paste": {"n_features": 128, "iterations": 50, "overlap_assumption": "full_overlap"},
        "spagcn": {"n_features": 128, "components": 8, "clusters": 6, "max_epochs": 100},
        "spatial_spectral": {"n_features": 256, "components": 8, "clusters": 6},
        "moran_svg": {"permutations": 99},
        "spatial_lr": {"permutations": 199},
        "progeny": {"top_targets": 100, "min_targets": 5, "license": "academic"},
    }
    failures = 0
    for method in args.methods:
        inputs = [roi["object_id"]]
        if method in {"nnls", "cell2location"}:
            inputs.append(ref["object_id"])
        if method in {"harmony", "paste"}:
            inputs.append(moved["object_id"])
        spec = {"method": method, "input_ids": inputs, "parameters": params[method], "seed": 20260920}
        print("START", method, flush=True)
        start = time.perf_counter()
        job = submit(p, spec, launch=False, force=True)
        job = execute(p, job["id"])
        receipt = {
            **job,
            "elapsed_seconds": time.perf_counter() - start,
            "validation_scope": "public_count_execution; transformed alignment/batch control not independent sample; cell prior 8 is an explicitly hypothetical sensitivity setting; short training is not convergence or biological validation",
        }
        with (args.output / "executions.jsonl").open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(receipt, ensure_ascii=False) + "\n")
        print(method, job["status"], job.get("error"), flush=True)
        if job["status"] == "succeeded":
            result = objects.get(p, job["result_id"], "analysisresult")
            result_dir = args.output / "results"
            result_dir.mkdir(exist_ok=True)
            (result_dir / (job["result_id"] + ".json")).write_text(
                json.dumps(result, ensure_ascii=False, allow_nan=False), encoding="utf-8"
            )
        failures += job["status"] != "succeeded"
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()
