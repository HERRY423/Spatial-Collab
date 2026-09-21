"""Portable public-data entrypoint. Download is explicit; existing projects are reused.

Run: python scripts/reproduce_v06.py --workspace projects/v06-public --download
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import time
import urllib.request

import h5py
import numpy as np

from spatial_collab.identity import qualify_observation_id
from spatial_collab.integration import run_baseline, compare_results
from spatial_collab.proteomics import register_protein, list_assays
from spatial_collab.store import Project

COMMIT = "9d59651d78149520661c57bab8a14d42446e6654"
URL = f"https://raw.githubusercontent.com/cmhimself/SMOPCA/{COMMIT}/data/RealDataSample/SpatialCITEseq/Humantonsil_filtered.h5"
SHA256 = "67ee4f78382c53580be49a2ba65e74a3c350834cf7ee08cc3a29e6bd8b11a336"


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def import_tonsil(path, destination):
    if digest(path) != SHA256:
        raise ValueError("Public data checksum does not match the pinned manifest.")
    if (destination / "project.sqlite3").exists():
        p = Project(destination)
        sources = p.summary()["metadata"].get("source_files", [])
        if not any(s.get("sha256") == SHA256 for s in sources):
            raise ValueError("Existing project is not the pinned tonsil source.")
        if list_assays(p)["assays"]:
            return p
        raise ValueError("Existing tonsil project has no protein layer; inspect interrupted import before retry.")
    with h5py.File(path, "r") as data:
        ids = data["cell"].asstr()[:].tolist()
        genes = data["gene"].asstr()[:].tolist()
        proteins = data["protein"].asstr()[:].tolist()
        if len(ids) != len(set(ids)) or len(genes) != len(set(genes)) or len(proteins) != len(set(proteins)):
            raise ValueError("Source axes are ambiguous; no renaming or silent aggregation.")
        xy = np.asarray(data["pos"])
        if data["raw_gene_count"].shape != (len(ids), len(genes)) or data["raw_protein_count"].shape != (len(ids), len(proteins)):
            raise ValueError("Published raw matrices do not match declared axes.")
        sample = "spatial-cite-tonsil"
        cells = []
        for i, cid in enumerate(ids):
            row = data["raw_gene_count"][i]
            if (row < 0).any() or not np.isfinite(row).all() or np.any(row % 1):
                raise ValueError("RNA source is not raw count data.")
            cells.append({"cell_id": qualify_observation_id(sample, cid), "source_cell_id": cid, "sample_id": sample,
                "x": float(xy[i, 0]), "y": float(xy[i, 1]), "label": "Unannotated",
                "counts": {genes[j]: int(row[j]) for j in np.flatnonzero(row)}})
        project = Project.create(destination, cells, {"name": "Public spatial CITE-seq human tonsil",
            "sample_id": sample, "slice_id": sample, "coordinate_system": "published_pos", "units": "array_index",
            "source_kind": "public_spatial_cite_seq", "platform": "spatial-CITE-seq", "observation_unit": "spot",
            "label_semantics": "spot_annotation", "identity_scope": "sample_qualified", "panel_genes": genes,
            "features": [{"feature_id": g, "symbol": g} for g in genes],
            "source_files": [{"path": str(path.resolve()), "sha256": SHA256}], "source_url": URL,
            "biological_replicates": 0, "limitations": ["Author-provided filtered public data, no independent tissue annotation.",
            "Coordinate units are uncalibrated array units; cannot infer micrometers from filenames.",
            "Shared HDF5 cell axis links both raw matrices; not an independent validation of experimental correspondence."]})
        protein_csv = destination / "protein-source.csv"
        with protein_csv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["cell_id", "x", "y", *proteins])
            for i, cid in enumerate(ids):
                writer.writerow([cid, float(xy[i, 0]), float(xy[i, 1]), *data["raw_protein_count"][i].tolist()])
        register_protein(project, protein_csv, sample_id=sample, source_sha256=project.context()["source_sha256"],
            name="Spatial CITE-seq measured ADT", measurement_type="antibody_count", coordinate_system="published_pos",
            units="array_index", registration_note=f"Exact shared cell axis from published HDF5 {SHA256}; raw_protein_count only.",
            format_id="csv", storage="chunked")
    if digest(path) != SHA256:
        raise ValueError("Source changed during import.")
    return project


def evaluate(project, backend="balanced_pca"):
    assays = list_assays(project)["assays"]
    revision = project.context()["head_revision"]
    started = time.perf_counter()
    result = run_baseline(project, revision, assays[0]["assay_id"], components=8, clusters=6, seed=20260920, backend=backend)
    comparison = compare_results(project, result["object_id"], left_partition="rna_only", right_partition="joint", limit=20)
    return {"result_id": result["object_id"], "source_sha256": project.context()["source_sha256"],
            "sample_id": project.summary()["metadata"]["sample_id"], "observations": len(result["observation_ids"]),
            "rna_features_used": len(result["input"]["rna_features"]), "protein_features_used": len(result["input"]["protein_features"]),
            "elapsed_seconds": time.perf_counter() - started, "method": result["method"],
            "rna_joint_comparison": comparison, "evidence_scope": "public_data_engineering_execution_not_independent_biological_or_user_validation"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--backend", choices=["balanced_pca", "smopca"], default="balanced_pca")
    parser.add_argument("--additional-project", action="append", type=Path, default=[])
    args = parser.parse_args()
    args.workspace.mkdir(parents=True, exist_ok=True)
    source = args.source or args.workspace / "Humantonsil_filtered.h5"
    if not source.exists():
        if not args.download:
            parser.error("Supply --source or explicitly use --download for public data.")
        temporary = source.with_suffix(".pending")
        urllib.request.urlretrieve(URL, temporary)
        if digest(temporary) != SHA256:
            raise ValueError("Downloaded file checksum mismatch; not imported.")
        temporary.rename(source)
    project = import_tonsil(source, args.workspace / "tonsil")
    results = []
    for p in [project, *[Project(path) for path in args.additional_project]]:
        print(f"Running {p.root.name}: {args.backend}", flush=True)
        try:
            results.append({"status": "succeeded", **evaluate(p, args.backend)})
        except Exception as exc:
            failure = {"status": "failed", "project": str(p.root), "backend": args.backend,
                       "error_type": type(exc).__name__, "error": str(exc),
                       "timestamp_ns": time.time_ns(), "automatic_fallback": False}
            with (args.workspace / "execution-failures.jsonl").open("a", encoding="utf-8") as stream:
                stream.write(json.dumps(failure, ensure_ascii=False) + "\n")
            results.append(failure)
        (args.workspace / f"validation-{args.backend}.json").write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps([{k: r.get(k) for k in ("status", "sample_id", "observations", "result_id", "elapsed_seconds", "error")} for r in results], indent=2))
    if any(r["status"] == "failed" for r in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
