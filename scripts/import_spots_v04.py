"""Bounded, read-only-source SPOTS import and exact RNA preservation receipt.

Never assign physical calibration or biological cell labels to these array spots.
Run samples sequentially to keep the working memory bounded.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gc
import hashlib
import json
from pathlib import Path
import time

import h5py
import numpy as np
from anndata.io import read_elem

from spatial_collab.identity import qualify_observation_id
from spatial_collab.importers import import_h5ad, source_digest
from spatial_collab.store import Project, _json


ROOT = Path(__file__).resolve().parents[1]
SOURCES = Path(r"C:\Plugin\ExoWarrant\output\sc-spatial-multiomics-20260911-01\inputs")
OUTPUT = ROOT / "projects" / "v04-real-data"


def preflight(sample: int) -> dict:
    path = SOURCES / f"spleen{sample}_adata_RNA.h5ad"
    namespace = f"spots-spleen-{sample}"
    with h5py.File(path, "r") as handle:
        node = handle["X"]
        obs, var = read_elem(handle["obs"]), read_elem(handle["var"])
        ids, genes, symbols = list(obs.index), var["gene_ids"].tolist(), var.index.tolist()
        xy = np.asarray(read_elem(handle["obsm/spatial"]))
        matrix = read_elem(node).tocsr()
        assert matrix.shape == (len(ids), len(genes))
        assert len(set(ids)) == len(ids) and len(set(genes)) == len(genes)
        assert np.isfinite(matrix.data).all() and (matrix.data >= 0).all()
        assert np.equal(matrix.data, np.floor(matrix.data)).all()
        payload_bytes = 2
        cells_digest = hashlib.sha256()
        for i in range(len(ids)):
            lo, hi = matrix.indptr[i:i+2]
            counts = {genes[int(j)]: float(v) for j, v in zip(matrix.indices[lo:hi], matrix.data[lo:hi]) if v != 0}
            cell = {"cell_id": qualify_observation_id(namespace, ids[i]), "sample_id": namespace,
                    "source_cell_id": ids[i], "x": float(xy[i, 0]), "y": float(xy[i, 1]),
                    "label": "Unannotated", "included": True, "region": "", "counts": counts}
            encoded = _json(cell).encode("utf-8")
            payload_bytes += len(encoded) + (1 if i else 0)
            cells_digest.update(encoded + b"\n")
        features = [{"feature_id": gene, "symbol": symbol} for gene, symbol in zip(genes, symbols)]
        # The small reserved remainder covers provenance, limitations and labels.
        upper_bound = payload_bytes + len(_json(features).encode()) + len(_json(genes).encode()) + 100_000
        return {"sample_id": namespace, "source": source_digest(path), "shape": list(map(int, matrix.shape)),
                "stored_nnz": int(matrix.nnz), "nonzero_count": int(np.count_nonzero(matrix.data)),
                "sum_raw_counts": int(matrix.data.sum(dtype=np.float64)), "unique_feature_ids": len(set(genes)),
                "duplicate_symbol_rows": len(symbols) - len(set(symbols)),
                "ambiguous_symbols": {name: count for name, count in Counter(symbols).items() if count > 1},
                "cell_payload_json_bytes": payload_bytes, "source_snapshot_json_upper_bound_bytes": upper_bound,
                "source_snapshot_file_limit_bytes": 256 * 1024**2,
                "source_order_full_cell_payload_sha256": cells_digest.hexdigest(),
                "units": "array_index", "coordinate_system": "source_spatial_array_xy",
                "observation_unit": "spot", "calibration": "unverified; no physical scale or transform in source"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sample", choices=[1, 2], type=int)
    parser.add_argument("--import-data", action="store_true")
    parser.add_argument("--export-only", action="store_true")
    parser.add_argument("--existing-export", type=Path)
    options = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    receipt_path = OUTPUT / f"spots-spleen-{options.sample}-receipt.json"
    if options.export_only:
        from spatial_collab.replay import verify_bundle
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        project = Project(receipt["project_path"])
        start = time.perf_counter()
        if options.existing_export:
            manifest_path = options.existing_export / "manifest.json"
            bundle = {"export_path": str(options.existing_export), "manifest_path": str(manifest_path),
                      "manifest": json.loads(manifest_path.read_text(encoding="utf-8"))}
        else:
            bundle = project.export_bundle(compact=True)
        del project
        gc.collect()
        result = verify_bundle(bundle["export_path"])
        receipt["compact_export"] = {"bundle": bundle, "verification": result,
                                      "elapsed_seconds": round(time.perf_counter() - start, 3)}
        receipt_path.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(json.dumps(receipt["compact_export"], ensure_ascii=False), flush=True)
        return
    receipt = preflight(options.sample)
    if options.import_data:
        assert receipt["source_snapshot_json_upper_bound_bytes"] < receipt["source_snapshot_file_limit_bytes"]
        start = time.perf_counter()
        destination = OUTPUT / f"spots-spleen-{options.sample}"
        project = import_h5ad(receipt["source"]["path"], destination,
            sample_id=receipt["sample_id"], slice_id=receipt["sample_id"],
            coordinate_system=receipt["coordinate_system"], units=receipt["units"],
            counts_layer="X", allow_unannotated=True, feature_id_key="gene_ids",
            max_features=40_000, max_nnz=10_000_000, platform="SPOTS", observation_unit="spot",
            name=f"SPOTS mouse spleen {options.sample}: full RNA, array coordinates")
        receipt["import_elapsed_seconds"] = round(time.perf_counter() - start, 3)
        receipt["project_path"] = str(destination)
        receipt["project_summary"] = project.summary()
        del project
        gc.collect()
        # Freshly reopen the snapshot, verify every count vector against source,
        # and retain the original raw barcode/coordinate correspondence.
        project = Project(destination)
        cells = project.cells()
        by_original = {cell["source_cell_id"]: cell for cell in cells}
        digest = hashlib.sha256()
        with h5py.File(receipt["source"]["path"], "r") as handle:
            obs = read_elem(handle["obs"])
            for raw_id in obs.index:
                digest.update(_json(by_original[raw_id]).encode("utf-8") + b"\n")
        receipt["imported_source_order_full_cell_payload_sha256"] = digest.hexdigest()
        receipt["full_cell_payload_exact_match"] = digest.hexdigest() == receipt["source_order_full_cell_payload_sha256"]
        assert receipt["full_cell_payload_exact_match"]
        receipt["source_sha256_unchanged"] = source_digest(receipt["source"]["path"])["sha256"] == receipt["source"]["sha256"]
        assert receipt["source_sha256_unchanged"]
        receipt["evidence_ceiling"] = "Local full RNA/coordinate/identity preservation, no biological labels, physical calibration or independent efficacy validation."
    receipt_path.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in receipt.items() if key not in {"project_summary", "ambiguous_symbols"}}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
