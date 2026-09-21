"""Verify full tonsil raw axes/values and replay a real SPOTS integration bundle."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import h5py
import numpy as np

from spatial_collab.proteomics import get_assay, list_assays
from spatial_collab.replay import verify_bundle
from spatial_collab.store import Project


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--tonsil-project", type=Path, required=True)
    parser.add_argument("--replay-project", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    start = time.perf_counter()
    p = Project(args.tonsil_project)
    cells = {c["source_cell_id"]: c for c in p.cells(p.context()["head_revision"])}
    protein = get_assay(p, list_assays(p)["assays"][0]["assay_id"])
    with args.source.open("rb") as f:
        digest = hashlib.file_digest(f, "sha256").hexdigest()
    with h5py.File(args.source) as raw:
        ids, genes, proteins = [raw[k].asstr()[:].tolist() for k in ("cell", "gene", "protein")]
        assert set(cells) == set(ids)
        # Project metadata canonicalizes the panel by ID; never compare a sorted
        # feature axis positionally to the source matrix's original column order.
        assert set(p.summary()["metadata"]["panel_genes"]) == set(genes)
        assert [f["feature_id"] for f in protein["features"]] == proteins
        for i, cid in enumerate(ids):
            cell = cells[cid]
            assert np.array_equal([cell["x"], cell["y"]], raw["pos"][i])
            assert np.array_equal([cell["counts"].get(g, 0) for g in genes], raw["raw_gene_count"][i])
            assert np.array_equal(protein["values"][cell["cell_id"]], raw["raw_protein_count"][i])
        receipt = {"source_sha256": digest, "observations": len(ids), "rna_features": len(genes),
                   "protein_features": len(proteins), "all_raw_ids_coordinates_rna_and_protein_values_equal": True}
    print(json.dumps(receipt), flush=True)
    project = Project(args.replay_project)
    bundle = project.export_bundle(compact=True)
    print("EXPORTED " + bundle["export_path"], flush=True)
    replay = verify_bundle(bundle["export_path"])
    receipt.update(replay=replay, elapsed_seconds=time.perf_counter()-start,
                   evidence_scope="engineering_consistency_not_independent_biology_or_user_efficiency")
    args.output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(replay, ensure_ascii=False), flush=True)
    if not replay["integrity_verified"] or not replay["recompute_matches"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
