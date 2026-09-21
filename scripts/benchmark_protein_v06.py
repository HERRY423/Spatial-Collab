"""Synthetic scale benchmark; measures engineering, never researcher productivity."""
import argparse
import csv
import json
from pathlib import Path
import statistics
import time

from spatial_collab.proteomics import register_protein, get_assay, _records, inspect_protein
from spatial_collab.store import Project


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("workspace", type=Path)
    args = parser.parse_args()
    n = 12001
    p = Project.create(args.workspace, [{"cell_id": str(i), "x": i % 200, "y": i // 200, "label": "Working", "counts": {"R": 1}} for i in range(n)],
        {"name": "SYNTHETIC scale benchmark", "source_kind": "synthetic", "sample_id": "scale", "slice_id": "scale",
         "coordinate_system": "synthetic_array", "units": "array_index", "panel_genes": ["R"]})
    source = p.root / "protein.csv"
    with source.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cell_id", "x", "y", *[f"P{i}" for i in range(21)]])
        w.writerows([[str(i), i % 200, i // 200, *[(i+j) % 11 for j in range(21)]] for i in range(n)])
    assays = [register_protein(p, source, sample_id="scale", source_sha256=p.context()["source_sha256"], name=f"Assay {i}",
        measurement_type="antibody_count", coordinate_system="synthetic_array", units="array_index", registration_note="Synthetic performance test",
        format_id="csv", storage="chunked") for i in range(6)]
    def timed(fn):
        values = []
        for _ in range(3):
            start = time.perf_counter()
            fn()
            values.append(time.perf_counter()-start)
        return statistics.median(values)
    target = assays[-1]["assay_id"]
    old = timed(lambda: [r for r in _records(p) if r["assay_id"] == target][0])
    new = timed(lambda: get_assay(p, target))
    a = inspect_protein(p, target, p.context()["head_revision"], "P0")
    b = inspect_protein(p, target, p.context()["head_revision"], "P0", offset=a["next_offset"])
    assert len(a["points"]) + len(b["points"]) == n and a["summary"] == b["summary"]
    result = {"input_kind": "synthetic", "observations": n, "protein_features": 21, "assays": 6,
              "all_assay_scan_median_seconds": old, "target_assay_median_seconds": new,
              "pages": [len(a["points"]), len(b["points"])], "exact_statistics_equal_across_pages": True,
              "array_storage_bytes": sum(f.stat().st_size for f in (p.root / "assays").glob("*.npy")),
              "scope": "One machine, warm filesystem, three repetitions; no user-efficiency or broad scalability claim."}
    (p.root / "performance.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
