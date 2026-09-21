"""Million-object bounded-memory index and native image-pyramid control."""
import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np
import psutil

from spatial_collab.atlas import build, view, extract
from spatial_collab.demo import create_demo
from spatial_collab.pyramid import register, tile
from spatial_collab.store import Project


def write(path, fields, values):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        writer.writerows(values)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cells", type=int, default=1000000)
    parser.add_argument("--reuse-inputs", action="store_true", help="Reuse explicitly generated raw controls; build a new immutable index")
    args = parser.parse_args()
    args.workspace.mkdir(parents=True, exist_ok=args.reuse_inputs)
    p = Project(args.workspace / "project") if args.reuse_inputs else create_demo(args.workspace / "project")
    n, genes = args.cells, 30001
    start = time.perf_counter()
    if not args.reuse_inputs:
        write(args.workspace / "cells.csv", ["id", "x", "y"], ((str(i), i % 1000, i // 1000) for i in range(n)))
        write(args.workspace / "features.csv", ["id", "symbol"], ((f"G{i}", f"G{i}") for i in range(genes)))
        write(args.workspace / "counts.csv", ["cell_id", "feature", "value"], ((str(i), f"G{(i+j)%genes}", j+1) for i in range(n) for j in range(2)))
        write(args.workspace / "transcripts.csv", ["id", "feature", "x", "y", "z", "qv", "cell_id"], ((str(i), f"G{(i//2)%genes}", (i//2) % 1000 + .1*(i % 2), (i//2)//1000 + .1*(i % 2), i % 3, 30, str(i//2)) for i in range(2*n)))
        write(args.workspace / "boundaries.csv", ["id", "cell_id", "kind", "vertex", "x", "y"], ((f"{i}-{kind}", str(i), kind, j, i % 1000 + x, i // 1000 + y) for i in range(min(n, 1000)) for kind, radius in (("cell", .4), ("nucleus", .2)) for j, (x, y) in enumerate([(-radius, -radius), (radius, -radius), (radius, radius), (-radius, radius)])))
    spec = {"metadata": {"name": "Synthetic million-object scale control", "sample_id": "synthetic-million", "coordinate_system": "known_grid", "units": "micrometer", "platform": "synthetic_scale_control", "observation_unit": "bin", "species": "mouse"}}
    for kind in ("cells", "features", "counts", "transcripts", "boundaries"):
        path = args.workspace / (kind + ".csv")
        with path.open() as f:
            fields = next(csv.reader(f))
        spec[kind] = {"path": str(path.resolve()), "columns": {k: k for k in fields}}
    spec["transcripts"]["projection"] = "xy_projection"
    (args.workspace / "atlas-spec.json").write_text(json.dumps(spec, indent=2))
    generated = time.perf_counter()
    r = build(p, spec)
    imported = time.perf_counter()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.with_suffix(".index.json").write_text(json.dumps({"atlas_id": r["object_id"], "totals": r["totals"], "import_seconds": imported-generated, "database_bytes": r["size"], "scope": "index built; downstream viewport and image validation pending"}, indent=2))
    timings, answers = {}, {}
    for name, kwargs in (("whole_cells", {}), ("whole_transcripts", {"layer": "transcripts"}), ("exact_roi", {"bounds": [-.5, -.5, 31.5, 31.5]}), ("molecule_feature", {"layer": "transcripts", "feature": "G0"}), ("boundaries", {"bounds": [-.5, -.5, 31.5, .5], "layer": "boundaries"})):
        before = time.perf_counter()
        result = view(p, r["object_id"], **kwargs)
        timings[name] = time.perf_counter()-before
        answers[name] = {k: v for k, v in result.items() if k != "records"}
    assert answers["whole_cells"]["total"] == n
    assert answers["whole_transcripts"]["total"] == 2*n
    assert answers["whole_cells"].get("aggregate_total", n) == n
    assert answers["whole_transcripts"].get("aggregate_total", 2*n) == 2*n
    before = time.perf_counter()
    roi = extract(p, r["object_id"], [-.5, -.5, 31.5, 31.5])
    timings["exact_roi_freeze"] = time.perf_counter()-before
    import tifffile
    image = np.indices((4096, 4096), dtype=np.uint16).sum(axis=0).astype(np.uint8)
    image_path = args.workspace / "synthetic.ome.tif"
    if not image_path.exists():
        with tifffile.TiffWriter(image_path, ome=True) as writer:
            writer.write(image, tile=(256, 256), compression="deflate", subifds=4, metadata={"axes": "YX"})
            for step in (2, 4, 8, 16):
                writer.write(image[::step, ::step], tile=(256, 256), compression="deflate", subfiletype=1)
    image_record = register(p, image_path, atlas_id=r["object_id"], pixel_to_world=[[1000/4096, 0, 0], [0, 1000/4096, 0], [0, 0, 1]])
    before = time.perf_counter()
    tiles = [tile(p, image_record["object_id"], level, 0, 0) for level in range(5)]
    timings["five_pyramid_levels"] = time.perf_counter()-before
    memory = psutil.Process().memory_info()
    report = {"scope": "synthetic scale and numerical identity control, not biological or platform validation", "version_unchanged": "0.7.0a1", "project": str(p.root), "atlas_id": r["object_id"], "image_id": image_record["object_id"], "totals": r["totals"], "database_bytes": r["size"], "generation_seconds": generated-start, "import_seconds": imported-generated, "query_seconds": timings, "answers": answers, "frozen_roi_shape": roi["shape"], "pyramid_levels": len(tiles), "peak_working_set_bytes": getattr(memory, "peak_wset", None), "rss_after_bytes": memory.rss, "sources": r["sources"], "database_sha256": r["database_sha256"]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k not in {"answers", "sources"}}, indent=2))


if __name__ == "__main__":
    main()
