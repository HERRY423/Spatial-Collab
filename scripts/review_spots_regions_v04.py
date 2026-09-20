"""Fixed geometry, full RNA region summaries and real duplicate-symbol checks."""
from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path
import time

import h5py
import numpy as np
from anndata.io import read_elem

from spatial_collab.exploration import compare_regions
from spatial_collab.replay import verify_bundle
from spatial_collab.server import ToolService
from spatial_collab.store import SpatialError


ROOT = Path(__file__).resolve().parents[1] / "projects" / "v04-real-data"
SOURCES = Path(r"C:\Plugin\ExoWarrant\output\sc-spatial-multiomics-20260911-01\inputs")
MARKERS = ["Ptprc", "Cd3e", "Ms4a1", "Cd79a", "Ptp4a1", "feature_id:ENSMUSG00000026064", "feature_id:ENSMUSG00000117310"]


def response_bytes(value):
    return len(json.dumps(value, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sample", choices=[1, 2], type=int)
    options = parser.parse_args()
    start = time.perf_counter()
    with h5py.File(SOURCES / f"spleen{options.sample}_adata_RNA.h5ad", "r") as handle:
        xy = np.asarray(read_elem(handle["obsm/spatial"]))
    xmin, ymin = xy.min(axis=0)
    xmax, ymax = np.median(xy[:, 0]), xy[:, 1].max()
    polygon = [[float(xmin), float(ymin)], [float(xmax), float(ymin)],
               [float(xmax), float(ymax)], [float(xmin), float(ymax)]]
    service = ToolService(ROOT / f"spots-spleen-{options.sample}")
    project = service.project
    revision = project.summary()["head_revision"]
    selection = project.set_selection(revision, polygon=polygon,
                                     name="Geometry-only x <= median source x; array scale unverified")
    opened = service.call("open_project")
    opened_size = response_bytes(opened)
    assert "features" not in opened["metadata"] and "panel_genes" not in opened["metadata"]
    del opened
    catalog = service.call("get_feature_catalog", {"query": "Ptp4a1", "limit": 10})
    inspection = service.call("inspect_selection", {"genes": MARKERS})
    assert inspection["expression"]["Ptp4a1"]["status"] == "ambiguous"
    for gene in MARKERS[-2:]:
        assert inspection["expression"][gene]["status"] == "measured"
    try:
        service.call("query_observations", {"revision_id": revision, "gene": "Ptp4a1", "min_count": 1.0})
    except SpatialError as exc:
        ambiguity_error = str(exc)
        assert "Ambiguous" in ambiguity_error
    else:
        raise AssertionError("Duplicate symbols must not silently resolve or aggregate.")
    exact_query = service.call("query_observations", {"revision_id": revision,
        "gene": "feature_id:ENSMUSG00000026064", "min_count": 1.0, "limit": 3})
    result = compare_regions(project, revision, revision, selection["selection_id"], genes=MARKERS)
    assert result["comparison"]["status"] == "unchanged"
    bundle = project.export_bundle(result["run_id"], compact=True)
    del project, service
    gc.collect()
    verification = verify_bundle(bundle["export_path"])
    assert verification["integrity_verified"] and verification["recompute_matches"]
    receipt = {"sample": options.sample, "selection_rule": "All centroids with source x <= source median x; no expression-based selection.",
        "polygon": polygon, "selected_count": selection["cell_count"], "markers_declared": MARKERS,
        "units": "array_index", "source_scope": "All 32285 measured RNA features retained; no cell/gene sampling.",
        "host_response_bytes": {"open_project": opened_size, "inspect_selection": response_bytes(inspection),
            "feature_catalog": response_bytes(catalog), "explicit_feature_query": response_bytes(exact_query)},
        "feature_catalog": catalog, "inspection": inspection, "duplicate_symbol_query_error": ambiguity_error,
        "explicit_feature_query": exact_query, "region_result": result,
        "compact_export": bundle, "bundle_verification": verification,
        "elapsed_seconds": round(time.perf_counter() - start, 3),
        "evidence_ceiling": "Local descriptive RNA/region and tool-service rehearsal. No biological identity adjudication, physical-radius inference, external host adoption or independent scientific validation."}
    path = ROOT / f"spots-spleen-{options.sample}-region-review.json"
    path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"sample": options.sample, "selected_count": selection["cell_count"],
                      "response_bytes": receipt["host_response_bytes"],
                      "foreground_markers": result["before"]["foreground"]["markers"],
                      "recompute_matches": verification["recompute_matches"],
                      "elapsed_seconds": receipt["elapsed_seconds"], "receipt": str(path)}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
