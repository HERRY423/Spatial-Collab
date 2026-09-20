"""HER2ST A1 expression/ROI rehearsal with unknown coordinate units.

Only the prepared input H5AD is opened. No sealed truth or private labels.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from anndata.io import read_elem

from spatial_collab.analysis import compare
from spatial_collab.exploration import compare_regions, query_observations
from spatial_collab.identity import qualify_observation_id
from spatial_collab.importers import import_h5ad, source_digest
from spatial_collab.replay import verify_bundle
from spatial_collab.store import SpatialError, _json


ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(r"C:\Spatial Transcriptomics\histoweave-her2st-data\prepared\inputs\A1.h5ad")
DESTINATION = ROOT / "projects" / "v04-real-data" / "her2st-A1"


def main():
    original = source_digest(SOURCE)
    with h5py.File(SOURCE, "r") as handle:
        obs, var = read_elem(handle["obs"]), read_elem(handle["var"])
        matrix = read_elem(handle["layers/counts"]).tocsr()
        xy = np.asarray(read_elem(handle["obsm/spatial"]))
        genes, ids = list(var.index), list(obs.index)
        assert set(obs["sample_id"]) == {"A1"}
        geometry = [float(xy[:, 0].min()), float(xy[:, 1].min()),
                    float(np.median(xy[:, 0])), float(xy[:, 1].max())]
        expected = hashlib.sha256()
        for i, raw_id in enumerate(ids):
            lo, hi = matrix.indptr[i:i+2]
            counts = {genes[int(j)]: float(v) for j, v in zip(matrix.indices[lo:hi], matrix.data[lo:hi]) if v != 0}
            cell = {"cell_id": qualify_observation_id("A1", raw_id), "sample_id": "A1", "source_cell_id": raw_id,
                    "x": float(xy[i, 0]), "y": float(xy[i, 1]), "counts": counts,
                    "label": "Unannotated", "included": True, "region": ""}
            expected.update(_json(cell).encode() + b"\n")
        shape, nnz, total = list(map(int, matrix.shape)), int(matrix.nnz), int(matrix.data.sum(dtype=np.float64))
    project = import_h5ad(SOURCE, DESTINATION, sample_id="A1", slice_id="her2st-A1",
        coordinate_system="prepared_native_xy_unverified", units="unknown", counts_layer="counts",
        allow_unannotated=True, observation_unit="spot", platform="HER2ST",
        name="HER2ST A1: prepared RNA, unverified coordinate scale")
    summary = project.summary()
    actual_cells = {cell["source_cell_id"]: cell for cell in project.cells()}
    actual = hashlib.sha256()
    for raw_id in ids:
        actual.update(_json(actual_cells[raw_id]).encode() + b"\n")
    assert actual.hexdigest() == expected.hexdigest()
    revision = summary["head_revision"]
    xmin, ymin, xmax, ymax = geometry
    selection = project.set_selection(revision, polygon=[[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]],
                                      name="Geometry-only left half by median source x; scale unknown")
    markers = ["ERBB2", "EPCAM", "KRT8", "PTPRC", "MS4A1", "CD79A"]
    inspection = project.inspect_selection(markers)
    result = compare_regions(project, revision, revision, selection["selection_id"], genes=markers)
    assert result["comparison"]["status"] == "unchanged"
    measured = next((gene for gene in markers if gene in genes), None)
    query = query_observations(project, revision, gene=measured, min_count=1 if measured else 0, limit=5)
    try:
        compare(project, revision, revision, selection["selection_id"], radius_um=35)
    except SpatialError as exc:
        radius_error = str(exc)
        assert "micrometer" in radius_error
    else:
        raise AssertionError("Unknown coordinates must not permit micrometre-radius analysis.")
    bundle = project.export_bundle(result["run_id"], compact=True)
    verification = verify_bundle(bundle["export_path"])
    assert verification["integrity_verified"] and verification["recompute_matches"]
    receipt = {"source": original, "source_unchanged": source_digest(SOURCE) == original,
        "project_path": str(DESTINATION), "shape": shape, "raw_nnz": nnz, "raw_count_sum": total,
        "private_truth_opened": False, "labels_inferred": False, "units": "unknown",
        "source_order_full_cell_payload_sha256": expected.hexdigest(), "full_cell_payload_exact_match": True,
        "selection_rule": "source x <= median source x, inclusive bounding polygon; selected from geometry before expression inspection",
        "selection": selection, "inspection": inspection, "region_comparison": result,
        "query_receipt": query, "micrometer_radius_rejection": radius_error,
        "compact_export": bundle, "bundle_verification": verification,
        "evidence_ceiling": "A local expression/ROI software rehearsal on prepared real input, no diagnosis, biological identity claim, physical distance conclusion or independent biological validation."}
    assert receipt["source_unchanged"]
    path = DESTINATION.parent / "her2st-A1-review-receipt.json"
    path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"project": str(DESTINATION), "shape": shape, "nnz": nnz,
                      "selected": selection["cell_count"], "raw_count_sum": total,
                      "source_unchanged": True, "full_payload_match": True,
                      "radius_rejection": radius_error, "recompute_matches": verification["recompute_matches"],
                      "receipt": str(path)}), flush=True)


if __name__ == "__main__":
    main()
