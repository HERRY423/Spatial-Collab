"""Frozen real-data rehearsal, with explicit automated working-label inputs.

Run from repository root with --input-project, --output and --plan. Outputs must
not already exist. No raw data or existing project is edited. No reviewer or
annotation approval is invented. The scientific plan is read, hashed and kept.
"""
from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
import time

from spatial_collab.replay import verify_bundle
from spatial_collab.server import ToolService
from spatial_collab.store import Project, SpatialError


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-project", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SpatialError("Output exists; use a new review directory to preserve prior evidence.")
    plan = json.loads(args.plan.read_text(encoding="utf-8"))
    source = Project(args.input_project)
    original = source.summary()
    # A previous rehearsal exposed a wrongly supplied tissue identifier. This
    # study refuses that project rather than patching its immutable metadata.
    if original["metadata"]["slice_id"] != "xenium-prime-reactive-lymph-node":
        raise SpatialError("Study requires the provenance-checked reactive lymph-node slice identity.")
    scope = original["metadata"].get("import_scope", {})
    if scope.get("bounds") != plan["import_window_um"] or scope.get("sampling") != "none":
        raise SpatialError("Imported window differs from the frozen geometry-only study plan.")
    rules = plan["marker_rule"]
    cells = source.cells(original["head_revision"])
    panel = set(original["metadata"]["panel_genes"])
    measured_markers = rules["B_program"] + rules["T_program"]
    if not set(measured_markers) <= panel:
        raise SpatialError("A declared working-group marker was not measured; do not turn missing into zero.")
    if any(c["label"] != "Unannotated" for c in cells):
        raise SpatialError("Expected native unannotated baseline; avoid silently replacing existing annotations.")
    for cell in cells:
        supports = [sum(cell["counts"].get(g, 0) >= rules["minimum_raw_count_per_positive_marker"] for g in rules[k])
                    for k in ("B_program", "T_program")]
        b, t = [n >= rules["minimum_positive_markers"] for n in supports]
        cell["label"] = rules["both_programs"] if b and t else "B_program" if b else "T_program" if t else rules["neither"]
        cell.setdefault("attributes", {}).update(b_marker_support=supports[0], t_marker_support=supports[1])
    meta = deepcopy(original["metadata"])
    meta.update(name="Xenium reactive lymph node | declared marker-program hypotheses", biological_replicates=1)
    meta["working_annotation"] = {"kind": "automated_exploratory_marker_rule", "rules": rules,
                                  "source_project_sha256": original["source_sha256"],
                                  "source_revision": original["head_revision"],
                                  "plan_sha256": hashlib.sha256(args.plan.read_bytes()).hexdigest(),
                                  "researcher_approval": "NOT_ESTABLISHED", "validated_cell_identity": False}
    meta["limitations"] += ["Working labels are deterministic marker programs, not expert cell identities or label-transfer truth.",
                            "A single fixed continuous window and one donor; no whole-tissue representativeness or biological replication."]
    args.output.mkdir(parents=True)
    write(args.output / "frozen-plan.json", plan)
    project = Project.create(args.output / "working-project", cells, meta)
    service = ToolService(project.root)
    revision = project.summary()["head_revision"]
    timings = []
    def call(operation, **arguments):
        start = time.perf_counter()
        result = service.call(operation, arguments)
        timings.append({"operation": operation, "seconds": time.perf_counter() - start})
        print(json.dumps({"completed": operation, "seconds": round(timings[-1]["seconds"], 3)}), flush=True)
        return result
    view = call("open_project")
    assert view["view"]["complete"] and view["view"]["returned_cells"] == len(cells)
    context = call("get_context")
    assert not call("get_context", after_cursor=context["cursor"])["changed"]
    xmin, ymin, xmax, ymax = plan["analysis_roi_um"]
    polygon = [[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]]
    roi = call("set_selection", expected_revision=revision, polygon=polygon, name="Frozen 400 um ROI; 100 um halo")
    roi_ids = set(roi["cell_ids"])
    write(args.output / "quality-full-window.json", call("inspect_quality", revision_id=revision))
    write(args.output / "quality-roi.json", call("inspect_quality", revision_id=revision, selection_id=roi["selection_id"]))
    markers = call("inspect_selection", genes=measured_markers + plan["unavailable_usual_markers"])
    write(args.output / "marker-availability.json", markers)
    for gene in plan["unavailable_usual_markers"]:
        assert markers["expression"][gene]["status"] == "unmeasured"
    # Exercise agent-visible exact-ID discovery rather than private ad hoc QC
    # selection. Null attributes must never be interpreted as zero or unequal.
    qc_ids, cursor = [], None
    while True:
        page = call("query_observations", revision_id=revision,
                    attribute_filter={"field": "nucleus_count", "operator": "ne", "value": 1}, limit=500, cursor=cursor)
        qc_ids.extend(row["cell_id"] for row in page["observations"])
        cursor = page["next_cursor"]
        if cursor is None:
            break
    expected_qc = {c["cell_id"] for c in cells if c["attributes"].get("nucleus_count") is not None and c["attributes"]["nucleus_count"] != 1}
    assert set(qc_ids) == expected_qc and len(qc_ids) == len(expected_qc)
    weak_ids = [c["cell_id"] for c in cells if (c["label"] == "B_program" and c["attributes"]["b_marker_support"] == 2)
                or (c["label"] == "T_program" and c["attributes"]["t_marker_support"] == 2)]
    variants = [{"name": "non_single_nucleus_uncertain", "cell_ids": qc_ids, "changes": {"label": "Uncertain"}, "rationale": plan["variants"][0]},
                {"name": "non_single_nucleus_excluded", "cell_ids": qc_ids, "changes": {"included": False}, "rationale": plan["variants"][1]},
                {"name": "require_three_markers", "cell_ids": weak_ids, "changes": {"label": "Uncertain"},
                 "rationale": plan["variants"][2] + "; applies to B_program/T_program single-program groups only; Mixed_program remains unchanged"}]
    write(args.output / "declared-variants.json", variants)
    # This is a computational hypothesis rehearsal. No proposal is committed.
    results = {}
    for graph_scope in ("whole_slice", "roi_induced"):
        results[graph_scope] = call("run_sensitivity", revision_id=revision, selection_id=roi["selection_id"],
                                    variants=variants, radii_um=[float(r) for r in plan["radii_um"]],
                                    source_label="T_program", target_label="B_program", graph_scope=graph_scope,
                                    min_effect=plan["min_effect"])
        write(args.output / f"sensitivity-{graph_scope}.json", results[graph_scope])
    fg_ids = [c["cell_id"] for c in cells if c["cell_id"] in roi_ids and c["label"] == "B_program"]
    background = call("set_selection", expected_revision=revision, cell_ids=sorted(roi_ids - set(fg_ids)), name="Other working groups within the same frozen ROI")
    foreground = call("set_selection", expected_revision=revision, cell_ids=fg_ids, name="B-program observations in frozen ROI")
    region = call("run_region_comparison", base_revision=revision, target_revision=revision,
                  selection_id=foreground["selection_id"], background_selection_id=background["selection_id"],
                  genes=measured_markers + ["MZB1", "MARCO", "PECAM1"] + plan["unavailable_usual_markers"])
    write(args.output / "region-description.json", region)
    # Restore the actual analysis ROI in the shared workbench, without labels.
    call("set_selection", expected_revision=revision, polygon=polygon, name="Frozen analysis ROI | working hypotheses only")
    exports = []
    for graph_scope, result in results.items():
        exported = call("export_review_bundle", run_id=result["run_id"])
        start = time.perf_counter()
        verification = verify_bundle(exported["export_path"])
        timings.append({"operation": "offline_replay_" + graph_scope, "seconds": time.perf_counter() - start})
        assert verification["recompute_matches"]
        exports.append({"scope": graph_scope, "export": exported, "verification": verification})
        print(json.dumps({"replay": graph_scope, "matches": True}), flush=True)
    write(args.output / "exports-and-replay.json", exports)
    rows = []
    for graph_scope, result in results.items():
        for entry in result["results"]:
            rows.append({"graph_scope": graph_scope, "variant": entry["variant"], "radius_um": entry["radius_um"],
                         "changed_all": entry["actual_changed_count"], "changed_roi": entry["changed_in_analysis_roi"],
                         "before_effect": entry["before"]["excess_over_abundance"], "after_effect": entry["after"]["excess_over_abundance"],
                         "before_classification": entry["before"]["classification"], "after_classification": entry["after"]["classification"],
                         "effect_delta": entry["comparison"]["effect_delta"], "status": entry["comparison"]["status"],
                         "neighbor_coverage": entry["neighbor_coverage"]})
    with (args.output / "sensitivity-summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    final = project.summary()
    assert final["head_revision"] == revision and final["revision_count"] == 1
    assert source.summary()["source_sha256"] == original["source_sha256"]
    assert source.summary()["head_revision"] == original["head_revision"]
    receipt = {"study_id": plan["study_id"], "plan_sha256": meta["working_annotation"]["plan_sha256"],
               "raw_project": str(source.root), "working_project": str(project.root),
               "source_sha256": original["source_sha256"], "working_sha256": final["source_sha256"],
               "source_files": meta["source_files"], "working_annotation": meta["working_annotation"],
               "window_cells": len(cells), "roi_cells": len(roi_ids), "gene_count": len(panel),
               "working_groups_window": dict(Counter(c["label"] for c in cells)),
               "working_groups_roi": dict(Counter(c["label"] for c in cells if c["cell_id"] in roi_ids)),
               "variant_counts": [{"name": v["name"], "all": len(v["cell_ids"]), "roi": len(set(v["cell_ids"]) & roi_ids)} for v in variants],
               "sensitivity_results": rows, "run_ids": {k: v["run_id"] for k, v in results.items()},
               "region_run_id": region["run_id"], "head_unchanged": True, "persistent_annotation_revisions_created": 0,
               "both_sensitivity_exports_recompute": True, "ui_full_window_complete": view["view"]["complete"],
               "timings": timings, "scientific_authorization": "NOT_ESTABLISHED",
               "limitations": plan["out_of_scope"]}
    write(args.output / "study-receipt.json", receipt)
    print(json.dumps({"receipt": str(args.output / "study-receipt.json"), "window_cells": len(cells), "roi_cells": len(roi_ids),
                      "group_counts": receipt["working_groups_roi"], "head_unchanged": True}), flush=True)


if __name__ == "__main__":
    main()
