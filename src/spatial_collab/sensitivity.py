"""Non-committing, explicit what-if calculations on a pinned real revision.

Hypotheses are calculation parameters, never researcher-approved revisions.
Every variant is independently applied to the same baseline; none changes head.
"""
from __future__ import annotations

from copy import deepcopy
import math

from .analysis import _evaluate, _number, _runtime_provenance, MAX_PROJECT_CELLS
from .exploration import semantics
from .store import SpatialError, _hash


def run_sensitivity(project, revision_id, selection_id, variants, radii_um,
                    source_label, target_label, graph_scope="whole_slice", min_effect=0.1):
    if not isinstance(variants, list) or not 1 <= len(variants) <= 6:
        raise SpatialError("Provide 1..6 explicit independent variants.")
    if not isinstance(radii_um, list) or not 1 <= len(radii_um) <= 6:
        raise SpatialError("Provide 1..6 declared radii, not an outcome-selected optimum.")
    radii = [_number(r, "radius_um") for r in radii_um]
    if len(set(radii)) != len(radii):
        raise SpatialError("Radii must be distinct.")
    threshold = _number(min_effect, "min_effect", maximum=1.0)
    if graph_scope not in {"roi_induced", "whole_slice"}:
        raise SpatialError("Unknown graph scope.")
    for value in (source_label, target_label):
        if not isinstance(value, str) or not value.strip() or len(value) > 512 or value == "Unannotated":
            raise SpatialError("Source/target need explicit working labels, not Unannotated.")
    summary = project.summary()
    if summary["cell_count"] > MAX_PROJECT_CELLS or summary["metadata"]["units"] != "micrometer":
        raise SpatialError("Sensitivity requires a bounded, calibrated micrometer project.")
    revision = project.get_revision(revision_id)
    selection = project.get_selection(selection_id)
    if not selection or selection["revision_id"] != revision_id:
        raise SpatialError("Selection must be captured at the specified revision.")
    if selection.get("source_sha256") != summary["source_sha256"]:
        raise SpatialError("Selection source hash differs from pinned source.")
    if any(selection.get(k) != summary["metadata"].get(k) for k in ("slice_id", "coordinate_system", "units")):
        raise SpatialError("Selection coordinate frame differs from source.")
    cells = project.cells(revision_id)
    universe = {c["cell_id"] for c in cells}
    selected = set(selection["cell_ids"])
    if not selected or not selected <= universe:
        raise SpatialError("Select a nonempty known observation set.")
    clean, names = [], set()
    for variant in variants:
        if not isinstance(variant, dict) or set(variant) != {"name", "cell_ids", "changes", "rationale"}:
            raise SpatialError("Each variant requires exactly name, cell_ids, changes and rationale.")
        name, ids, changes = variant["name"], variant["cell_ids"], variant["changes"]
        if not isinstance(name, str) or not name.strip() or len(name) > 128 or name in names:
            raise SpatialError("Variant names must be nonempty, distinct and at most 128 characters.")
        if not isinstance(variant["rationale"], str) or not variant["rationale"].strip() or len(variant["rationale"]) > 4000:
            raise SpatialError("Each variant needs a bounded scientific rationale.")
        if (not isinstance(ids, list) or not ids or len(ids) > MAX_PROJECT_CELLS
                or any(not isinstance(x, str) for x in ids) or len(set(ids)) != len(ids) or not set(ids) <= universe):
            raise SpatialError("Variant IDs must be distinct, nonempty and belong to this imported universe.")
        if not isinstance(changes, dict) or not changes or set(changes) - {"label", "included"}:
            raise SpatialError("What-if changes support label and included only.")
        if "included" in changes and type(changes["included"]) is not bool:
            raise SpatialError("included must be a boolean.")
        if "label" in changes and (not isinstance(changes["label"], str) or not changes["label"].strip() or len(changes["label"]) > 512):
            raise SpatialError("label must be a nonempty string of at most 512 characters.")
        names.add(name)
        clean.append({**variant, "cell_ids": sorted(ids)})
    baseline = {str(radius): _evaluate(cells, selected, radius=radius, scope=graph_scope,
                                     source=source_label, target=target_label, threshold=threshold) for radius in radii}
    results = []
    import_scope = summary["metadata"].get("import_scope", {})
    if not isinstance(import_scope, dict):
        import_scope = {}
    bounds = import_scope.get("bounds")
    declared_complete = (import_scope.get("kind") == "spatial_window" and import_scope.get("sampling") == "none"
                         and import_scope.get("selection_rule") == "all cell centroids with xmin <= x <= xmax and ymin <= y <= ymax"
                         and import_scope.get("units") == "micrometer"
                         and import_scope.get("selected_observation_count") == len(cells)
                         and isinstance(bounds, list) and len(bounds) == 4
                         and all(type(v) in (int, float) and math.isfinite(v) for v in bounds)
                         and bounds[0] < bounds[2] and bounds[1] < bounds[3]
                         and all(bounds[0] <= c["x"] <= bounds[2] and bounds[1] <= c["y"] <= bounds[3] for c in cells))
    for variant in clean:
        affected = set(variant["cell_ids"])
        altered = [dict(c, **variant["changes"]) if c["cell_id"] in affected else c for c in cells]
        changed = sum(any(c[k] != v for k, v in variant["changes"].items()) for c in cells if c["cell_id"] in affected)
        for radius in radii:
            before = baseline[str(radius)]
            after = _evaluate(altered, selected, radius=radius, scope=graph_scope,
                              source=source_label, target=target_label, threshold=threshold)
            sufficient = before["classification"] != "inconclusive" and after["classification"] != "inconclusive"
            coverage = "ROI_induced_neighbors_intentionally_truncated" if graph_scope == "roi_induced" else "not_established"
            if declared_complete and graph_scope == "whole_slice":
                candidates = [c for version in (cells, altered) for c in version
                              if c["cell_id"] in selected and c["included"] and c["label"] == source_label]
                if candidates:
                    margin = min(min(c["x"] - bounds[0], bounds[2] - c["x"], c["y"] - bounds[1], bounds[3] - c["y"]) for c in candidates)
                    coverage = "complete_for_radius_within_declared_window" if margin >= radius else "window_boundary_may_truncate_neighbors"
            results.append({"variant": variant["name"], "radius_um": radius, "actual_changed_count": changed,
                            "changed_in_analysis_roi": sum(c["cell_id"] in selected and c["cell_id"] in affected
                                                           and any(c[k] != v for k, v in variant["changes"].items()) for c in cells),
                            "before": before, "after": after, "neighbor_coverage": coverage,
                            "comparison": {"status": ("stable" if before["classification"] == after["classification"] else "changed") if sufficient else "indeterminate",
                                           "effect_delta": after["excess_over_abundance"] - before["excess_over_abundance"] if sufficient else None}})
    # Base and target bind the *same* persistent revision. Counterfactual states
    # are named and hashed parameters; no fabricated revision or reviewer exists.
    result = {"analysis_schema": "spatial-collab.hypothesis-sensitivity.v1",
              "base_revision": revision_id, "target_revision": revision_id, "selection_id": selection_id,
              "selection": {k: selection[k] for k in ("selection_id", "revision_id", "cell_ids", "slice_id", "coordinate_system", "units")},
              "parameters": {"variants": deepcopy(clean), "radii_um": radii, "source_label": source_label,
                             "target_label": target_label, "graph_scope": graph_scope, "min_effect": threshold},
              "variants_sha256": _hash(clean), "results": results, "semantics": semantics(summary["metadata"]),
              "import_scope": import_scope, "revision_created": False, "researcher_approval": "NOT_REQUESTED_HYPOTHESIS_ONLY",
              "evidence_ceiling": "single_slice_counterfactual_descriptive_sensitivity", "scientific_authorization": "NOT_ESTABLISHED",
              "data_hashes": {"source_sha256": summary["source_sha256"], "base_revision_sha256": revision["revision_sha256"],
                              "target_revision_sha256": revision["revision_sha256"]},
              "provenance": {**_runtime_provenance(), "algorithm_version": "independent-hypothesis-grid.v1"},
              "limitations": ["No source, annotation revision or head was changed; this is a computational what-if, not researcher approval.",
                              "All variants independently start from the same revision. Radius/threshold choices must not be optimized then presented as confirmation.",
                              "Whole_slice means the imported observation universe, which may be a crop; background abundance is not full-original-slice abundance.",
                              "Window halo coverage addresses missing radius neighbors only, not tissue sampling bias or the background definition.",
                              "QC attributes and marker coexpression do not establish erroneous segmentation, doublets or true cell identity.",
                              "No significance, independent biological replication or causal inference."]}
    return project.save_run(result)
