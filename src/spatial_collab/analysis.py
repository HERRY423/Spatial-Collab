"""Bounded, descriptive sensitivity analysis of immutable annotation revisions.

The metric is a source-cell-edge-weighted neighbor fraction, minus the target
label's abundance in the eligible non-self neighbor universe. It is *not* a
permutation test, evidence of communication, or inference across patients.
"""

from __future__ import annotations

from collections import Counter
from importlib.metadata import PackageNotFoundError, version
import math
from numbers import Real
import platform
from typing import Any

import numpy as np
import scipy
from scipy.spatial import cKDTree

from . import __version__
from .store import SpatialError


# Check project size before materializing cells and graph size before requesting
# neighbor index lists. Only one source's neighbor list is retained at a time.
MAX_PROJECT_CELLS = 250_000
MAX_DIRECTED_SOURCE_EDGES = 5_000_000
MAX_NEIGHBORS_PER_SOURCE = 50_000
MAX_REPORTED_CHANGED_CELLS = 100
ALGORITHM_VERSION = "directed-radius-abundance-sensitivity.v1"


def _runtime_provenance() -> dict:
    try:
        installed_version = version("spatial-collab")
    except PackageNotFoundError:
        installed_version = None
    return {
        "algorithm_version": ALGORITHM_VERSION,
        "python_version": platform.python_version(),
        "software_versions": {
            "spatial_collab_source": __version__,
            "spatial_collab_installed_distribution": installed_version,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "package_installation_status": "installed_distribution_present" if installed_version else "source_checkout_distribution_unknown",
        "runtime_scope": "Versions describe this computation, not independent reproducibility or scientific validation.",
    }


def _text(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > 4000:
        raise SpatialError(f"{name} must be a nonempty string of at most 4000 characters.")
    return value


def _number(value: Any, name: str, *, maximum: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise SpatialError(f"{name} must be a finite number, not a boolean.")
    try:
        value = float(value)
    except (OverflowError, ValueError):
        raise SpatialError(f"{name} must be a finite number.") from None
    if not math.isfinite(value) or value <= 0:
        raise SpatialError(f"{name} must be finite and greater than zero.")
    if maximum is not None and value > maximum:
        raise SpatialError(f"{name} must be at most {maximum}.")
    return value


def _composition(cells: list[dict], selected: set[str]) -> dict:
    selected_cells = [cell for cell in cells if cell["cell_id"] in selected]
    active = [cell for cell in selected_cells if cell["included"]]
    labels = dict(sorted(Counter(cell["label"] for cell in active).items()))
    return {
        "selected_total": len(selected_cells),
        "selected_included": len(active),
        "selected_excluded": len(selected_cells) - len(active),
        "label_counts": labels,
        "label_fractions": {label: count / len(active) for label, count in labels.items()},
        "denominator": "included cells in the frozen selection",
    }


def _evaluate(
    cells: list[dict], selected: set[str], *, radius: float, scope: str,
    source: str, target: str, threshold: float,
) -> dict:
    active = [cell for cell in cells if cell["included"]]
    nodes = active if scope == "whole_slice" else [
        cell for cell in active if cell["cell_id"] in selected
    ]
    sources = [index for index, cell in enumerate(nodes)
               if cell["cell_id"] in selected and cell["label"] == source]
    target_count = sum(cell["label"] == target for cell in nodes)
    eligible_target_count = target_count - int(source == target)
    eligible_neighbor_count = max(0, len(nodes) - 1)
    graph = {
        "node_count": len(nodes),
        "source_count": len(sources),
        "target_count": target_count,
        "eligible_nonself_neighbors_per_source": eligible_neighbor_count,
        "eligible_target_neighbors_per_source": max(0, eligible_target_count),
        "directed_source_edges": 0,
        "target_edges": 0,
        "isolated_sources": len(sources),
        "max_source_degree": 0,
        "edge_definition": "distinct non-self cells at Euclidean distance <= radius_um",
        "edge_denominator": "all directed neighbor edges from included selected source-label cells",
    }
    result = {
        "classification": "inconclusive",
        "reason": None,
        "composition": _composition(cells, selected),
        "graph": graph,
        "observed_fraction": None,
        "null_fraction": None,
        "expected_target_edges": None,
        "excess_over_abundance": None,
    }
    if not sources:
        result["reason"] = "source_label_absent_from_included_selection"
        return result

    coordinates = np.asarray([(cell["x"], cell["y"]) for cell in nodes], dtype=np.float64)
    if not np.isfinite(coordinates).all():
        raise SpatialError("Cell coordinates must be finite.")
    try:
        tree = cKDTree(coordinates)
        # Length-only query avoids potentially quadratic neighbor-list allocation.
        degrees = tree.query_ball_point(
            coordinates[sources], radius, return_length=True, workers=1,
        ) - 1
    except (ValueError, OverflowError) as exc:
        raise SpatialError("Coordinates or radius exceed the supported numeric range.") from exc
    edge_count = int(np.sum(degrees, dtype=np.int64))
    max_degree = int(np.max(degrees)) if len(degrees) else 0
    if edge_count > MAX_DIRECTED_SOURCE_EDGES:
        raise SpatialError(
            f"Graph resource budget exceeded: {edge_count} directed source edges; "
            f"limit {MAX_DIRECTED_SOURCE_EDGES}. Reduce radius or selection."
        )
    if max_degree > MAX_NEIGHBORS_PER_SOURCE:
        raise SpatialError(
            f"Graph resource budget exceeded: one source has {max_degree} neighbors; "
            f"limit {MAX_NEIGHBORS_PER_SOURCE}. Reduce radius or selection."
        )
    graph.update(
        directed_source_edges=edge_count,
        isolated_sources=int(np.count_nonzero(degrees == 0)),
        max_source_degree=max_degree,
    )
    target_flags = np.asarray([cell["label"] == target for cell in nodes], dtype=bool)
    target_edges = 0
    for source_index in sources:
        neighbors = tree.query_ball_point(
            coordinates[source_index], radius, workers=1, return_sorted=True,
        )
        target_edges += sum(bool(target_flags[index]) for index in neighbors if index != source_index)
    graph["target_edges"] = target_edges

    # Missing groups are unavailable evidence, never an observed negative claim.
    if not target_count:
        result["reason"] = "target_label_absent_from_neighbor_universe"
    elif not edge_count:
        result["reason"] = "no_eligible_source_neighbor_edges"
    elif eligible_neighbor_count <= 0 or eligible_target_count <= 0:
        result["reason"] = "nonself_target_abundance_null_unavailable"
    else:
        observed = target_edges / edge_count
        null = eligible_target_count / eligible_neighbor_count
        excess = observed - null
        result.update(
            classification="descriptive_supported" if excess >= threshold else "descriptive_not_supported",
            reason="excess_meets_declared_threshold" if excess >= threshold else "excess_below_declared_threshold",
            observed_fraction=observed,
            null_fraction=null,
            expected_target_edges=edge_count * null,
            excess_over_abundance=excess,
        )
    return result


def _changes(before: list[dict], after: list[dict], selected: set[str]) -> dict:
    previous = {cell["cell_id"]: cell for cell in before}
    totals = Counter()
    selected_totals = Counter()
    examples = []
    changed_count = 0
    for cell in after:
        old = previous[cell["cell_id"]]
        changes = {}
        for field, counter in (("label", "label_changes"), ("included", "exclusion_changes"),
                               ("region", "region_changes")):
            if old.get(field, "") != cell.get(field, ""):
                totals[counter] += 1
                if cell["cell_id"] in selected:
                    selected_totals[counter] += 1
                changes[field] = {"before": old.get(field, ""), "after": cell.get(field, "")}
        if changes:
            changed_count += 1
            if len(examples) < MAX_REPORTED_CHANGED_CELLS:
                examples.append({"cell_id": cell["cell_id"], "in_selection": cell["cell_id"] in selected,
                                 "changes": changes})
    names = ("label_changes", "exclusion_changes", "region_changes")
    return {
        **{name: totals[name] for name in names},
        "selected_changes": {name: selected_totals[name] for name in names},
        "changed_cell_count": changed_count,
        "changed_cells_preview": examples,
        "changed_cells_preview_omitted": max(0, changed_count - len(examples)),
    }


def compare(
    project: Any, base_revision: str, target_revision: str,
    selection_id: str | None = None, radius_um: float = 35.0,
    graph_scope: str = "roi_induced", source_label: str = "T cell",
    target_label: str = "Myeloid", min_effect: float = 0.1,
) -> dict:
    """Save a descriptive before/after comparison using a frozen cell universe.

    ``min_effect`` is an absolute fraction-point threshold for observed neighbor
    fraction minus non-self target abundance. A positive threshold is required.
    It is a researcher-declared descriptive rule, not a significance cutoff.
    """
    _text(base_revision, "base_revision")
    _text(target_revision, "target_revision")
    _text(source_label, "source_label")
    _text(target_label, "target_label")
    if source_label == "Unannotated" or target_label == "Unannotated":
        raise SpatialError("Unannotated is an exploration placeholder, not an established analysis label; inspect and annotate first.")
    radius = _number(radius_um, "radius_um")
    threshold = _number(min_effect, "min_effect", maximum=1.0)
    if not isinstance(graph_scope, str) or graph_scope not in {"roi_induced", "whole_slice"}:
        raise SpatialError("graph_scope must be roi_induced or whole_slice.")
    if selection_id is not None:
        _text(selection_id, "selection_id")

    summary = project.summary()
    if summary["cell_count"] > MAX_PROJECT_CELLS:
        raise SpatialError(
            f"Project resource budget exceeded: {summary['cell_count']} cells; "
            f"limit {MAX_PROJECT_CELLS}. Import a smaller scientifically declared dataset."
        )
    metadata = summary["metadata"]
    if metadata.get("units") != "micrometer":
        raise SpatialError("Radius analysis requires micrometer coordinates.")
    base_info = project.get_revision(base_revision)
    target_info = project.get_revision(target_revision)
    before_cells = sorted(project.cells(base_revision), key=lambda cell: cell["cell_id"])
    after_cells = sorted(project.cells(target_revision), key=lambda cell: cell["cell_id"])
    before_ids = {cell["cell_id"] for cell in before_cells}
    after_ids = {cell["cell_id"] for cell in after_cells}
    if before_ids != after_ids:
        raise SpatialError("Revision cell universes differ; this alpha requires immutable imported cells.")
    for before_cell, after_cell in zip(before_cells, after_cells):
        if before_cell["x"] != after_cell["x"] or before_cell["y"] != after_cell["y"]:
            raise SpatialError("Revision coordinates differ; coordinate changes are outside this analysis contract.")

    selection = project.get_selection(selection_id) if selection_id is not None else None
    if selection is not None:
        if selection["revision_id"] not in {base_revision, target_revision}:
            raise SpatialError("Selection belongs to neither comparison revision; create a fresh selection or compare its revision.")
        for key in ("slice_id", "coordinate_system", "units"):
            if selection.get(key) != metadata.get(key):
                raise SpatialError(f"Selection {key} does not match this project.")
        selected = set(selection["cell_ids"])
        if not selected <= before_ids:
            raise SpatialError("Selection contains unknown cell IDs.")
        selection_record = {
            key: selection.get(key) for key in
            ("selection_id", "revision_id", "slice_id", "coordinate_system", "units", "name", "method", "polygon")
        }
        selection_record["cell_count"] = len(selected)
    elif selection_id is not None:
        raise SpatialError("Unknown selection ID.")
    else:
        selected = before_ids
        selection_record = {"selection_id": None, "method": "all_imported_cells", "cell_count": len(selected)}

    arguments = dict(radius=radius, scope=graph_scope, source=source_label, target=target_label, threshold=threshold)
    before = _evaluate(before_cells, selected, **arguments)
    after = _evaluate(after_cells, selected, **arguments)
    if metadata.get("observation_unit", "cell") != "cell":
        for item in (before, after):
            for key in ("edge_definition", "edge_denominator"):
                item["graph"][key] = item["graph"][key].replace("cells", "observations")
            item["composition"]["denominator"] = "included observations in the frozen selection"
    complete = before["classification"] != "inconclusive" and after["classification"] != "inconclusive"
    comparison = {
        "status": ("stable" if before["classification"] == after["classification"] else "changed") if complete else "indeterminate",
        "effect_delta": after["excess_over_abundance"] - before["excess_over_abundance"] if complete else None,
        "observed_fraction_delta": after["observed_fraction"] - before["observed_fraction"] if complete else None,
        "null_fraction_delta": after["null_fraction"] - before["null_fraction"] if complete else None,
        "selected_cell_count": len(selected),
        "source_edge_denominator_delta": after["graph"]["directed_source_edges"] - before["graph"]["directed_source_edges"],
        "graph_node_count_delta": after["graph"]["node_count"] - before["graph"]["node_count"],
        "source_count_delta": after["graph"]["source_count"] - before["graph"]["source_count"],
        "target_edge_count_delta": after["graph"]["target_edges"] - before["graph"]["target_edges"],
        **_changes(before_cells, after_cells, selected),
    }
    limitations = [
        "Descriptive sensitivity of one slice to specified annotation/inclusion revisions; no p-values or patient-level inference.",
        "Threshold status refers only to the declared excess-over-abundance rule; it does not adjudicate the original biological conclusion.",
        "Baseline abundance is a non-self candidate-label fraction, not a spatial permutation null or control for tissue architecture.",
        "Source-edge weighting gives more weight to high-degree source cells; source-neighbor edges are not independent replicates.",
        "ROI-induced graphs drop neighbors outside the selection; whole-slice graphs retain them. Scope changes the question.",
        "Cell centroids and Euclidean proximity do not establish contact, signaling, causality, or correct segmentation.",
        "A missing source, target, eligible edge, or non-self null makes evidence inconclusive, never a negative claim.",
        "Named reviewer records and computational reproducibility do not establish authenticated review or biological validity.",
    ]
    if metadata.get("source_kind") == "synthetic":
        limitations.insert(0, "SYNTHETIC DATA: demonstration and software verification only; no real-tissue finding.")
    limitations.extend(str(item) for item in metadata.get("limitations", []))
    from .exploration import semantics
    observation_semantics = semantics(metadata)
    if observation_semantics["observation_unit"] != "cell":
        limitations.append("Observations are spots/bins, not single cells; proximity compares annotated capture regions, not cell contacts.")
    result = {
        "analysis_schema": "spatial-collab.descriptive-sensitivity.v1",
        "semantics": observation_semantics,
        "base_revision": base_revision,
        "target_revision": target_revision,
        "selection_id": selection_id,
        "selection": selection_record,
        "parameters": {
            "radius_um": radius, "graph_scope": graph_scope, "source_label": source_label,
            "target_label": target_label, "min_effect": threshold,
            "metric": "directed source-to-target neighbor fraction minus non-self target abundance",
            "threshold_rule": "descriptive_supported iff excess_over_abundance >= min_effect",
        },
        "before": before,
        "after": after,
        "comparison": comparison,
        "evidence_ceiling": "single_slice_descriptive_annotation_sensitivity",
        "scientific_authorization": "NOT_ESTABLISHED",
        "limitations": list(dict.fromkeys(limitations)),
        "data_hashes": {
            "source_sha256": summary.get("source_sha256"),
            "base_revision_sha256": base_info.get("revision_sha256"),
            "target_revision_sha256": target_info.get("revision_sha256"),
        },
        "provenance": _runtime_provenance(),
        "resource_limits": {
            "max_project_cells": MAX_PROJECT_CELLS,
            "max_directed_source_edges": MAX_DIRECTED_SOURCE_EDGES,
            "max_neighbors_per_source": MAX_NEIGHBORS_PER_SOURCE,
        },
    }
    return project.save_run(result)
