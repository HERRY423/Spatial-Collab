"""Descriptive inspection of source-reported QC; never an exclusion decision.

Reported Xenium totals include different source feature classes than a selected
gene-expression import. The two denominators are kept distinct. Missing source
fields, explicit nulls, measured zeros and unavailable computed totals differ.
"""
from __future__ import annotations

from collections import Counter
import math

import numpy as np

from .store import SpatialError

NUMERIC_FIELDS = (
    "transcript_counts", "total_counts", "control_probe_counts", "genomic_control_counts",
    "control_codeword_counts", "unassigned_codeword_counts", "deprecated_codeword_counts",
    "cell_area", "nucleus_area", "nucleus_count",
)
MAX_EXAMPLES = 100
MAX_CATEGORIES = 20


def _numeric_summary(values, total, *, absent=0, null=0, invalid_type=0, invalid_range=0, overflow=0):
    count = len(values)
    quantiles = np.quantile(values, [0.25, 0.5, 0.75]).tolist() if values else [None, None, None]
    return {
        "status": ("no_observations" if total == 0 else "available" if count else
                   "unavailable" if overflow or invalid_type or invalid_range else "not_recorded"),
        "observation_count": total, "valid_count": count, "missing_count": absent + null,
        "absent_count": absent, "null_count": null, "invalid_type_count": invalid_type,
        "invalid_range_count": invalid_range, "overflow_count": overflow,
        "zero_count": sum(value == 0 for value in values) if count else None,
        "min": min(values) if count else None, "max": max(values) if count else None,
        "mean": math.fsum(value / count for value in values) if count else None,
        "q25": quantiles[0], "median": quantiles[1], "q75": quantiles[2],
        "quantile_method": "linear interpolation across all valid observed values",
    }


def _recorded_numeric(cells, field):
    values = []
    absent = null = invalid_type = invalid_range = 0
    for cell in cells:
        attributes = cell.get("attributes", {})
        if field not in attributes:
            absent += 1
            continue
        value = attributes[field]
        if value is None:
            null += 1
        elif type(value) not in (int, float):
            invalid_type += 1
        elif value < 0:
            invalid_range += 1
        else:
            values.append(float(value))
    return _numeric_summary(values, len(cells), absent=absent, null=null,
                            invalid_type=invalid_type, invalid_range=invalid_range)


def _segmentation_summary(cells):
    categories = Counter()
    absent = null = invalid_type = 0
    for cell in cells:
        attributes = cell.get("attributes", {})
        if "segmentation_method" not in attributes:
            absent += 1
        elif attributes["segmentation_method"] is None:
            null += 1
        elif not isinstance(attributes["segmentation_method"], str) or not attributes["segmentation_method"].strip():
            invalid_type += 1
        else:
            categories[attributes["segmentation_method"]] += 1
    ordered = sorted(categories.items(), key=lambda item: (-item[1], item[0]))
    shown = ordered[:MAX_CATEGORIES]
    return {"status": "available" if categories else "no_observations" if not cells else "not_recorded",
            "observation_count": len(cells), "valid_count": sum(categories.values()),
            "missing_count": absent + null, "absent_count": absent, "null_count": null,
            "invalid_type_count": invalid_type, "category_count": len(categories),
            "categories": [{"value": key, "count": count} for key, count in shown],
            "omitted_category_count": max(0, len(categories) - MAX_CATEGORIES),
            "omitted_observation_count": sum(count for _, count in ordered[MAX_CATEGORIES:]),
            "interpretation": "Source-reported segmentation method, not a correctness or confidence score."}


def inspect_quality(project, revision_id, selection_id=None):
    """Inspect all frozen selected observations, including excluded ones.

    The method neither creates a run nor changes a selection or annotation. No
    universal QC threshold, doublet score or pass/fail classification is applied.
    """
    if not isinstance(revision_id, str) or not revision_id.strip():
        raise SpatialError("revision_id must identify an exact saved revision.")
    summary = project.summary()
    revision = project.get_revision(revision_id)
    metadata = summary["metadata"]
    cells = project.cells(revision_id)
    selection = None
    if selection_id is not None:
        selection = project.get_selection(selection_id)
        if selection is None or selection["revision_id"] != revision_id:
            raise SpatialError("Quality selection must belong to the exact inspected revision.")
        if selection["source_sha256"] != summary["source_sha256"] or any(
                selection[key] != metadata[key] for key in ("slice_id", "coordinate_system", "units")):
            raise SpatialError("Quality selection does not match the source and coordinate frame.")
        selected = set(selection["cell_ids"])
        if not selected <= {cell["cell_id"] for cell in cells}:
            raise SpatialError("Quality selection contains unknown observation IDs.")
        cells = [cell for cell in cells if cell["cell_id"] in selected]
    panel_count = len(metadata["panel_genes"])
    libraries, detected, zero_examples = [], [], []
    zero_count = overflow = 0
    if panel_count:
        for cell in cells:
            detected.append(sum(value > 0 for value in cell["counts"].values()))
            try:
                library = math.fsum(cell["counts"].values())
            except OverflowError:
                overflow += 1
                continue
            libraries.append(library)
            if library == 0:
                zero_count += 1
                if len(zero_examples) < MAX_EXAMPLES:
                    zero_examples.append({"cell_id": cell["cell_id"], "included": cell["included"],
                                          "label": cell["label"], "imported_panel_count": 0.0})
    imported_library = _numeric_summary(libraries, len(cells),
                                        absent=len(cells) if not panel_count else 0, overflow=overflow)
    detected_features = _numeric_summary(detected, len(cells), absent=len(cells) if not panel_count else 0)
    if not panel_count:
        imported_library["status"] = detected_features["status"] = "no_imported_features"
    recorded = {field: _recorded_numeric(cells, field) for field in NUMERIC_FIELDS}
    return {
        "inspection_schema": "spatial-collab.quality-inspection.v1",
        "revision_id": revision_id, "head_revision": summary["head_revision"],
        "historical": revision_id != summary["head_revision"],
        "source_sha256": summary["source_sha256"], "revision_sha256": revision["revision_sha256"],
        "selection_id": selection_id, "observation_count": len(cells),
        "included_count": sum(cell["included"] for cell in cells),
        "excluded_count": sum(not cell["included"] for cell in cells),
        "observation_unit": metadata.get("observation_unit", "cell"),
        "scope": "all frozen selected observations, including excluded observations",
        "selection_scope": "exact selection" if selection is not None else "all imported observations",
        "recorded_numeric_fields": recorded, "segmentation_method": _segmentation_summary(cells),
        "expression": {"imported_feature_count": panel_count,
                       "library_counts": imported_library, "detected_features": detected_features,
                       "zero_library_observation_count": zero_count if panel_count else None,
                       "zero_library_examples": zero_examples,
                       "zero_library_examples_omitted": max(0, zero_count - len(zero_examples)),
                       "denominator": "only imported measured expression features; controls and omitted genes are not inferred"},
        "changes_applied": False, "selection_changed": False,
        "scientific_authorization": "NOT_ESTABLISHED",
        "limitations": [
            "Source QC values are retained measurements, not verified segmentation quality or scientific confidence.",
            "Absent/null source QC is unknown and is never replaced with zero; malformed types and negative QC values are excluded from numeric summaries and counted separately.",
            "Vendor transcript/total counts and computed imported-panel counts may cover different feature classes; they are not interchangeable denominators.",
            "A zero imported-panel library does not imply no molecules in unimported features or a failed cell.",
            "No calibrated doublet probability, universal QC cutoff, automatic exclusion, or biological validation is provided.",
            "All selected observations, including exclusions, remain in these inspection denominators.",
        ],
    }
