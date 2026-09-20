"""Small, revision-pinned research operations shared by all host interfaces."""
from __future__ import annotations

import base64
from collections import Counter
import json
import math

from .analysis import _runtime_provenance
from .store import SpatialError, _hash


def semantics(metadata):
    unit = metadata.get("observation_unit", "cell")
    return {"observation_unit": unit, "platform": metadata.get("platform", "unknown"),
            "label_semantics": metadata.get("label_semantics", "cell_type"),
            "interpretation": ("Working cell labels; centroid proximity does not establish contact."
                               if unit == "cell" else
                               "Spot/bin annotations describe captured regions, not pure cell identities or deconvolution.")}


def query_observations(project, revision_id, *, labels=None, included=None, bounds=None,
                       gene=None, min_count=0.0, limit=100, cursor=None, attribute_filter=None):
    """Explicit paging without sampling or altering the shared selection."""
    if type(limit) is not int or not 1 <= limit <= 500:
        raise SpatialError("limit must be an integer from 1 to 500.")
    if labels is not None and (not isinstance(labels, list) or len(labels) > 100
                               or any(not isinstance(x, str) for x in labels)):
        raise SpatialError("labels must contain at most 100 label strings.")
    if included is not None and type(included) is not bool:
        raise SpatialError("included must be true, false or null.")
    if bounds is not None and (len(bounds) != 4 or any(isinstance(x, bool) or not isinstance(x, (int, float))
                                                     or not math.isfinite(x) for x in bounds)
                               or bounds[0] >= bounds[2] or bounds[1] >= bounds[3]):
        raise SpatialError("bounds must be finite [xmin,ymin,xmax,ymax] with positive area.")
    if isinstance(min_count, bool) or not isinstance(min_count, (int, float)) or not math.isfinite(min_count) or min_count < 0:
        raise SpatialError("min_count must be finite and nonnegative.")
    summary = project.summary()
    project.get_revision(revision_id)
    resolved_gene = gene
    if gene is not None and ("features" in summary["metadata"] or isinstance(gene, str) and gene.startswith(("feature_id:", "symbol:"))):
        from .identity import resolve_feature
        try:
            resolved_gene = resolve_feature(summary["metadata"], gene)
        except ValueError as exc:
            raise SpatialError(str(exc)) from exc
    if gene is not None and resolved_gene not in summary["metadata"]["panel_genes"]:
        raise SpatialError("Requested gene was not measured in the imported feature scope; it cannot be filtered as zero.")
    if gene is None and min_count != 0:
        raise SpatialError("min_count requires a measured gene.")
    if attribute_filter is not None:
        if not isinstance(attribute_filter, dict) or set(attribute_filter) != {"field", "operator", "value"}:
            raise SpatialError("attribute_filter requires exactly field, operator, value.")
        field, op, value = (attribute_filter[k] for k in ("field", "operator", "value"))
        if not isinstance(field, str) or not field.strip() or len(field) > 128 or not isinstance(op, str) or op not in {"eq", "ne", "lt", "le", "gt", "ge"}:
            raise SpatialError("Use a recorded attribute field and eq/ne/lt/le/gt/ge operator.")
        if type(value) not in (str, int, float, bool) or (type(value) is float and not math.isfinite(value)):
            raise SpatialError("Attribute comparisons require a finite number, string or boolean; missing is unknown.")
        if op not in {"eq", "ne"} and type(value) not in (int, float):
            raise SpatialError("Ordered attribute comparisons require a finite number.")
    binding = _hash({"revision_id": revision_id, "source": summary["source_sha256"],
                     "labels": sorted(set(labels)) if labels is not None else None, "included": included,
                     "bounds": bounds, "gene": gene, "min_count": float(min_count), "attribute_filter": attribute_filter})
    offset = 0
    if cursor is not None:
        try:
            if not isinstance(cursor, str) or len(cursor) > 1024:
                raise ValueError()
            decoded = json.loads(base64.urlsafe_b64decode(cursor))
            offset = decoded["offset"]
            if decoded["binding"] != binding or type(offset) is not int or offset < 0:
                raise ValueError()
        except (ValueError, TypeError, KeyError):
            raise SpatialError("Cursor does not match this revision and query; restart paging.") from None
    rows = []
    attribute_unavailable = 0
    for c in project.cells(revision_id):
        if labels is not None and c["label"] not in labels:
            continue
        if included is not None and c["included"] != included:
            continue
        if bounds is not None and not (bounds[0] <= c["x"] <= bounds[2] and bounds[1] <= c["y"] <= bounds[3]):
            continue
        if gene is not None and c["counts"].get(resolved_gene, 0) < min_count:
            continue
        if attribute_filter is not None:
            actual = c.get("attributes", {}).get(field)
            compatible = (type(actual) in (int, float) and type(value) in (int, float)) or type(actual) is type(value)
            if actual is None or not compatible:
                attribute_unavailable += 1
                continue
            matches = {"eq": lambda: actual == value, "ne": lambda: actual != value,
                       "lt": lambda: actual < value, "le": lambda: actual <= value,
                       "gt": lambda: actual > value, "ge": lambda: actual >= value}[op]()
            if not matches:
                continue
        item = {k: c[k] for k in ("cell_id", "x", "y", "label", "included", "region")}
        if "sample_id" in c:
            item.update(sample_id=c["sample_id"], source_cell_id=c["source_cell_id"])
        if gene is not None:
            item["marker_count"] = c["counts"].get(resolved_gene, 0)
            if "features" in summary["metadata"]:
                item["feature_id"] = resolved_gene
        if attribute_filter is not None:
            item["matched_attribute"] = {field: actual}
        rows.append(item)
    rows.sort(key=lambda c: c["cell_id"])
    if offset > len(rows):
        raise SpatialError("Cursor offset exceeds matching observations.")
    page = rows[offset:offset + limit]
    next_offset = offset + len(page)
    next_cursor = (base64.urlsafe_b64encode(json.dumps({"binding": binding, "offset": next_offset}).encode()).decode()
                   if next_offset < len(rows) else None)
    return {"revision_id": revision_id, "head_revision": summary["head_revision"],
            "historical": revision_id != summary["head_revision"], "source_sha256": summary["source_sha256"],
            "semantics": semantics(summary["metadata"]), "observations": page, "matching_count": len(rows),
            "returned_count": len(page), "next_cursor": next_cursor, "complete": next_cursor is None,
            "selection_changed": False, "attribute_filter": attribute_filter,
            "attribute_unavailable_count": attribute_unavailable,
            "attribute_interpretation": "Missing, null and incompatible values are unknown, never zero or automatically unequal; count applies after other filters."}


def _region_summary(cells, ids, genes, panel, resolutions=None):
    selected = [c for c in cells if c["cell_id"] in ids]
    active = [c for c in selected if c["included"]]
    labels = dict(sorted(Counter(c["label"] for c in active).items()))
    markers = {}
    libraries = []
    for c in active:
        largest = max(c["counts"].values(), default=0)
        if largest > 0:
            libraries.append((c, largest, math.fsum(v / largest for v in c["counts"].values())))
    for gene in genes:
        resolution = resolutions.get(gene) if resolutions else None
        feature_id = resolution["feature_id"] if resolution else gene
        if feature_id not in panel:
            markers[gene] = {"status": "unmeasured", "mean_raw_count": None, "detection_fraction": None,
                             "mean_panel_counts_per_10000": None}
            if resolution:
                markers[gene].update(status=resolution["status"], feature_resolution=resolution)
            continue
        values = [c["counts"].get(feature_id, 0) for c in active]
        # Scale before summing to avoid overflow from individually finite values.
        normalized = [(c["counts"].get(feature_id, 0) / largest) / denominator * 10000 for c, largest, denominator in libraries]
        markers[gene] = {"status": "measured" if active else "empty_region",
                         "mean_raw_count": math.fsum(v / len(values) for v in values) if values else None,
                         "detection_fraction": sum(v > 0 for v in values) / len(values) if values else None,
                         "mean_panel_counts_per_10000": math.fsum(normalized) / len(normalized) if normalized else None,
                         "normalized_denominator": len(normalized), "zero_library_observations": len(active) - len(normalized)}
        if resolution:
            markers[gene]["feature_resolution"] = resolution
    return {"selected_count": len(selected), "included_count": len(active),
            "excluded_count": len(selected) - len(active), "label_counts": labels,
            "label_fractions": {k: v / len(active) for k, v in labels.items()}, "markers": markers}


def compare_regions(project, base_revision, target_revision, selection_id, *, background_selection_id=None, genes=None):
    """Frozen foreground and disjoint background; descriptive changes only."""
    genes = [] if genes is None else genes
    if not isinstance(genes, list) or len(genes) > 100 or any(not isinstance(g, str) or not g.strip() for g in genes):
        raise SpatialError("genes must contain at most 100 nonempty names.")
    genes = sorted(set(genes))
    summary = project.summary()
    meta = summary["metadata"]
    resolutions = None
    if "features" in meta or any(gene.startswith(("feature_id:", "symbol:")) for gene in genes):
        from .identity import describe_feature_query
        resolutions = {gene: describe_feature_query(meta, gene) for gene in genes}
    before_rev, after_rev = project.get_revision(base_revision), project.get_revision(target_revision)
    before, after = project.cells(base_revision), project.cells(target_revision)
    universe = {c["cell_id"] for c in before}
    if universe != {c["cell_id"] for c in after}:
        raise SpatialError("Revision observation universes differ.")
    def selection(sid):
        record = project.get_selection(sid)
        if record is None or record["revision_id"] not in {base_revision, target_revision}:
            raise SpatialError("Each selection must belong to either comparison revision.")
        if any(record.get(k) != meta.get(k) for k in ("slice_id", "coordinate_system", "units")):
            raise SpatialError("Selection coordinate frame differs from project.")
        ids = set(record["cell_ids"])
        if not ids or not ids <= universe:
            raise SpatialError("Selections must contain nonempty known observation IDs.")
        return record, ids
    foreground, fg_ids = selection(selection_id)
    if background_selection_id is not None:
        background, bg_ids = selection(background_selection_id)
        background_name = background["name"]
    else:
        bg_ids = universe - fg_ids
        background_name = "All other imported observations (possibly a window, not the complete original slice)"
    if fg_ids & bg_ids:
        raise SpatialError("Foreground and background must be disjoint; overlapping denominators are not silently adjusted.")
    if not bg_ids:
        raise SpatialError("Background is empty; select a smaller foreground or explicit disjoint background.")
    def evaluate(cells):
        fg, bg = (_region_summary(cells, ids, genes, set(meta["panel_genes"]), resolutions) for ids in (fg_ids, bg_ids))
        enough = bool(fg["included_count"] and bg["included_count"])
        labels = sorted(set(fg["label_counts"]) | set(bg["label_counts"]))
        deltas = {label: (fg["label_fractions"].get(label, 0) - bg["label_fractions"].get(label, 0))
                  if enough else None for label in labels}
        markers = {}
        for gene in genes:
            a, b = fg["markers"][gene], bg["markers"][gene]
            markers[gene] = {key + "_difference": a[key] - b[key] if a[key] is not None and b[key] is not None else None
                             for key in ("mean_raw_count", "detection_fraction", "mean_panel_counts_per_10000")}
        return {"status": "descriptive" if enough else "indeterminate", "foreground": fg, "background": bg,
                "label_fraction_difference": deltas, "marker_contrast": markers}
    a, b = evaluate(before), evaluate(after)
    label_union = sorted(set(a["label_fraction_difference"]) | set(b["label_fraction_difference"]))
    enough = a["status"] == b["status"] == "descriptive"
    result = {"analysis_schema": "spatial-collab.region-contrast.v1", "base_revision": base_revision,
              "target_revision": target_revision, "selection_id": selection_id,
              "selection": {"foreground_selection_id": selection_id, "foreground_name": foreground["name"],
                            "foreground_count": len(fg_ids), "background_selection_id": background_selection_id,
                            "background_name": background_name, "background_count": len(bg_ids), "overlap_count": 0},
              "parameters": {"background_selection_id": background_selection_id, "genes": genes,
                             "normalization": "counts per 10000 measured panel counts; zero-library observations excluded only from normalized mean"},
              "semantics": semantics(meta), "before": a, "after": b,
              "comparison": {"status": ("unchanged" if a == b else "descriptive_values_changed") if enough else "indeterminate",
                             "label_contrast_delta": {label: b["label_fraction_difference"].get(label, 0) - a["label_fraction_difference"].get(label, 0)
                                                      if enough else None for label in label_union}},
              "data_hashes": {"source_sha256": summary["source_sha256"], "base_revision_sha256": before_rev["revision_sha256"],
                              "target_revision_sha256": after_rev["revision_sha256"]},
              "evidence_ceiling": "single_slice_descriptive_region_contrast", "scientific_authorization": "NOT_ESTABLISHED",
              "provenance": {**_runtime_provenance(), "algorithm_version": "disjoint-region-contrast.v1"},
              "limitations": ["Selected objects and background are frozen across revisions; excluded observations leave each denominator.",
                              "Most abundant does not mean most enriched. Differences are descriptive fractions, not enrichment significance.",
                              "Normalization uses only imported measured features, not full-transcriptome library size unless those coincide.",
                              "Unmeasured genes remain unknown. Marker contrast is not differential expression or a cell identity verdict.",
                              "One slice; objects and ROIs are not biological replicates. No patient, causal or segmentation-validity inference.",
                              *meta.get("limitations", [])]}
    return project.save_run(result)
