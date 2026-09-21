"""Study design first: subject-level summaries with explicit ROI comparability."""
from collections import defaultdict

import numpy as np

from . import objects
from .store import SpatialError, _text


def register_study(project, spec):
    _text(spec.get("study_id"), "study_id")
    rows = spec.get("sections")
    if not isinstance(rows, list) or not rows or len(rows) > 10000:
        raise SpatialError("Study needs 1..10000 section records.")
    seen, subjects, batches = set(), defaultdict(set), defaultdict(set)
    for row in rows:
        for key in ("subject_id", "sample_id", "section_id", "condition", "batch", "sample_source", "roi_definition", "roi_class", "source_sha256"):
            _text(row.get(key), key)
        if row["section_id"] in seen:
            raise SpatialError("Study section IDs must be unique.")
        seen.add(row["section_id"])
        if row.get("roi_origin") not in {"prespecified", "posthoc_exploratory", "independent_annotation"}:
            raise SpatialError("Declare ROI origin, including posthoc exploration.")
        if not isinstance(row.get("assay_ids"), list) or not row["assay_ids"] or not isinstance(row.get("run_ids"), list):
            raise SpatialError("Study section requires assay_ids and explicit run_ids (may be empty).")
        if row.get("pair_id") is not None:
            _text(row["pair_id"], "pair_id")
        subjects[row["condition"]].add(row["subject_id"])
        batches[row["batch"]].add(row["condition"])
    conditions = set(subjects)
    confounded = len(conditions) > 1 and all(len(v) == 1 for v in batches.values())
    return objects.put(project, "study", {**spec, "subjects_per_condition": {k: len(v) for k, v in subjects.items()},
        "condition_batch_complete_confounding": confounded,
        "inference_status": "blocked_condition_batch_confounding" if confounded else "design_registered_no_model_run",
        "inference_tool": "run_study_inference",
        "analysis_unit": "subject", "registration_status": "researcher_declared_design_not_verified_sample_identity"})


def compare_study(project, study_id, records):
    design = objects.get(project, study_id, "study")
    sections = {r["section_id"]: r for r in design["sections"]}
    groups = defaultdict(list)
    seen = set()
    for result in records:
        sid = result.get("section_id")
        if sid not in sections or sid in seen:
            raise SpatialError("Unknown or repeated study section result.")
        seen.add(sid)
        row = sections[sid]
        if result.get("source_sha256") != row["source_sha256"] or result.get("run_id") not in row["run_ids"]:
            raise SpatialError("Study summary does not match declared source/run.")
        value = result.get("value")
        if type(value) not in (float, int) or not np.isfinite(value):
            raise SpatialError("Study summary needs a finite value.")
        for key in ("metric", "units", "method_version"):
            _text(result.get(key), key)
        groups[(row["subject_id"], row["condition"], row["roi_class"], result["metric"], result["units"], result["method_version"])].append(value)
    return {"study_id": study_id, "subject_summaries": [
        {"subject_id": key[0], "condition": key[1], "roi_class": key[2], "metric": key[3], "units": key[4],
         "method_version": key[5], "mean_of_section_summaries": float(np.mean(values)), "section_count": len(values)}
        for key, values in sorted(groups.items())],
        "missing_sections": sorted(set(sections) - seen), "condition_batch_complete_confounding": design["condition_batch_complete_confounding"],
        "inference_status": design["inference_status"], "weighting": "Sections equally weighted within each subject; spots never counted as biological replicates."}
