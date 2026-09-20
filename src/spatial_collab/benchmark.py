"""Developer-authored synthetic task replay; never a model or scientific benchmark score.

Fixtures and their expected answers are public. The scripted reference driver
exercises persisted product operations without network calls or an agent model.
"""
from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import time
import uuid

from . import analysis
from .replay import verify_bundle
from .store import Project, SpatialError


SUITE_VERSION = "spatial-collab.synthetic-collaboration.v2"
EVIDENCE_CEILING = "developer_authored_synthetic_reference_replay"

_TASKS = (
    ("panel_evidence", "cell_typing",
     "For the selected T-cell candidate, distinguish a measured marker with zero counts from a marker "
     "outside the panel. Keep the sparse zeros and missing measurements distinct.",
     {"measured_zero_status": "measured", "measured_zero_mean": 0,
      "outside_panel_status": "unmeasured", "outside_panel_mean": None, "selected_count": 1}),
    ("spatial_object_identity", "spatial_analysis",
     "Identify the observations whose centroids lie in the supplied polygon, including its boundary. "
     "Preserve the exact identifiers, physical units, coordinate system and annotation version.",
     {"selected_ids": ["near-other", "near-target", "source"], "units": "micrometer",
      "coordinate_system": "synthetic_xy", "revision_matches": True, "source_matches": True}),
    ("annotation_sensitivity", "cell_typing",
     "A researcher changes the near-target label from Myeloid to Stromal. Determine whether the "
     "original descriptive neighbor result changes, keeping its radius and observation universe fixed. "
     "A separate region-name-only revision should leave that metric unchanged.",
     {"before_effect": 0.5, "after_effect": -0.25, "status": "changed",
      "effect_delta": -0.75, "label_changes": 1, "region_only_status": "stable",
      "region_only_effect_delta": 0, "old_result_historical": True}),
    ("neighbor_context", "spatial_analysis",
     "For the ROI containing source and near-target, compare an ROI-induced graph with one retaining "
     "outside-ROI observations as potential neighbors and background. Report the two backgrounds "
     "and effects without calling this a change caused by annotation.",
     {"roi_nodes": 2, "context_nodes": 5, "roi_background": 1, "context_background": 0.5,
      "roi_effect": 0, "context_effect": 0.5, "selected_count": 2}),
    ("missing_population", "cell_typing",
     "The researcher revises both Myeloid labels to Stromal. Is a source-to-Myeloid neighbor effect "
     "still evaluable? Distinguish unavailable evidence from a negative biological finding.",
     {"status": "indeterminate", "after_classification": "inconclusive", "effect_delta": None,
      "reason": "target_label_absent_from_neighbor_universe"}),
    ("resume_and_recompute", "collaboration",
     "After a completed label correction, reopen the workspace, recover its exact current version, "
     "and independently recompute the exported before/after comparison. Retain the original source.",
     {"revision_matches": True, "source_matches": True, "old_label": "Myeloid",
      "current_label": "Stromal", "recompute_status": "matched", "integrity": "verified",
      "scientific_authorization": "NOT_ESTABLISHED"}),
    ("unavailable_units_and_partial_data", "qc",
     "One snapshot has coordinates with unknown units, and a separate snapshot has calibrated "
     "centroids but no expression panel. Determine which operations remain meaningful. "
     "Do not interpret unknown coordinate distances as micrometers or missing counts as zeros.",
     {"uncalibrated_radius_rejected": True, "uncalibrated_exploration_available": True, "partial_selection_count": 3,
      "partial_marker_status": "unmeasured", "partial_marker_mean": None}),
    ("spot_observation_semantics", "cell_typing",
     "The same synthetic five observations are declared as spatial spots. Recover the observation "
     "unit from the saved dataset. Count selected spots without presenting them as five measured cells.",
     {"observation_unit": "spot", "observation_count": 5, "selected_observation_count": 3,
      "scientific_authorization": "NOT_ESTABLISHED"}),
)


def task_catalog() -> list[dict]:
    """Public task descriptions without grader targets; inputs are synthetic."""
    return [{"task_id": task_id, "category": category, "prompt": prompt,
             "input_kind": "synthetic", "suite_version": SUITE_VERSION}
            for task_id, category, prompt, _ in _TASKS]


def _matches(observed, expected) -> bool:
    if expected is None:
        return observed is None
    if isinstance(expected, bool):
        return type(observed) is bool and observed == expected
    if isinstance(expected, (int, float)):
        return (type(observed) in (int, float) and math.isfinite(observed)
                and math.isclose(observed, expected, rel_tol=1e-10, abs_tol=1e-12))
    if isinstance(expected, list):
        return (isinstance(observed, list) and all(isinstance(item, str) for item in observed)
                and len(set(observed)) == len(observed) and set(observed) == set(expected))
    return type(observed) is type(expected) and observed == expected


def grade_answer(task_id: str, answer: dict) -> dict:
    """Grade only these public synthetic fixtures; no claim of independent truth."""
    expected = next((spec[3] for spec in _TASKS if spec[0] == task_id), None)
    if expected is None:
        raise SpatialError(f"Unknown synthetic task: {task_id}")
    if not isinstance(answer, dict):
        raise SpatialError("A task answer must be a JSON object.")
    checks = [{"field": name, "passed": name in answer and _matches(answer[name], value),
               "expected": value, "observed": answer.get(name), "present": name in answer}
              for name, value in expected.items()]
    return {"task_id": task_id, "passed": all(item["passed"] for item in checks), "checks": checks,
            "evidence_ceiling": EVIDENCE_CEILING}


def _snapshot(*, panel: bool = True, units: str = "micrometer", observation_unit: str = "cell") -> dict:
    specifications = [("source", 0, "T cell"), ("near-target", 1, "Myeloid"),
                      ("near-other", 2, "Stromal"), ("far-target", 20, "Myeloid"),
                      ("far-other", 30, "Stromal")]
    cells = [{"cell_id": identifier, "x": x, "y": 0, "label": label, "included": True,
              "counts": ({"CD3D": 7} if identifier == "source" else {}) if panel else {}}
             for identifier, x, label in specifications]
    return {"cells": cells, "metadata": {
        "name": "Synthetic collaboration task; not tissue", "slice_id": "synthetic-section",
        "coordinate_system": "synthetic_xy", "units": units,
        "panel_genes": ["CD3D", "LST1", "EPCAM"] if panel else [],
        "source_kind": "synthetic", "source_files": [], "biological_replicates": 0,
        "observation_unit": observation_unit,
        "label_semantics": "spot_annotation" if observation_unit == "spot" else "cell_type",
        "limitations": ["Developer-authored synthetic task, no biological ground truth.",
                        "No tissue images, segmentation or independent donor evidence."],
    }}


def _save(path: Path, payload: dict | list) -> None:
    # Every result lives under an exclusively created run directory.
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")


def _new_project(folder: Path, snapshot: dict, *, name: str = "project") -> Project:
    _save(folder / f"{name}-input.json", snapshot)
    return Project.create(folder / name, snapshot["cells"], snapshot["metadata"])


def _revision(project: Project) -> str:
    return project.summary()["head_revision"]


def _revise(project: Project, ids: list[str], changes: dict, trace: list[dict]) -> str:
    base = _revision(project)
    selection = project.set_selection(base, cell_ids=ids, name="Synthetic researcher correction")
    proposal = project.propose_revision(base, selection["selection_id"], changes,
                                        "Scripted synthetic correction; no human review claim", "reference_driver")
    committed = project.apply_revision(proposal["proposal_id"], base,
                                       "SYNTHETIC REFERENCE DRIVER", True)
    trace.append({"operation": "select_preview_apply", "selection_id": selection["selection_id"],
                  "proposal_id": proposal["proposal_id"], "base_revision": base,
                  "target_revision": committed["revision_id"], "cell_ids": ids, "changes": changes,
                  "actor_kind": "scripted_synthetic_driver"})
    return committed["revision_id"]


def _compare(project: Project, before: str, after: str, trace: list[dict], **kwargs) -> dict:
    result = analysis.compare(project, before, after, radius_um=1.1, min_effect=0.2, **kwargs)
    trace.append({"operation": "compare", "run_id": result["run_id"],
                  "base_revision": before, "target_revision": after, "parameters": result["parameters"],
                  "source_sha256": result["data_hashes"]["source_sha256"]})
    return result


def _execute(task_id: str, folder: Path, trace: list[dict] | None = None) -> tuple[dict, list[dict]]:
    trace = [] if trace is None else trace
    observation_unit = "spot" if task_id == "spot_observation_semantics" else "cell"
    project = _new_project(folder, _snapshot(observation_unit=observation_unit))
    base = _revision(project)
    trace.append({"operation": "import", "project_id": project.summary()["project_id"],
                  "revision_id": base, "source_sha256": project.summary()["source_sha256"]})

    if task_id == "panel_evidence":
        selected = project.set_selection(base, cell_ids=["source"])
        inspected = project.inspect_selection(["LST1", "KDR"])
        trace.append({"operation": "inspect_selection", "selection_id": selected["selection_id"],
                      "genes": ["LST1", "KDR"]})
        return {"measured_zero_status": inspected["expression"]["LST1"]["status"],
                "measured_zero_mean": inspected["expression"]["LST1"]["mean"],
                "outside_panel_status": inspected["expression"]["KDR"]["status"],
                "outside_panel_mean": inspected["expression"]["KDR"]["mean"],
                "selected_count": inspected["cell_count"]}, trace

    if task_id in {"spatial_object_identity", "spot_observation_semantics"}:
        selected = project.set_selection(base, polygon=[[-1, -1], [2, -1], [2, 1], [-1, 1]])
        trace.append({"operation": "polygon_selection", "selection": selected})
        if task_id == "spot_observation_semantics":
            summary = project.summary()
            return {"observation_unit": summary["metadata"].get("observation_unit"),
                    "observation_count": summary["cell_count"],
                    "selected_observation_count": selected["cell_count"],
                    "scientific_authorization": summary["scientific_authorization"]}, trace
        return {"selected_ids": selected["cell_ids"], "units": selected["units"],
                "coordinate_system": selected["coordinate_system"], "revision_matches": selected["revision_id"] == base,
                "source_matches": selected["source_sha256"] == project.summary()["source_sha256"]}, trace

    if task_id == "annotation_sensitivity":
        old = _compare(project, base, base, trace)
        target = _revise(project, ["near-target"], {"label": "Stromal"}, trace)
        compared = _compare(project, base, target, trace)
        region = _revise(project, ["near-other"], {"region": "candidate niche"}, trace)
        unchanged = _compare(project, target, region, trace)
        return {"before_effect": compared["before"]["excess_over_abundance"],
                "after_effect": compared["after"]["excess_over_abundance"],
                "status": compared["comparison"]["status"],
                "effect_delta": compared["comparison"]["effect_delta"],
                "label_changes": compared["comparison"]["label_changes"],
                "region_only_status": unchanged["comparison"]["status"],
                "region_only_effect_delta": unchanged["comparison"]["effect_delta"],
                "old_result_historical": project.get_run(old["run_id"])["stale"]}, trace

    if task_id == "neighbor_context":
        selected = project.set_selection(base, cell_ids=["source", "near-target"])
        roi = _compare(project, base, base, trace, selection_id=selected["selection_id"], graph_scope="roi_induced")
        context = _compare(project, base, base, trace, selection_id=selected["selection_id"], graph_scope="whole_slice")
        return {"roi_nodes": roi["before"]["graph"]["node_count"],
                "context_nodes": context["before"]["graph"]["node_count"],
                "roi_background": roi["before"]["null_fraction"],
                "context_background": context["before"]["null_fraction"],
                "roi_effect": roi["before"]["excess_over_abundance"],
                "context_effect": context["before"]["excess_over_abundance"],
                "selected_count": roi["comparison"]["selected_cell_count"]}, trace

    if task_id == "missing_population":
        target = _revise(project, ["near-target", "far-target"], {"label": "Stromal"}, trace)
        compared = _compare(project, base, target, trace)
        return {"status": compared["comparison"]["status"],
                "after_classification": compared["after"]["classification"],
                "effect_delta": compared["comparison"]["effect_delta"],
                "reason": compared["after"]["reason"]}, trace

    if task_id == "resume_and_recompute":
        source_hash = project.summary()["source_sha256"]
        target = _revise(project, ["near-target"], {"label": "Stromal"}, trace)
        result = _compare(project, base, target, trace)
        reopened = Project(project.root)
        exported = reopened.export_bundle(result["run_id"])
        verified = verify_bundle(exported["export_path"])
        trace.append({"operation": "reopen_export_recompute", "export_path": exported["export_path"],
                      "integrity": verified["integrity"], "recompute_status": verified["recompute_status"]})
        before = {cell["cell_id"]: cell for cell in reopened.cells(base)}
        after = {cell["cell_id"]: cell for cell in reopened.cells(target)}
        return {"revision_matches": _revision(reopened) == target,
                "source_matches": reopened.summary()["source_sha256"] == source_hash,
                "old_label": before["near-target"]["label"], "current_label": after["near-target"]["label"],
                "recompute_status": verified["recompute_status"], "integrity": verified["integrity"],
                "scientific_authorization": verified["scientific_authorization"]}, trace

    if task_id == "unavailable_units_and_partial_data":
        rejected = False
        uncalibrated = _new_project(folder, _snapshot(units="unknown"), name="uncalibrated")
        native = uncalibrated.set_selection(_revision(uncalibrated), cell_ids=["source"])
        try:
            analysis.compare(uncalibrated, _revision(uncalibrated), _revision(uncalibrated), native["selection_id"])
        except SpatialError as exc:
            rejected = True
            trace.append({"operation": "uncalibrated_physical_radius", "error": str(exc)})
        partial = _new_project(folder, _snapshot(panel=False), name="partial")
        selected = partial.set_selection(_revision(partial), cell_ids=["source", "near-target", "near-other"])
        inspected = partial.inspect_selection(["CD3D"])
        trace.append({"operation": "partial_data_selection", "selection_id": selected["selection_id"],
                      "expression": inspected["expression"]})
        return {"uncalibrated_radius_rejected": rejected, "uncalibrated_exploration_available": native["units"] == "unknown", "partial_selection_count": selected["cell_count"],
                "partial_marker_status": inspected["expression"]["CD3D"]["status"],
                "partial_marker_mean": inspected["expression"]["CD3D"]["mean"]}, trace
    raise SpatialError(f"No reference driver for {task_id}")


def run_benchmark(output_dir: str | Path) -> dict:
    """Run eight local scripted tasks in a new child folder; retain failures.

    This does not grade a language model, test real-host collaboration or estimate
    researcher time savings. It never opens or revises an existing user project.
    """
    root = Path(output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    run_dir = root / f"synthetic-collaboration-{uuid.uuid4().hex[:12]}"
    run_dir.mkdir(exist_ok=False)
    tasks = task_catalog()
    _save(run_dir / "tasks.json", tasks)
    results = []
    for task in tasks:
        task_id = task["task_id"]
        folder = run_dir / task_id
        folder.mkdir()
        _save(folder / "task.json", task)
        started = time.monotonic()
        trace = []
        try:
            answer, trace = _execute(task_id, folder, trace)
            graded = grade_answer(task_id, answer)
            result = {**graded, "answer": answer, "trace": trace, "error": None}
            json.dumps(result, allow_nan=False)
        except Exception as exc:
            # This is a test harness: retain a failed task and continue other tasks.
            result = {"task_id": task_id, "passed": False, "checks": [], "answer": None,
                      "trace": trace, "error": {"type": type(exc).__name__, "message": str(exc)},
                      "evidence_ceiling": EVIDENCE_CEILING}
        result["elapsed_seconds"] = round(time.monotonic() - started, 6)
        _save(folder / "result.json", result)
        results.append({**result, "result_path": str(folder / "result.json")})
    passed = sum(result["passed"] for result in results)
    report = {"suite_version": SUITE_VERSION, "generated_at": datetime.now(timezone.utc).isoformat(),
              "execution_mode": "scripted_reference_driver", "data_kind": "synthetic",
              "evidence_ceiling": EVIDENCE_CEILING, "scientific_authorization": "NOT_ESTABLISHED",
              "external_validation": "NOT_ESTABLISHED", "model_evaluation": "NOT_PERFORMED",
              "human_efficiency_evaluation": "NOT_PERFORMED", "spatialbench_evaluation": "NOT_PERFORMED",
              "task_count": len(results), "passed_count": passed, "failed_count": len(results) - passed,
              "all_passed": passed == len(results), "output_dir": str(run_dir), "results": results,
              "limitations": ["Public developer-authored inputs and expected answers; not held-out validation.",
                              "Reference operations are scripted; no autonomous agent or real researcher was evaluated.",
                              "Recorded latency is local software runtime, not researcher or agent task time."]}
    _save(run_dir / "report.json", report)
    # Bytes of task inputs, observed outputs and reports can be checked independently.
    files = {}
    for path in sorted(run_dir.rglob("*.json")):
        content = path.read_bytes()
        files[path.relative_to(run_dir).as_posix()] = {"sha256": hashlib.sha256(content).hexdigest(),
                                                      "bytes": len(content)}
    _save(run_dir / "manifest.json", {"algorithm": "sha256", "files": files,
                                      "meaning": "Internal byte integrity; not scientific validity."})
    return deepcopy(report)
