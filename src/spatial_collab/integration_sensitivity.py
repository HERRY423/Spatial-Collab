"""Extend existing frozen hypotheses with preprocessing/backend/refit comparisons."""
from . import objects
from .hypotheses import get_hypothesis_plan, PLAN_SCHEMA
from .integration_jobs import submit, get_job
from .store import SpatialError


def run_plan(project, plan_id, version, assay_id, alternatives):
    frozen = get_hypothesis_plan(project, plan_id, version)
    if frozen["schema"] != PLAN_SCHEMA or frozen["status"] != "frozen":
        raise SpatialError("Integration sensitivity requires an existing frozen hypothesis plan.")
    if not isinstance(alternatives, list) or not 1 <= len(alternatives) <= 4:
        raise SpatialError("Declare 1..4 complete backend/preprocessing alternatives.")
    allowed = {"backend", "protein_transform", "rna_features", "protein_features", "components", "clusters", "seed"}
    for alternative in alternatives:
        if not isinstance(alternative, dict) or set(alternative) - allowed or not {"backend", "protein_transform", "seed"} <= set(alternative):
            raise SpatialError("Every alternative requires backend, protein_transform and seed; other keys must be fitting parameters.")
    revision = frozen["spec"]["revision_id"]
    baseline = [c["cell_id"] for c in project.cells(revision) if c["included"]]
    variants = [{"name": "baseline", "changes": {}, "cell_ids": [], "unknown_cell_ids": []}, *frozen["resolved"]["variants"]]
    # Freeze every intended alternative before submitting any job.
    experiment = objects.put(project, "integrationexperiment", {"frozen_plan": frozen, "assay_id": assay_id,
        "alternatives": alternatives, "revision_created": False, "researcher_approval": "not_requested_computational_plan",
        "fit_scope": "whole imported population with explicit hypothetical exclusions"})
    rows = []
    for i, alternative in enumerate(alternatives):
        for variant in variants:
            row = {"alternative": i, "variant": variant["name"]}
            if variant["unknown_cell_ids"]:
                rows.append({**row, "status": "unknown", "reason": "selector_has_unknown_objects"})
                continue
            changes = variant["changes"]
            if changes.get("included") is True and not set(variant["cell_ids"]) <= set(baseline):
                rows.append({**row, "status": "unknown", "reason": "hypothetical_reinclusion_requires_explicit_revision"})
                continue
            ids = [cid for cid in baseline if changes.get("included") is not False or cid not in set(variant["cell_ids"])]
            recipe = {"revision_id": revision, "assay_id": assay_id, **alternative}
            if ids != baseline:
                recipe["observation_ids"] = ids
            try:
                job = submit(project, recipe)
                rows.append({**row, "status": job["status"], "job_id": job["id"], "observation_count": len(ids),
                    "semantics": "hypothetical_population_refit" if ids != baseline else "unchanged_unsupervised_inputs_reuse_fit",
                    "label_changes_applied_to_model": False})
            except (SpatialError, ValueError) as exc:
                rows.append({**row, "status": "failed", "reason": str(exc)})
    return objects.put(project, "integrationexperimentrun", {"experiment_id": experiment["object_id"], "rows": rows,
        "revision_created": False, "all_alternatives_retained": True})


def inspect_plan_run(project, run_id):
    result = objects.get(project, run_id, "integrationexperimentrun")
    for row in result["rows"]:
        if "job_id" in row:
            state = get_job(project, row["job_id"])
            row.update(status=state["status"], result_id=state["result_id"], error=state["error"])
    return result
