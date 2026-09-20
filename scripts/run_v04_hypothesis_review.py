"""Rehearse frozen v2 plans on an unchanged, previously audited real workspace.

Run from the repository with PYTHONPATH=src. Existing destinations are rejected.
This performs no label fitting, downloads, or human approval attribution.
"""
from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sqlite3
from time import perf_counter

import numpy as np

from spatial_collab import __version__
from spatial_collab.hypotheses import (
    create_hypothesis_plan, freeze_hypothesis_plan, hypothesis_context,
    revise_hypothesis_plan, run_hypothesis_plan,
)
from spatial_collab.replay import verify_bundle
from spatial_collab.store import Project, _hash


def read_source(folder):
    with sqlite3.connect((folder / "project.sqlite3").as_uri() + "?mode=ro", uri=True) as db:
        payload, checksum = db.execute("SELECT payload,sha256 FROM source WHERE id=1").fetchone()
        head = db.execute("SELECT value FROM state WHERE key='head_revision'").fetchone()[0]
    source = json.loads(payload)
    if _hash(source) != checksum:
        raise ValueError("Predecessor source integrity failed.")
    return source, checksum, head


def file_hash(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write(folder, name, value):
    path = folder / name
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    return path


def brute_force(cells, selected, radii, source="T_program", target="B_program"):
    """Independent direct pairwise distances; no plugin graph/evaluator used."""
    active = [c for c in cells if c["included"]]
    answers = {}
    for scope in ("whole_slice", "roi_induced"):
        nodes = active if scope == "whole_slice" else [c for c in active if c["cell_id"] in selected]
        xy = np.array([[c["x"], c["y"]] for c in nodes])
        sources = [i for i, c in enumerate(nodes) if c["cell_id"] in selected and c["label"] == source]
        targets = np.array([c["label"] == target for c in nodes])
        totals = {radius: {"edges": 0, "target_edges": 0, "isolated": 0, "ratios": []} for radius in radii}
        for i in sources:
            difference = xy - xy[i]
            distances_squared = difference[:, 0] ** 2 + difference[:, 1] ** 2
            notself = np.arange(len(nodes)) != i
            for radius, values in totals.items():
                neighbor = (distances_squared <= radius ** 2) & notself
                degree = int(neighbor.sum())
                hits = int((neighbor & targets).sum())
                values["edges"] += degree
                values["target_edges"] += hits
                values["isolated"] += int(degree == 0)
                if degree:
                    values["ratios"].append(hits / degree)
        for radius, values in totals.items():
            null = int(targets.sum()) / (len(nodes) - 1)
            edge = values["target_edges"] / values["edges"]
            equal = sum(values["ratios"]) / len(values["ratios"])
            answers[scope, radius] = {
                "source_count": len(sources), "node_count": len(nodes), "target_count": int(targets.sum()),
                "isolated_sources": values["isolated"], "sources_with_neighbors": len(values["ratios"]),
                "directed_source_edges": values["edges"], "target_edges": values["target_edges"],
                "methods": {"edge_weighted": {"observed_fraction": edge, "null_fraction": null,
                                               "excess_over_abundance": edge - null},
                            "source_equal_weighted": {"observed_fraction": equal, "null_fraction": null,
                                                       "excess_over_abundance": equal - null}},
            }
    return answers


def main():
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-study", type=Path,
                        default=repo / "projects/real-data-review/study-v03-complete")
    parser.add_argument("--destination", type=Path, default=repo / "projects/v04-real-data/lymph-node-study")
    args = parser.parse_args()
    source_study = args.source_study.resolve()
    destination = args.destination.resolve()
    if destination.exists():
        raise ValueError("Destination already exists; choose a fresh project to preserve previous runs.")
    if destination == source_study or source_study in destination.parents:
        raise ValueError("Create this rehearsal outside the predecessor study.")
    original, old_hash, old_head = read_source(source_study / "working-project")
    inherited = json.loads((source_study / "declared-variants.json").read_text())
    previous_plan = json.loads((source_study / "frozen-plan.json").read_text())
    previous_run = json.loads((source_study / "sensitivity-whole_slice.json").read_text())
    selected = set(previous_run["selection"]["cell_ids"])
    if "CD3D" in original["metadata"]["panel_genes"]:
        raise ValueError("The declared unmeasured-feature control no longer applies to this input.")
    source_files = [{**entry, "sha256_current": file_hash(Path(entry["path"]))}
                    for entry in original["metadata"]["source_files"]]
    if any(entry["sha256"] != entry["sha256_current"] for entry in source_files):
        raise ValueError("An original raw input changed since the audited study.")
    metadata = deepcopy(original["metadata"])
    metadata["name"] = "v0.4 real lymph-node frozen hypothesis review"
    metadata["continued_from"] = {"project_path": str(source_study / "working-project"),
                                  "source_sha256": old_hash, "revision_id": old_head,
                                  "label_evidence": "Inherited deterministic expression working groups; not validated cell identities."}
    started = perf_counter()
    project = Project.create(destination, deepcopy(original["cells"]), metadata)
    head = project.summary()["head_revision"]
    selection = project.set_selection(head, cell_ids=sorted(selected))
    variants = [{"name": v["name"], "rationale": v["rationale"],
                 "selector": {"cell_ids": v["cell_ids"]}, "changes": v["changes"]} for v in inherited]
    variants[2]["rationale"] += " Actual inherited ID set includes only B_program/T_program groups; Mixed_program is unchanged."
    spec = {"name": "Real lymph-node annotation and weighting sensitivity",
            "rationale": "Predeclare inherited local alternatives and compare edge versus source weighting; retain an unmeasured-feature unknown control.",
            "revision_id": head, "selection_id": selection["selection_id"],
            "source_label": "T_program", "target_label": "B_program", "radii_um": [15, 35, 75],
            "graph_scopes": ["whole_slice", "roi_induced"], "min_effect": 0.1,
            "background": {"kind": "neighbor_universe"}, "variants": variants}
    draft = create_hypothesis_plan(project, spec)
    contexts = [{"stage": "initial_draft", **hypothesis_context(project)}]
    revised_spec = deepcopy(spec)
    revised_spec["variants"].append({"name": "unmeasured_CD3D_negative_control",
                                     "rationale": "CD3D is absent from this measured panel; membership must remain unknown rather than imputed zero or dropped.",
                                     "selector": {"predicate": {"all": [
                                         {"field": "counts.CD3D", "op": "ge", "value": 1}]}},
                                     "changes": {"label": "Uncertain"}})
    revised = revise_hypothesis_plan(project, draft["plan_id"], draft["version"], revised_spec)
    contexts.append({"stage": "revised_draft_before_execution", **hypothesis_context(project)})
    frozen = freeze_hypothesis_plan(project, revised["plan_id"], revised["version"])
    contexts.append({"stage": "frozen_before_execution", **hypothesis_context(project)})
    write(destination, "declared-spec.json", revised_spec)
    write(destination, "frozen-plan.json", frozen)
    write(destination, "plan-context-transitions.json", contexts)
    result = run_hypothesis_plan(project, frozen["plan_id"], frozen["version"])
    write(destination, "hypothesis-result.json", result)
    write(destination, "complete-row-example.json", result["results"][0])

    # Independent numerical implementation uses only copied input values and the
    # exact declared alternatives, not intermediate plugin statistics.
    expected = {"baseline": brute_force(original["cells"], selected, spec["radii_um"])}
    for variant in inherited:
        ids = set(variant["cell_ids"])
        changed = [dict(c, **variant["changes"]) if c["cell_id"] in ids else c for c in original["cells"]]
        expected[variant["name"]] = brute_force(changed, selected, spec["radii_um"])
    numeric_mismatches = []
    legacy_mismatches = []
    comparisons = []
    for row in result["results"]:
        if row["status"] != "computed":
            continue
        for phase, state in (("before", "baseline"), ("after", row["variant"])):
            actual = row[phase]
            wanted = expected[state][row["graph_scope"], row["radius_um"]]
            for key, value in wanted.items():
                if key != "methods" and actual["graph"][key] != value:
                    numeric_mismatches.append({"row": row["row_id"], "phase": phase, "field": key})
            for method in wanted["methods"]:
                for key, value in wanted["methods"][method].items():
                    if not np.isclose(actual["methods"][method][key], value, atol=1e-12, rtol=0):
                        numeric_mismatches.append({"row": row["row_id"], "phase": phase, "method": method, "field": key})
            comparisons.append({"row_id": row["row_id"], "phase": phase, "expected": wanted})
        legacy = json.loads((source_study / f"sensitivity-{row['graph_scope']}.json").read_text())
        previous = next(r for r in legacy["results"] if r["variant"] == row["variant"] and r["radius_um"] == row["radius_um"])
        for phase in ("before", "after"):
            if not np.isclose(previous[phase]["excess_over_abundance"],
                              row[phase]["methods"]["edge_weighted"]["excess_over_abundance"], atol=1e-12, rtol=0):
                legacy_mismatches.append({"row_id": row["row_id"], "phase": phase})
    numeric = {"implementation": "Independent direct source-to-candidate Euclidean squared-distance loops and scalar equal-source averaging; no plugin evaluator or spatial index.",
               "status": "PASS" if not numeric_mismatches and not legacy_mismatches else "FAIL",
               "comparisons": comparisons, "numeric_mismatches": numeric_mismatches,
               "legacy_edge_weighted_mismatches": legacy_mismatches,
               "external_validation": "NOT_ESTABLISHED"}
    write(destination, "independent-numeric-receipt.json", numeric)
    unknown = [r for r in result["results"] if r["variant"] == "unmeasured_CD3D_negative_control"]
    unknown_valid = len(unknown) == 6 and all(r["status"] == "unknown" and r["selector_unknown_count"] == len(original["cells"])
                                           and r["after"]["reason"] == "hypothesis_selector_has_unknown_observations" for r in unknown)
    exported = project.export_bundle(result["run_id"])
    try:
        replay = verify_bundle(exported["export_path"])
    except Exception as exc:
        replay = {"recompute_matches": False, "error_type": type(exc).__name__, "error": str(exc)}
    write(destination, "bundle-replay.json", replay)
    _, old_hash_after, old_head_after = read_source(source_study / "working-project")
    summary = project.summary()
    checks = {"24_declared_grid_rows_retained": len(result["results"]) == 24,
              "all_rows_in_difference_order": set(result["difference_order"]) == {r["row_id"] for r in result["results"]},
              "unknown_negative_control_retained": unknown_valid,
              "two_estimands_on_every_row": all(set(r["after"]["methods"]) == {"edge_weighted", "source_equal_weighted"} for r in result["results"]),
              "all_independent_numeric_checks_match": numeric["status"] == "PASS",
              "new_scientific_head_unchanged": summary["head_revision"] == head and summary["revision_count"] == 1,
              "old_source_and_head_unchanged": old_hash == old_hash_after and old_head == old_head_after,
              "all_context_tokens_changed": len({c["state_sha256"] for c in contexts}) == 3,
              "offline_bundle_replay_matches": replay.get("recompute_matches") is True}
    differences = [r["reference_comparison"][phase]["source_equal_minus_edge_weighted"]
                   for r in result["results"] if r["status"] == "computed" for phase in ("before", "after")]
    threshold_changes = [{"row_id": r["row_id"], "method": method} for r in result["results"]
                         for method, comparison in r["comparison"].items() if comparison["threshold_crossed"]]
    receipt = {"status": "PASS" if all(checks.values()) else "REVIEW_NEEDED", "package_version": __version__,
               "script_sha256": file_hash(Path(__file__)), "source_study": str(source_study),
               "predecessor_source_sha256": old_hash, "raw_file_checks": source_files,
               "working_label_counts": dict(Counter(c["label"] for c in original["cells"])),
               "window_observations": len(original["cells"]), "roi_observations": len(selected),
               "plan_id": frozen["plan_id"], "plan_version": frozen["version"], "plan_sha256": frozen["plan_sha256"],
               "run_id": result["run_id"], "run_sha256": result["run_sha256"], "status_counts": result["status_counts"],
               "checks": checks, "reference_minus_primary_excess_range": [min(differences), max(differences)],
               "threshold_changes": threshold_changes, "export_path": exported["export_path"],
               "elapsed_seconds": perf_counter() - started,
               "scientific_authorization": "NOT_ESTABLISHED", "expert_review": "NOT_PERFORMED",
               "interpretation": "Local descriptive comparison of explicit annotation alternatives and weighting estimands. Unknown CD3D membership is not a negative biological result. Threshold stability is limited to this grid and cannot validate cell identities or absence of biological relationships.",
               "inherited_plan": {"study_id": previous_plan["study_id"], "role": "Historical explanatory context, not an independent preregistration."}}
    write(destination, "acceptance-receipt.json", receipt)
    lines = ["# Real lymph-node frozen-plan review", "", f"Status: {receipt['status']}", "",
             f"The new workspace contains {len(original['cells'])} observations and a fixed {len(selected)}-observation ROI. Original expression-program labels are inherited hypotheses, not confirmed identities.", "",
             f"All 24 grid rows are retained: {result['status_counts']}. The CD3D control remains unknown in all six radius/scope cells because the feature was not measured.", "",
             f"Independent direct-distance checks: {numeric['status']}; old edge-weighted results match: {not legacy_mismatches}. Offline replay: {replay.get('recompute_matches')}. Both old and new scientific heads are preserved.", "",
             f"Equal-source minus edge-weighted excess ranges from {min(differences):.8f} to {max(differences):.8f} over computed before/after states. Threshold crossings: {len(threshold_changes)}. Differences describe weighting sensitivity, not which method is biologically correct.", "",
             "| Scope | Variant | Radius | Primary delta | Equal-source delta | Status |", "|---|---|---:|---:|---:|---|"]
    for row in result["results"]:
        delta = [row["comparison"][m]["effect_delta"] for m in ("edge_weighted", "source_equal_weighted")]
        display = ["unknown" if v is None else f"{v:.8f}" for v in delta]
        lines.append(f"| {row['graph_scope']} | {row['variant']} | {row['radius_um']:g} | {display[0]} | {display[1]} | {row['status']} |")
    (destination / "REVIEW.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({"destination": str(destination), "status": receipt["status"], "checks": checks,
                      "status_counts": result["status_counts"], "reference_difference_range": receipt["reference_minus_primary_excess_range"]}))


if __name__ == "__main__":
    main()
