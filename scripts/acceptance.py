"""Reproducible synthetic threshold-flip + export/replay acceptance.

Usage: python scripts/acceptance.py <new-output-directory>
"""
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from spatial_collab.analysis import compare
from spatial_collab.replay import verify_bundle
from spatial_collab.store import Project, SpatialError


def acceptance(destination):
    project = Project.create(destination, [
        {"cell_id": key, "x": x, "y": 0, "label": label, "counts": {"CD3D": count}}
        for key, x, label, count in [("source", 0, "T cell", 8), ("near-target", 1, "Myeloid", 0),
            ("near-other", 2, "Stromal", 0), ("far-target", 20, "Myeloid", 0),
            ("far-other", 30, "Stromal", 0)]
    ], {"name": "Synthetic hand-computable acceptance", "slice_id": "synthetic-one",
        "coordinate_system": "synthetic_xy", "units": "micrometer", "panel_genes": ["CD3D"],
        "source_kind": "synthetic", "source_files": [], "biological_replicates": 0,
        "limitations": ["Synthetic software acceptance, not measured tissue or scientific confirmation."]})
    base = project.summary()["head_revision"]
    baseline = compare(project, base, base, radius_um=1.1, min_effect=.2)
    selection = project.set_selection(base, cell_ids=["near-target"])
    evidence = project.inspect_selection(["CD3D", "UNMEASURED"])
    assert evidence["expression"]["CD3D"]["mean"] == 0
    assert evidence["expression"]["UNMEASURED"]["measured"] is False
    proposal = project.propose_revision(base, selection["selection_id"], {"label": "Stromal"},
                                         "Synthetic correction demonstrates descriptive sensitivity.", "test")
    assert project.summary()["head_revision"] == base
    target = project.apply_revision(proposal["proposal_id"], base, "Synthetic acceptance test", True)["revision_id"]
    assert project.get_run(baseline["run_id"])["stale"] is True
    try:
        project.apply_revision(proposal["proposal_id"], base, "Synthetic acceptance test", True)
        raise AssertionError("Stale commit was accepted")
    except SpatialError:
        pass
    result = compare(project, base, target, radius_um=1.1, min_effect=.2)
    assert result["before"]["excess_over_abundance"] == .5
    assert result["after"]["excess_over_abundance"] == -.25
    assert result["comparison"]["status"] == "changed"
    exported = project.export_bundle(result["run_id"])
    replay = verify_bundle(exported["export_path"])
    assert replay["recompute_matches"] is True
    reverted = project.revert_revision(base, target, "Synthetic acceptance test", True)["revision_id"]
    assert project.cells(reverted) == project.cells(base)
    report = {"status": "passed", "data_kind": "synthetic", "project": str(project.root),
              "source_sha256": project.summary()["source_sha256"], "before_excess": .5, "after_excess": -.25,
              "effect_delta": result["comparison"]["effect_delta"], "comparison_status": "changed",
              "stale_commit_rejected": True, "baseline_result_marked_stale": True,
              "original_revision_preserved": True, "revert_creates_new_revision": reverted not in (base, target),
              "export_path": exported["export_path"], "replay": replay,
              "scientific_authorization": "NOT_ESTABLISHED"}
    (project.root / "acceptance.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python scripts/acceptance.py <new-output-directory>")
    print(json.dumps(acceptance(sys.argv[1]), ensure_ascii=False, indent=2))

