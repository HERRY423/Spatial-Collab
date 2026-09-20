from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import sqlite3

import pytest

from spatial_collab.store import Project, SpatialError


@pytest.fixture
def cells():
    return [
        {"cell_id": "c", "x": 2, "y": 0, "label": "B", "counts": {"CD3D": 0}},
        {"cell_id": "a", "x": 0, "y": 0, "label": "T", "counts": {"CD3D": 3}},
        {"cell_id": "b", "x": 1, "y": 1, "label": "T", "counts": {}},
    ]


@pytest.fixture
def metadata():
    return {"name": "synthetic fixture", "slice_id": "slice1", "coordinate_system": "slice1-micrometer",
            "units": "micrometer", "panel_genes": ["CD3D", "ZERO"], "source_kind": "synthetic",
            "source_files": [], "biological_replicates": 1, "limitations": ["Synthetic fixture only"]}


@pytest.fixture
def project(tmp_path, cells, metadata):
    return Project.create(tmp_path / "project", cells, metadata)


def prepare(project, ids=None, changes=None):
    head = project.summary()["head_revision"]
    selection = project.set_selection(head, cell_ids=ids or ["a"])
    preview = project.propose_revision(head, selection["selection_id"], changes or {"label": "B"}, "marker review")
    return head, selection, preview


def test_canonical_source_and_no_overwrite(project, tmp_path, cells, metadata):
    assert [c["cell_id"] for c in project.cells()] == ["a", "b", "c"]
    other = Project.create(tmp_path / "other", list(reversed(cells)), metadata)
    assert other.summary()["source_sha256"] == project.summary()["source_sha256"]
    with pytest.raises(SpatialError, match="already exists"):
        Project.create(project.root, cells, metadata)
    assert len(project.cells()) == 3


@pytest.mark.parametrize("field,value", [("x", math.nan), ("y", math.inf), ("x", True), ("included", 1), ("label", "")])
def test_invalid_cell_values(tmp_path, cells, metadata, field, value):
    cells[0][field] = value
    with pytest.raises(SpatialError):
        Project.create(tmp_path / "bad", cells, metadata)
    assert not (tmp_path / "bad" / "project.sqlite3").exists()


@pytest.mark.parametrize("counts", [{"CD3D": -1}, {"CD3D": float("nan")}, {"NOT_IN_PANEL": 1}, {"CD3D": True}])
def test_counts_require_measured_nonnegative_finite_panel(tmp_path, cells, metadata, counts):
    cells[0]["counts"] = counts
    with pytest.raises(SpatialError):
        Project.create(tmp_path / "bad", cells, metadata)


def test_duplicate_ids_and_uncalibrated_units(tmp_path, cells, metadata):
    cells[1]["cell_id"] = cells[0]["cell_id"]
    with pytest.raises(SpatialError, match="Duplicate"):
        Project.create(tmp_path / "bad", cells, metadata)
    cells[1]["cell_id"] = "new"
    metadata["units"] = "invented_scale"
    with pytest.raises(SpatialError, match="micrometer"):
        Project.create(tmp_path / "bad", cells, metadata)


def test_units_must_be_explicit_and_replicates_never_assumed(tmp_path, cells, metadata):
    metadata.pop("units")
    with pytest.raises(SpatialError, match="micrometer"):
        Project.create(tmp_path / "bad", cells, metadata)
    metadata["units"] = "micrometer"
    metadata.pop("biological_replicates")
    project = Project.create(tmp_path / "unknown", cells, metadata)
    assert project.summary()["metadata"]["biological_replicates"] == 0
    metadata["biological_replicates"] = 2
    with pytest.raises(SpatialError, match="single-slice"):
        Project.create(tmp_path / "replicate_bad", cells, metadata)


def test_exact_selection_unknown_ids_and_empty(project):
    head = project.summary()["head_revision"]
    selection = project.set_selection(head, cell_ids=["c", "a"])
    assert selection["cell_ids"] == ["a", "c"]
    assert selection["slice_id"] == "slice1"
    with pytest.raises(SpatialError, match="Unknown cell"):
        project.set_selection(head, cell_ids=["missing"])
    with pytest.raises(SpatialError, match="duplicates"):
        project.set_selection(head, cell_ids=["a", "a"])
    with pytest.raises(SpatialError, match="exactly one"):
        project.set_selection(head, cell_ids=[], polygon=[[0, 0], [1, 0], [0, 1]])
    empty = project.set_selection(head, cell_ids=[])
    assert empty["cell_count"] == 0
    assert project.inspect_selection(["CD3D"])["expression"]["CD3D"]["mean"] is None


def test_polygon_boundary_and_immutable_exact_universe(project):
    head = project.summary()["head_revision"]
    selection = project.set_selection(head, polygon=[[0, 0], [1, 0], [1, 1], [0, 1]])
    assert selection["cell_ids"] == ["a", "b"]
    assert selection["method"] == "centroid_containment_boundary_included"
    with pytest.raises(SpatialError, match="polygon"):
        project.set_selection(head, polygon=[[0, 0], [2, 2], [0, 2], [2, 0]])
    assert project.get_selection(selection["selection_id"])["cell_ids"] == ["a", "b"]


def test_polygon_small_area_survives_large_coordinate_origin(tmp_path, cells, metadata):
    for cell in cells:
        cell["x"] += 1e8
        cell["y"] += 1e8
    project = Project.create(tmp_path / "translated", cells, metadata)
    roi = [[1e8, 1e8], [1e8 + 1, 1e8], [1e8 + 1, 1e8 + 1], [1e8, 1e8 + 1]]
    assert project.set_selection(project.summary()["head_revision"], polygon=roi)["cell_ids"] == ["a", "b"]


def test_measured_zero_is_distinct_from_absent_gene(project):
    project.set_selection(project.summary()["head_revision"], cell_ids=["a", "b", "c"])
    inspection = project.inspect_selection(["CD3D", "ZERO", "UNMEASURED"])
    assert inspection["expression"]["CD3D"]["mean"] == 1
    assert inspection["expression"]["CD3D"]["zero_cells"] == 2
    assert inspection["expression"]["ZERO"]["mean"] == 0
    assert inspection["expression"]["UNMEASURED"]["mean"] is None
    assert inspection["expression"]["UNMEASURED"]["status"] == "unmeasured"
    assert inspection["label_counts"] == {"B": 1, "T": 2}


def test_proposal_does_not_mutate_and_commit_requires_attribution(project):
    head, selection, proposal = prepare(project)
    assert project.summary()["head_revision"] == head
    assert project.cells(head)[0]["label"] == "T"
    assert proposal["deltas"] == [{"cell_id": "a", "changes": {"label": {"before": "T", "after": "B"}}}]
    with pytest.raises(SpatialError, match="reviewer"):
        project.apply_revision(proposal["proposal_id"], head, "", True)
    with pytest.raises(SpatialError, match="confirmation"):
        project.apply_revision(proposal["proposal_id"], head, "Researcher A", False)
    with pytest.raises(SpatialError, match="confirmation"):
        project.apply_revision(proposal["proposal_id"], head, "Researcher A", 1)
    result = project.apply_revision(proposal["proposal_id"], head, "Researcher A", True)
    assert result["scientific_authorization"] == "NOT_ESTABLISHED"
    assert result["reviewer_identity"] == "user_supplied_attribution_not_authenticated"
    assert project.cells()[0]["label"] == "B"
    assert project.cells(head)[0]["label"] == "T"
    assert project.get_selection()["stale"]
    with pytest.raises(SpatialError, match="Stale selection"):
        project.inspect_selection()
    with pytest.raises(SpatialError, match="Stale"):
        project.propose_revision(result["revision_id"], selection["selection_id"], {"label": "Q"}, "again")


@pytest.mark.parametrize("changes", [{"x": 5}, {"counts": {}}, {"included": "false"}, {"label": " "}, {}])
def test_only_supported_overlays(project, changes):
    head = project.summary()["head_revision"]
    selection = project.set_selection(head, cell_ids=["a"])
    with pytest.raises(SpatialError):
        project.propose_revision(head, selection["selection_id"], changes, "test")


def test_noop_rejected(project):
    with pytest.raises(SpatialError, match="changes no"):
        prepare(project, changes={"label": "T"})


def test_revert_appends_and_restores_only_overlays(project):
    original = project.cells()
    head, _, proposal = prepare(project, changes={"label": "B", "included": False, "region": "ROI"})
    edit = project.apply_revision(proposal["proposal_id"], head, "A", True)
    restored = project.revert_revision(head, edit["revision_id"], "A", True)
    assert restored["revision_id"] not in (head, edit["revision_id"])
    assert restored["parent_revision"] == edit["revision_id"]
    assert project.cells() == original
    assert project.cells(edit["revision_id"])[0]["included"] is False
    assert project.summary()["revision_count"] == 3


def test_concurrent_compare_and_swap_only_one_commit(project):
    head, selection, first = prepare(project)
    second = project.propose_revision(head, selection["selection_id"], {"label": "Other"}, "independent review")
    def commit(preview):
        try:
            return Project(project.root).apply_revision(preview["proposal_id"], head, "A", True)
        except SpatialError as exc:
            return exc
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(commit, [first, second]))
    assert sum(isinstance(r, dict) for r in results) == 1
    assert sum(isinstance(r, SpatialError) for r in results) == 1
    assert project.summary()["revision_count"] == 2


def test_runs_are_immutable_historical_and_export_hashes_verify(project):
    head, selection, proposal = prepare(project)
    run = project.save_run({"base_revision": head, "target_revision": head,
                            "selection_id": selection["selection_id"], "effect": 0.2})
    assert run["stale"] is False
    edited = project.apply_revision(proposal["proposal_id"], head, "A", True)
    assert edited["invalidated_run_ids"] == [run["run_id"]]
    assert project.get_run(run["run_id"])["stale"] is True
    assert project.get_run(run["run_id"])["effect"] == 0.2
    bundle = project.export_bundle(run["run_id"])
    from pathlib import Path
    folder = Path(bundle["export_path"])
    manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
    for filename, entry in manifest["files"].items():
        assert hashlib.sha256((folder / filename).read_bytes()).hexdigest() == entry["sha256"]
    assert json.loads((folder / "result_status_at_export.json").read_text())["stale"]
    assert json.loads((folder / "run_selection.json").read_text())["cell_ids"] == ["a"]
    assert project.export_bundle(run["run_id"])["export_path"] != str(folder)


def test_late_finished_run_is_saved_as_historical(project):
    head, _, proposal = prepare(project)
    project.apply_revision(proposal["proposal_id"], head, "A", True)
    run = project.save_run({"base_revision": head, "target_revision": head})
    assert run["stale"]
    assert run["applicability"] == "historical_revision"
    assert run["run_sha256"] == project.get_run(run["run_id"])["run_sha256"]


def test_sqlite_append_only_guards(project):
    head, _, _ = prepare(project)
    project.save_run({"base_revision": head, "target_revision": head})
    with sqlite3.connect(project.db_path) as db:
        for table in ("source", "revisions", "selections", "proposals", "runs"):
            with pytest.raises(sqlite3.IntegrityError, match="Immutable"):
                db.execute(f"DELETE FROM {table}")


def test_unknown_revision_and_runs_rejected(project):
    with pytest.raises(SpatialError, match="Unknown"):
        project.cells("unknown")
    with pytest.raises(SpatialError):
        project.save_run({"base_revision": "unknown", "target_revision": "unknown"})
    with pytest.raises(SpatialError):
        project.get_run("unknown")


def test_inspection_cell_payload_is_bounded(tmp_path, metadata):
    large = [{"cell_id": str(i), "x": i, "y": 0, "label": "A"} for i in range(150)]
    project = Project.create(tmp_path / "large", large, metadata)
    project.set_selection(project.summary()["head_revision"], cell_ids=[c["cell_id"] for c in large])
    inspected = project.inspect_selection()
    assert inspected["cell_count"] == 150
    assert len(inspected["cells"]) == 100
    assert inspected["truncated"]


def test_finite_extreme_counts_do_not_export_infinity(tmp_path, cells, metadata):
    for cell in cells:
        cell["counts"] = {"CD3D": 1e308}
    project = Project.create(tmp_path / "extreme", cells, metadata)
    project.set_selection(project.summary()["head_revision"], cell_ids=[c["cell_id"] for c in cells])
    expression = project.inspect_selection(["CD3D"])["expression"]["CD3D"]
    assert math.isfinite(expression["mean"])
    assert expression["sum"] is None
    assert expression["sum_status"] == "overflow"
    json.dumps(expression, allow_nan=False)
