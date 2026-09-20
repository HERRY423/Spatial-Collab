import copy
import pytest

from spatial_collab.replay import verify_bundle
from spatial_collab.sensitivity import run_sensitivity
from spatial_collab.store import Project, SpatialError


@pytest.fixture
def project(tmp_path):
    return Project.create(tmp_path / "project", [
        {"cell_id": "s", "x": 0, "y": 0, "label": "T", "counts": {}},
        {"cell_id": "t", "x": 1, "y": 0, "label": "B", "counts": {}},
        {"cell_id": "b", "x": 20, "y": 0, "label": "B", "counts": {}},
        {"cell_id": "o", "x": 30, "y": 0, "label": "O", "counts": {}}],
        {"name": "test", "slice_id": "s", "coordinate_system": "native", "units": "micrometer", "source_kind": "synthetic",
         "panel_genes": [], "import_scope": {"kind": "spatial_window", "bounds": [-10, -10, 40, 40],
                                               "units": "micrometer", "sampling": "none", "selected_observation_count": 4,
                                               "selection_rule": "all cell centroids with xmin <= x <= xmax and ymin <= y <= ymax"}})


def variants():
    return [{"name": "relabel", "cell_ids": ["t"], "changes": {"label": "O"}, "rationale": "Hypothetical alternative"},
            {"name": "exclude", "cell_ids": ["t"], "changes": {"included": False}, "rationale": "Independent exclusion assumption"}]


def test_hypotheses_do_not_commit_or_fabricate_approval_and_replay(project):
    summary = project.summary()
    rev = summary["head_revision"]
    selected = project.set_selection(rev, cell_ids=["s"])
    cells = copy.deepcopy(project.cells())
    r = run_sensitivity(project, rev, selected["selection_id"], variants(), [2., 15.], "T", "B")
    assert project.summary()["head_revision"] == rev
    assert project.summary()["revision_count"] == 1
    assert project.cells() == cells and not r["revision_created"]
    assert r["researcher_approval"] == "NOT_REQUESTED_HYPOTHESIS_ONLY"
    first = r["results"][0]
    assert first["before"]["observed_fraction"] == 1
    assert first["after"]["observed_fraction"] == 0
    assert first["comparison"]["status"] == "changed"
    assert r["results"][2]["comparison"]["status"] == "indeterminate"
    assert first["neighbor_coverage"] == "complete_for_radius_within_declared_window"
    assert r["results"][1]["neighbor_coverage"] == "window_boundary_may_truncate_neighbors"
    assert verify_bundle(project.export_bundle(r["run_id"])["export_path"])["recompute_matches"]


@pytest.mark.parametrize("value", [[], [{"name": "x"}], [{"name": "x", "cell_ids": ["missing"], "changes": {"included": False}, "rationale": "x"}], [{"name": "x", "cell_ids": ["t"], "changes": {"x": 9}, "rationale": "x"}], [{"name": "x", "cell_ids": ["t"], "changes": {"included": "false"}, "rationale": "x"}]])
def test_invalid_hypotheses_leave_no_result(project, value):
    rev = project.summary()["head_revision"]
    sel = project.set_selection(rev, cell_ids=["s"])
    with pytest.raises(SpatialError):
        run_sensitivity(project, rev, sel["selection_id"], value, [2.], "T", "B")
    assert project.summary()["runs"] == [] and project.summary()["revision_count"] == 1
