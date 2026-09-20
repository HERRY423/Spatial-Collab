import math

import pytest

from spatial_collab.exploration import compare_regions, query_observations
from spatial_collab.replay import verify_bundle
from spatial_collab.server import ToolService
from spatial_collab.store import Project, SpatialError


@pytest.fixture
def project(tmp_path):
    cells = [{"cell_id": str(i), "x": float(i), "y": 0.0, "label": label,
              "counts": counts} for i, (label, counts) in enumerate([
                  ("A", {"G": 4, "H": 4}), ("B", {}), ("A", {"H": 8}), ("B", {"G": 2})])]
    return Project.create(tmp_path / "project", cells, {"name": "Fixture", "slice_id": "s1", "source_kind": "synthetic",
                         "coordinate_system": "slide", "units": "micrometer", "panel_genes": ["G", "H"]})


def edit(project, rev, ids, changes):
    sel = project.set_selection(rev, cell_ids=ids)
    prop = project.propose_revision(rev, sel["selection_id"], changes, "test", "test driver")
    new = project.apply_revision(prop["proposal_id"], rev, "Test fixture", True)
    return new["revision_id"], prop


def test_paging_is_pinned_across_edits_and_rejects_query_changes(project):
    rev = project.summary()["head_revision"]
    first = query_observations(project, rev, limit=1, labels=["A"])
    new, _ = edit(project, rev, ["2"], {"label": "B"})
    second = query_observations(project, rev, limit=1, labels=["A"], cursor=first["next_cursor"])
    assert [c["cell_id"] for c in first["observations"] + second["observations"]] == ["0", "2"]
    assert second["historical"] and second["complete"]
    assert query_observations(project, new, labels=["A"])["matching_count"] == 1
    for kwargs in ({"revision_id": new, "labels": ["A"]}, {"revision_id": rev, "labels": ["B"]}):
        with pytest.raises(SpatialError, match="Cursor"):
            query_observations(project, **kwargs, cursor=first["next_cursor"])


def test_unmeasured_gene_is_not_zero_filter(project):
    rev = project.summary()["head_revision"]
    assert query_observations(project, rev, gene="G", min_count=1)["matching_count"] == 2
    with pytest.raises(SpatialError, match="not measured"):
        query_observations(project, rev, gene="Missing")
    for kwargs in ({"limit": True}, {"min_count": math.nan}, {"bounds": [0, 0, 0, 2]}, {"cursor": "bad"}):
        with pytest.raises(SpatialError):
            query_observations(project, rev, **kwargs)


def test_context_and_cross_client_proposal_handoff(project):
    service = ToolService(project.root)
    original = service.call("get_context")
    assert not service.call("get_context", {"after_cursor": original["cursor"]})["changed"]
    rev = original["head_revision"]
    sel = project.set_selection(rev, cell_ids=["0"])
    prop = project.propose_revision(rev, sel["selection_id"], {"label": "C"}, "marker uncertainty")
    context = service.call("get_context", {"after_cursor": original["cursor"]})
    assert context["changed"] and context["latest_proposal_id"] == prop["proposal_id"]
    receipt = service.call("get_proposal", {"proposal_id": prop["proposal_id"]})
    assert receipt["deltas"] == prop["deltas"] and not receipt["stale"] and not receipt["applied_revision"]
    project.apply_revision(prop["proposal_id"], rev, "Human attribution", True)
    assert service.call("get_proposal", {"proposal_id": prop["proposal_id"]})["applied_revision"]
    with pytest.raises(SpatialError):
        service.call("query_observations", {"revision_id": rev, "limit": "100"})


def test_explicit_background_marker_denominators_and_replay(project):
    rev = project.summary()["head_revision"]
    foreground = project.set_selection(rev, cell_ids=["0", "1"], name="foreground")
    background = project.set_selection(rev, cell_ids=["2", "3"], name="matched background")
    result = compare_regions(project, rev, rev, foreground["selection_id"],
                             background_selection_id=background["selection_id"], genes=["G", "H", "Missing"])
    fg = result["before"]["foreground"]
    assert fg["markers"]["G"]["mean_raw_count"] == 2
    assert fg["markers"]["G"]["mean_panel_counts_per_10000"] == 5000
    assert fg["markers"]["G"]["normalized_denominator"] == 1
    assert fg["markers"]["G"]["zero_library_observations"] == 1
    assert fg["markers"]["G"]["detection_fraction"] == .5
    assert fg["markers"]["Missing"]["status"] == "unmeasured"
    assert result["before"]["marker_contrast"]["Missing"]["mean_raw_count_difference"] is None
    bundle = project.export_bundle(result["run_id"])
    assert verify_bundle(bundle["export_path"])["recompute_matches"]


def test_relabel_preserves_expression_but_exclusion_changes_denominators(project):
    rev = project.summary()["head_revision"]
    sel = project.set_selection(rev, cell_ids=["0", "1"])
    new, _ = edit(project, rev, ["0"], {"label": "B"})
    result = compare_regions(project, rev, new, sel["selection_id"], genes=["G"])
    assert result["comparison"]["status"] == "descriptive_values_changed"
    assert result["comparison"]["label_contrast_delta"] == {"A": -.5, "B": .5}
    assert result["before"]["marker_contrast"] == result["after"]["marker_contrast"]
    fresh = project.set_selection(new, cell_ids=["0", "1"])
    excluded, _ = edit(project, new, ["0"], {"included": False})
    result = compare_regions(project, new, excluded, fresh["selection_id"], genes=["G"])
    assert result["after"]["foreground"]["included_count"] == 1
    assert result["after"]["foreground"]["markers"]["G"]["normalized_denominator"] == 0
    assert result["after"]["marker_contrast"]["G"]["mean_panel_counts_per_10000_difference"] is None


def test_background_is_explicit_and_disjoint(project):
    rev = project.summary()["head_revision"]
    a = project.set_selection(rev, cell_ids=["0", "1"])
    b = project.set_selection(rev, cell_ids=["1", "2"])
    with pytest.raises(SpatialError, match="disjoint"):
        compare_regions(project, rev, rev, a["selection_id"], background_selection_id=b["selection_id"])
    all_ids = project.set_selection(rev, cell_ids=["0", "1", "2", "3"])
    with pytest.raises(SpatialError, match="empty"):
        compare_regions(project, rev, rev, all_ids["selection_id"])


def test_region_all_excluded_remains_indeterminate(project):
    rev = project.summary()["head_revision"]
    a = project.set_selection(rev, cell_ids=["0", "1"])
    new, _ = edit(project, rev, ["0", "1"], {"included": False})
    result = compare_regions(project, rev, new, a["selection_id"], genes=["G"])
    assert result["comparison"]["status"] == "indeterminate"
    assert result["after"]["foreground"]["markers"]["G"]["mean_raw_count"] is None
    assert verify_bundle(project.export_bundle(result["run_id"])["export_path"])["recompute_matches"]
