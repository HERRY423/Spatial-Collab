import pytest

from spatial_collab.analysis import compare
from spatial_collab.exploration import compare_regions, query_observations
from spatial_collab.replay import verify_bundle
from spatial_collab.sensitivity import run_sensitivity
from spatial_collab.server import ToolService
from spatial_collab.store import Project, SpatialError


@pytest.mark.parametrize("units", ["unknown", "array_index", "pixel"])
def test_uncalibrated_expression_selection_revision_and_export_remain_available(tmp_path, units):
    project = Project.create(tmp_path / "p", [
        {"cell_id": "a", "x": 10, "y": 15, "label": "A", "counts": {"G": 5}},
        {"cell_id": "b", "x": 13, "y": 17, "label": "B", "counts": {"G": 2}}],
        {"name": "uncalibrated test", "slice_id": "s", "coordinate_system": "recorded_array",
         "units": units, "source_kind": "synthetic", "panel_genes": ["G"]})
    rev = project.summary()["head_revision"]
    selected = project.set_selection(rev, polygon=[[9, 14], [11, 14], [11, 16], [9, 16]])
    assert selected["cell_ids"] == ["a"] and selected["units"] == units
    assert project.inspect_selection(["G"])["expression"]["G"]["mean"] == 5
    assert query_observations(project, rev, gene="G", min_count=4)["matching_count"] == 1
    run = compare_regions(project, rev, rev, selected["selection_id"], genes=["G"])
    assert run["before"]["marker_contrast"]["G"]["mean_raw_count_difference"] == 3
    assert verify_bundle(project.export_bundle(run["run_id"])["export_path"])["recompute_matches"]
    capabilities = ToolService(project.root).open_project()["coordinate_capabilities"]
    assert capabilities["expression_review"] and not capabilities["physical_radius_analysis"]
    with pytest.raises(SpatialError, match="micrometer"):
        compare(project, rev, rev, selected["selection_id"])
    with pytest.raises(SpatialError, match="micrometer"):
        run_sensitivity(project, rev, selected["selection_id"],
                        [{"name": "hyp", "cell_ids": ["b"], "changes": {"label": "C"}, "rationale": "test"}], [35.0], "A", "B")
    proposal = project.propose_revision(rev, selected["selection_id"], {"label": "Uncertain"}, "synthetic test")
    assert proposal["base_revision"] == rev and project.summary()["head_revision"] == rev
