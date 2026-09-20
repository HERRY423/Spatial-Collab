"""Hand-computable scientific/graph boundaries for annotation sensitivity."""

from __future__ import annotations

import json
import math

import pytest

from spatial_collab import analysis
from spatial_collab.store import Project, SpatialError


def metadata():
    return {
        "name": "Hand-computable synthetic graph", "slice_id": "slice-one",
        "coordinate_system": "local-centroids", "units": "micrometer",
        "panel_genes": ["CD3D", "LST1"], "source_kind": "synthetic",
        "source_files": [], "biological_replicates": 0, "limitations": [],
    }


def cell(cell_id, x, y, label):
    return {"cell_id": cell_id, "x": x, "y": y, "label": label,
            "included": True, "region": "", "counts": {"CD3D": 0}}


@pytest.fixture
def project(tmp_path):
    return Project.create(tmp_path / "project", [
        cell("source", 0, 0, "T cell"), cell("near-target", 1, 0, "Myeloid"),
        cell("near-other", 2, 0, "Stromal"), cell("far-target", 20, 0, "Myeloid"),
        cell("far-other", 30, 0, "Stromal"),
    ], metadata())


def head(project):
    return project.summary()["head_revision"]


def revise(project, ids, changes):
    base = head(project)
    selection = project.set_selection(expected_revision=base, cell_ids=ids)
    proposal = project.propose_revision(
        expected_revision=base, selection_id=selection["selection_id"],
        changes=changes, rationale="Synthetic test correction", actor="test",
    )
    return project.apply_revision(
        proposal_id=proposal["proposal_id"], expected_revision=base,
        reviewer="Synthetic test reviewer", confirmation=True,
    )["revision_id"]


def run(project, before=None, after=None, **kwargs):
    return analysis.compare(project, before or head(project), after or head(project),
                            radius_um=kwargs.pop("radius_um", 1.1), **kwargs)


def test_hand_computable_graph_and_immutable_saved_run(project):
    result = run(project)
    before = result["before"]
    assert before["observed_fraction"] == 1
    assert before["null_fraction"] == .5
    assert before["excess_over_abundance"] == .5
    assert before["graph"]["source_count"] == 1
    assert before["graph"]["directed_source_edges"] == 1
    assert before["graph"]["target_edges"] == 1
    assert before["classification"] == "descriptive_supported"
    assert result["comparison"]["status"] == "stable"
    assert result["comparison"]["effect_delta"] == 0
    assert result["scientific_authorization"] == "NOT_ESTABLISHED"
    assert "SYNTHETIC DATA" in " ".join(result["limitations"])
    assert result["data_hashes"]["source_sha256"]
    assert result["data_hashes"]["base_revision_sha256"]
    assert result["provenance"]["algorithm_version"] == analysis.ALGORITHM_VERSION
    assert result["provenance"]["python_version"]
    assert result["provenance"]["software_versions"]["numpy"]
    assert result["provenance"]["software_versions"]["scipy"]
    saved = project.get_run(result["run_id"])
    assert saved["before"] == before
    json.dumps(result, allow_nan=False)


def test_relabel_changes_threshold_and_preserves_original(project):
    base = head(project)
    target = revise(project, ["near-target"], {"label": "Stromal"})
    result = run(project, base, target, min_effect=.2)
    assert result["before"]["excess_over_abundance"] == .5
    assert result["after"]["observed_fraction"] == 0
    assert result["after"]["null_fraction"] == .25
    assert result["after"]["excess_over_abundance"] == -.25
    assert result["after"]["classification"] == "descriptive_not_supported"
    assert result["comparison"]["status"] == "changed"
    assert result["comparison"]["effect_delta"] == -.75
    assert result["comparison"]["label_changes"] == 1
    assert result["comparison"]["source_edge_denominator_delta"] == 0
    assert next(c for c in project.cells(base) if c["cell_id"] == "near-target")["label"] == "Myeloid"


def test_region_revision_can_have_no_metric_effect(project):
    base = head(project)
    target = revise(project, ["near-target"], {"region": "candidate niche"})
    result = run(project, base, target)
    assert result["before"] == result["after"]
    assert result["comparison"]["status"] == "stable"
    assert result["comparison"]["region_changes"] == 1


def test_exclusion_recomputes_edges_and_abundance(tmp_path):
    project = Project.create(tmp_path / "exclude", [
        cell("s", 0, 0, "T cell"), cell("t", 1, 0, "Myeloid"),
        cell("n", 0, 1, "Stromal"), cell("ft", 20, 0, "Myeloid"),
        cell("fn", 30, 0, "Stromal"),
    ], metadata())
    base = head(project)
    target = revise(project, ["n"], {"included": False})
    result = run(project, base, target)
    assert result["before"]["observed_fraction"] == .5
    assert result["before"]["null_fraction"] == .5
    assert result["after"]["observed_fraction"] == 1
    assert result["after"]["null_fraction"] == pytest.approx(2 / 3)
    assert result["after"]["excess_over_abundance"] == pytest.approx(1 / 3)
    assert result["comparison"]["status"] == "changed"
    assert result["comparison"]["source_edge_denominator_delta"] == -1
    assert result["comparison"]["graph_node_count_delta"] == -1
    assert result["comparison"]["exclusion_changes"] == 1
    assert result["after"]["composition"]["selected_included"] == 4
    assert result["after"]["composition"]["selected_excluded"] == 1


def test_frozen_roi_and_whole_slice_are_different_questions(project):
    base = head(project)
    selection = project.set_selection(expected_revision=base, cell_ids=["source", "near-target"])
    induced = run(project, selection_id=selection["selection_id"], graph_scope="roi_induced")
    whole = run(project, selection_id=selection["selection_id"], graph_scope="whole_slice")
    assert induced["before"]["graph"]["node_count"] == 2
    assert induced["before"]["null_fraction"] == 1
    assert induced["before"]["excess_over_abundance"] == 0
    assert whole["before"]["graph"]["node_count"] == 5
    assert whole["before"]["null_fraction"] == .5
    assert whole["before"]["classification"] == "descriptive_supported"
    assert induced["before"]["composition"] == whole["before"]["composition"]
    target = revise(project, ["far-target"], {"label": "Stromal"})
    frozen = run(project, base, target, selection_id=selection["selection_id"], graph_scope="whole_slice")
    assert frozen["comparison"]["selected_cell_count"] == 2
    assert frozen["comparison"]["selected_changes"]["label_changes"] == 0
    assert frozen["comparison"]["label_changes"] == 1
    assert frozen["after"]["null_fraction"] == .25


@pytest.mark.parametrize("source,target,reason", [
    ("Absent", "Myeloid", "source_label_absent_from_included_selection"),
    ("T cell", "Absent", "target_label_absent_from_neighbor_universe"),
])
def test_absent_groups_are_inconclusive(project, source, target, reason):
    result = run(project, source_label=source, target_label=target)
    assert result["before"]["classification"] == "inconclusive"
    assert result["before"]["reason"] == reason
    assert result["comparison"]["status"] == "indeterminate"
    assert result["comparison"]["effect_delta"] is None


def test_disappearing_target_does_not_become_negative_conclusion(project):
    base = head(project)
    target = revise(project, ["near-target", "far-target"], {"label": "Stromal"})
    result = run(project, base, target)
    assert result["before"]["classification"] == "descriptive_supported"
    assert result["after"]["classification"] == "inconclusive"
    assert result["comparison"]["status"] == "indeterminate"


def test_no_eligible_edges_is_inconclusive(project):
    result = run(project, radius_um=.01)
    assert result["before"]["reason"] == "no_eligible_source_neighbor_edges"
    assert result["before"]["graph"]["isolated_sources"] == 1
    assert result["before"]["observed_fraction"] is None


def test_same_label_null_excludes_self(tmp_path):
    project = Project.create(tmp_path / "same-label", [
        cell("a", 0, 0, "Myeloid"), cell("b", 1, 0, "Myeloid"),
        cell("c", 20, 0, "Stromal"), cell("d", 30, 0, "Stromal"),
    ], metadata())
    result = run(project, source_label="Myeloid", target_label="Myeloid")
    assert result["before"]["graph"]["directed_source_edges"] == 2
    assert result["before"]["graph"]["target_edges"] == 2
    assert result["before"]["null_fraction"] == pytest.approx(1 / 3)
    assert result["before"]["excess_over_abundance"] == pytest.approx(2 / 3)


def test_only_one_same_label_target_has_unavailable_null(project):
    result = run(project, source_label="T cell", target_label="T cell")
    assert result["before"]["reason"] == "nonself_target_abundance_null_unavailable"
    assert result["before"]["classification"] == "inconclusive"


def test_distinct_cells_at_same_coordinate_are_neighbors(tmp_path):
    project = Project.create(tmp_path / "duplicates", [
        cell("a", 0, 0, "T cell"), cell("b", 0, 0, "Myeloid"),
        cell("c", 30, 0, "Stromal"),
    ], metadata())
    result = run(project, radius_um=.01)
    assert result["before"]["graph"]["directed_source_edges"] == 1
    assert result["before"]["observed_fraction"] == 1


def test_radius_boundary_is_inclusive(project):
    assert run(project, radius_um=1)["before"]["graph"]["directed_source_edges"] == 1
    assert run(project, radius_um=math.nextafter(1, 0))["before"]["graph"]["directed_source_edges"] == 0


def test_threshold_boundary_is_inclusive(project):
    assert run(project, min_effect=.5)["before"]["classification"] == "descriptive_supported"
    assert run(project, min_effect=math.nextafter(.5, 1))["before"]["classification"] == "descriptive_not_supported"


def test_empty_selection_is_inconclusive(project):
    selection = project.set_selection(expected_revision=head(project), cell_ids=[])
    result = run(project, selection_id=selection["selection_id"])
    assert result["before"]["composition"]["selected_total"] == 0
    assert result["before"]["composition"]["label_fractions"] == {}
    assert result["before"]["classification"] == "inconclusive"


@pytest.mark.parametrize("parameter,value", [
    ("radius_um", True), ("radius_um", "35"), ("radius_um", float("nan")),
    ("radius_um", float("inf")), ("radius_um", 0), ("radius_um", -1),
    ("radius_um", 10 ** 1000),
    ("min_effect", True), ("min_effect", 0), ("min_effect", -.1),
    ("min_effect", float("nan")), ("min_effect", 1.01),
    ("graph_scope", "unknown"), ("graph_scope", []),
    ("source_label", ""), ("target_label", False), ("selection_id", ""),
])
def test_strict_parameters(project, parameter, value):
    with pytest.raises(SpatialError):
        run(project, **{parameter: value})


def test_unknown_selection_and_revision_are_rejected(project):
    with pytest.raises(SpatialError):
        run(project, selection_id="does-not-exist")
    with pytest.raises(SpatialError):
        run(project, before="does-not-exist")


def test_selection_from_neither_revision_is_rejected(project):
    selection = project.set_selection(expected_revision=head(project), cell_ids=["source"])
    middle = revise(project, ["far-other"], {"region": "one"})
    target = revise(project, ["far-other"], {"region": "two"})
    with pytest.raises(SpatialError, match="neither comparison revision"):
        run(project, middle, target, selection_id=selection["selection_id"])


@pytest.mark.parametrize("field,value", [
    ("slice_id", "different-slice"), ("coordinate_system", "different-frame"),
    ("units", "pixel"), ("cell_ids", ["absent-cell"]),
])
def test_selection_mismatch_is_rejected(project, monkeypatch, field, value):
    selection = project.set_selection(expected_revision=head(project), cell_ids=["source"])
    selection[field] = value
    monkeypatch.setattr(project, "get_selection", lambda _id: selection)
    with pytest.raises(SpatialError):
        run(project, selection_id=selection["selection_id"])


def test_project_budget_checked_before_materialization(project, monkeypatch):
    monkeypatch.setattr(analysis, "MAX_PROJECT_CELLS", 4)
    def forbidden(*_args):
        raise AssertionError("cells must not be materialized before project size guard")
    monkeypatch.setattr(project, "cells", forbidden)
    with pytest.raises(SpatialError, match="Project resource budget"):
        run(project)


def test_graph_edge_budget(project, monkeypatch):
    monkeypatch.setattr(analysis, "MAX_DIRECTED_SOURCE_EDGES", 2)
    with pytest.raises(SpatialError, match="directed source edges"):
        run(project, radius_um=100)
    assert project.summary()["runs"] == []


def test_graph_neighbor_budget(project, monkeypatch):
    monkeypatch.setattr(analysis, "MAX_NEIGHBORS_PER_SOURCE", 2)
    with pytest.raises(SpatialError, match="one source has"):
        run(project, radius_um=100)


def test_tree_edges_match_independent_brute_force(tmp_path):
    # Deliberately uneven degrees, duplicate positions and excluded cells.
    cells = [cell(str(i), (i * 7) % 11, (i * 3) % 13,
                  ("T cell", "Myeloid", "Stromal")[i % 3]) for i in range(40)]
    cells[8]["included"] = False
    cells[10]["x"], cells[10]["y"] = cells[0]["x"], cells[0]["y"]
    project = Project.create(tmp_path / "brute-force", cells, metadata())
    selected = {str(i) for i in range(19)}
    selection = project.set_selection(expected_revision=head(project), cell_ids=sorted(selected))
    for scope in ("roi_induced", "whole_slice"):
        result = run(project, selection_id=selection["selection_id"], radius_um=4.25, graph_scope=scope)
        nodes = [c for c in cells if c["included"] and (scope == "whole_slice" or c["cell_id"] in selected)]
        sources = [c for c in nodes if c["label"] == "T cell" and c["cell_id"] in selected]
        edges = [(a, b) for a in sources for b in nodes if a["cell_id"] != b["cell_id"]
                 and math.hypot(a["x"] - b["x"], a["y"] - b["y"]) <= 4.25]
        target_edges = sum(b["label"] == "Myeloid" for _, b in edges)
        assert result["before"]["graph"]["directed_source_edges"] == len(edges)
        assert result["before"]["graph"]["target_edges"] == target_edges
        assert result["before"]["observed_fraction"] == target_edges / len(edges)
        assert result["before"]["null_fraction"] == sum(c["label"] == "Myeloid" for c in nodes) / (len(nodes) - 1)
