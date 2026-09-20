from copy import deepcopy
import json
import sqlite3

import pytest

from spatial_collab import hypotheses as h
from spatial_collab.store import Project, SpatialError, _hash


@pytest.fixture
def project(tmp_path):
    return Project.create(tmp_path / "project", [
        {"cell_id": "s1", "x": 0, "y": 0, "label": "T", "counts": {"G": 2}, "attributes": {"nucleus_count": 1}},
        {"cell_id": "b1", "x": 1, "y": 0, "label": "B", "counts": {}},
        {"cell_id": "o1", "x": 2, "y": 0, "label": "O", "counts": {}, "attributes": {"nucleus_count": 2}},
        {"cell_id": "o2", "x": 3, "y": 0, "label": "O", "counts": {}},
        {"cell_id": "s2", "x": 10, "y": 0, "label": "T", "counts": {}},
        {"cell_id": "b2", "x": 11, "y": 0, "label": "B", "counts": {}},
        {"cell_id": "isolated", "x": 30, "y": 0, "label": "T", "counts": {}},
    ], {"name": "weighting test", "slice_id": "slice", "coordinate_system": "native_xy",
        "units": "micrometer", "source_kind": "synthetic", "panel_genes": ["G"]})


def spec(project):
    rev = project.summary()["head_revision"]
    selection = project.set_selection(rev, cell_ids=["s1", "s2", "isolated"])
    return {"name": "Declared alternatives", "rationale": "Inspect how source weighting changes a descriptive contrast.",
            "revision_id": rev, "selection_id": selection["selection_id"], "source_label": "T", "target_label": "B",
            "radii_um": [3.0], "graph_scopes": ["whole_slice"], "min_effect": 0.25,
            "background": {"kind": "neighbor_universe"}, "variants": [
                {"name": "relabel", "rationale": "A competing working annotation.",
                 "selector": {"cell_ids": ["o1"]}, "changes": {"label": "B"}},
                {"name": "exclude", "rationale": "Independent source exclusion assumption.",
                 "selector": {"cell_ids": ["s1"]}, "changes": {"included": False}},
            ]}


def freeze(project, specification=None):
    draft = h.create_hypothesis_plan(project, specification or spec(project))
    return h.freeze_hypothesis_plan(project, draft["plan_id"], draft["version"])


def test_versions_cas_frozen_history_and_scientific_head_unchanged(project):
    before = project.summary()
    values = spec(project)
    draft = h.create_hypothesis_plan(project, values)
    assert draft["status"] == "draft" and draft["version"] == 1
    with pytest.raises(SpatialError, match="frozen"):
        h.run_hypothesis_plan(project, draft["plan_id"], 1)
    frozen = h.freeze_hypothesis_plan(project, draft["plan_id"], 1)
    edited = deepcopy(values)
    edited["radii_um"] = [2.0, 3.0]
    revised = h.revise_hypothesis_plan(project, draft["plan_id"], 2, edited)
    assert revised["status"] == "draft" and revised["version"] == 3
    with pytest.raises(SpatialError, match="Stale"):
        h.revise_hypothesis_plan(project, draft["plan_id"], 2, values)
    with pytest.raises(SpatialError, match="Stale"):
        h.freeze_hypothesis_plan(project, draft["plan_id"], 2)
    assert h.get_hypothesis_plan(project, draft["plan_id"], 2) == frozen
    result = h.run_hypothesis_plan(project, frozen["plan_id"], 2)
    assert len(result["results"]) == 2
    assert project.summary()["head_revision"] == before["head_revision"]
    assert project.summary()["revision_count"] == before["revision_count"]
    assert result["scientific_authorization"] == "NOT_ESTABLISHED" and not result["revision_created"]
    assert len(h.list_hypothesis_plans(project)["plans"][0]["versions"]) == 3
    with sqlite3.connect(project.root / "hypotheses.sqlite3") as db:
        with pytest.raises(sqlite3.IntegrityError, match="Immutable"):
            db.execute("UPDATE versions SET payload='{}'")


def test_equal_source_reference_differs_for_unequal_degree_and_retains_isolation(project):
    frozen = freeze(project)
    result = h.run_hypothesis_plan(project, frozen_plan=frozen)
    row = result["results"][0]
    before = row["before"]
    edge, equal = (before["methods"][name] for name in h.METHODS)
    assert before["graph"]["source_count"] == 3
    assert before["graph"]["isolated_sources"] == 1
    assert before["graph"]["target_edges"] == 2
    assert before["graph"]["directed_source_edges"] == 4
    assert edge["observed_fraction"] == pytest.approx(0.5)
    assert equal["observed_fraction"] == pytest.approx(2 / 3)
    assert edge["null_fraction"] == equal["null_fraction"] == pytest.approx(1 / 3)
    assert edge["denominator"] == 4 and equal["denominator"] == 2
    assert edge["classification"] == "descriptive_not_supported"
    assert equal["classification"] == "descriptive_supported"
    assert row["reference_comparison"]["before"]["source_equal_minus_edge_weighted"] == pytest.approx(1 / 6)
    # Larger absolute difference may be inspected first; no rows are omitted.
    assert result["difference_order"] == ["row_0002", "row_0001"]
    assert set(result["difference_order"]) == {r["row_id"] for r in result["results"]}


def test_selected_background_excludes_self_per_source_with_matching_method_weights(project):
    values = spec(project)
    bg = project.set_selection(values["revision_id"], cell_ids=["s1", "b1", "b2"])
    values["background"] = {"kind": "selection", "selection_id": bg["selection_id"]}
    result = h.run_hypothesis_plan(project, frozen_plan=freeze(project, values))
    before = result["results"][0]["before"]
    assert before["graph"]["background_count"] == 3
    assert before["methods"]["edge_weighted"]["null_fraction"] == pytest.approx(11 / 12)
    assert before["methods"]["source_equal_weighted"]["null_fraction"] == pytest.approx(5 / 6)


def test_missing_predicate_values_stay_unknown_and_do_not_erase_other_variants(project):
    values = spec(project)
    values["variants"].append({"name": "missing attribute", "rationale": "Unknown nuclei require follow-up.",
                               "selector": {"predicate": {"all": [
                                   {"field": "attributes.nucleus_count", "op": "ne", "value": 1}]}},
                               "changes": {"included": False}})
    frozen = freeze(project, values)
    resolution = frozen["resolved"]["variants"][2]
    assert resolution["cell_ids"] == ["o1"] and resolution["unknown_count"] == 5
    result = h.run_hypothesis_plan(project, frozen_plan=frozen)
    assert result["status_counts"] == {"computed": 2, "unknown": 1}
    assert result["results"][2]["after"]["reason"] == "hypothesis_selector_has_unknown_observations"
    assert result["results"][2]["comparison"]["edge_weighted"]["effect_delta"] is None
    assert result["difference_order"][-1] == "row_0003"


def test_unmeasured_feature_is_unknown_but_measured_zero_and_false_dominance_are_explicit(project):
    values = spec(project)
    values["variants"] = [
        {"name": "absent panel", "rationale": "Test absent feature evidence.", "selector": {"predicate": {"all": [
            {"field": "counts.NOT_MEASURED", "op": "eq", "value": 0}]}}, "changes": {"label": "U"}},
        {"name": "measured zeros", "rationale": "Zero in a measured feature.", "selector": {"predicate": {"all": [
            {"field": "counts.G", "op": "eq", "value": 0}]}}, "changes": {"label": "U"}},
        {"name": "not applicable", "rationale": "False clause makes missing clause irrelevant.", "selector": {"predicate": {"all": [
            {"field": "label", "op": "eq", "value": "NOT_PRESENT"},
            {"field": "attributes.missing", "op": "ne", "value": 1}]}}, "changes": {"label": "U"}},
    ]
    frozen = freeze(project, values)
    resolved = frozen["resolved"]["variants"]
    assert resolved[0]["unknown_count"] == 7
    assert resolved[1]["matched_count"] == 6 and resolved[1]["unknown_count"] == 0
    assert resolved[2]["matched_count"] == resolved[2]["unknown_count"] == 0
    result = h.run_hypothesis_plan(project, frozen_plan=frozen)
    assert result["results"][2]["status"] == "computed"
    assert result["results"][2]["comparison"]["edge_weighted"]["effect_delta"] == 0


@pytest.mark.parametrize(("actual", "op", "expected", "truth"), [
    ("0", "ne", 0, None), ("0", "eq", 0, None),
    (True, "ne", 1, None), (False, "eq", 0, None),
    (1, "eq", 1.0, True), (1.0, "ne", 1, False),
    (True, "eq", True, True), (False, "ne", True, True),
    ("0", "ne", "1", True), ("0", "eq", "0", True),
    (None, "ne", 0, None), (1, "ne", None, None),
    (1, "in", ["1", 1.0], True), (1, "in", [2, "1"], None),
    (1, "in", [2, 3.0], False), (True, "in", [1], None),
    (True, "in", [1, True], True), (False, "in", [True], False),
    ("0", "in", [0, "1"], None), ("0", "in", [0, "0"], True),
])
def test_predicate_comparisons_preserve_unknown_type_compatibility(actual, op, expected, truth):
    clause = {"field": "attributes.qc", "op": op, "value": expected}
    assert h._clause({"attributes": {"qc": actual}}, clause, set()) is truth


def test_type_incompatible_qc_values_freeze_as_unknown_and_never_exclude(tmp_path):
    project = Project.create(tmp_path / "typed-qc", [
        {"cell_id": "string", "x": 0, "y": 0, "label": "T", "counts": {}, "attributes": {"nucleus_count": "0"}},
        {"cell_id": "boolean", "x": 1, "y": 0, "label": "B", "counts": {}, "attributes": {"nucleus_count": True}},
        {"cell_id": "number", "x": 2, "y": 0, "label": "B", "counts": {}, "attributes": {"nucleus_count": 0}},
        {"cell_id": "missing", "x": 3, "y": 0, "label": "B", "counts": {}},
    ], {"name": "Mixed QC types", "slice_id": "slice", "coordinate_system": "native_xy",
        "units": "micrometer", "source_kind": "synthetic", "panel_genes": ["G"]})
    revision = project.summary()["head_revision"]
    selection = project.set_selection(revision, cell_ids=["string"])
    values = {"name": "QC alternative", "rationale": "Do not coerce imported QC metadata.",
              "revision_id": revision, "selection_id": selection["selection_id"],
              "source_label": "T", "target_label": "B", "radii_um": [3],
              "graph_scopes": ["whole_slice"], "min_effect": 0.1,
              "background": {"kind": "neighbor_universe"}, "variants": [{
                  "name": "non-single nucleus", "rationale": "A conditional QC exclusion.",
                  "selector": {"predicate": {"all": [{"field": "attributes.nucleus_count", "op": "ne", "value": 1}]}},
                  "changes": {"included": False}}]}
    frozen = freeze(project, values)
    resolved = frozen["resolved"]["variants"][0]
    assert resolved["cell_ids"] == ["number"]
    assert resolved["unknown_cell_ids"] == ["boolean", "missing", "string"]
    result = h.run_hypothesis_plan(project, frozen_plan=frozen)
    assert result["status_counts"] == {"unknown": 1}
    assert result["results"][0]["after"]["reason"] == "hypothesis_selector_has_unknown_observations"
    assert all(cell["included"] for cell in project.cells())
    assert project.summary()["head_revision"] == revision


def test_resource_failures_preserved_per_radius_other_grid_cells_continue(project, monkeypatch):
    values = spec(project)
    values["radii_um"] = [1.1, 3.0]
    monkeypatch.setattr(h, "MAX_DIRECTED_SOURCE_EDGES", 3)
    result = h.run_hypothesis_plan(project, frozen_plan=freeze(project, values))
    assert len(result["results"]) == 4
    assert [r["status"] for r in result["results"]] == ["computed", "failed", "computed", "failed"]
    assert result["results"][1]["before"]["reason"] == "resource_limit"
    assert len(result["difference_order"]) == 4


def test_freeze_unknown_units_keeps_every_grid_cell_unknown_without_evaluating_graph(project, monkeypatch):
    # Adapter represents an uncalibrated imported project while older Project
    # contracts are also exercised by this test suite.
    class Uncalibrated:
        root = project.root

        def summary(self):
            result = project.summary()
            result["metadata"]["units"] = "unknown"
            return result

        def get_selection(self, sid):
            result = project.get_selection(sid)
            result["units"] = "unknown"
            return result

        get_revision = staticmethod(project.get_revision)
        cells = staticmethod(project.cells)
        save_run = staticmethod(project.save_run)

    values = spec(project)
    adapter = Uncalibrated()
    frozen = freeze(adapter, values)
    monkeypatch.setattr(h, "_evaluate_v2", lambda *a, **k: pytest.fail("uncalibrated geometry was evaluated"))
    result = h.run_hypothesis_plan(adapter, frozen_plan=frozen)
    assert result["status_counts"] == {"unknown": 2}
    assert all(r["before"]["reason"] == "physical_coordinate_calibration_unavailable" for r in result["results"])


def test_frozen_tampering_and_wrong_revision_background_rejected(project):
    frozen = freeze(project)
    corrupted = deepcopy(frozen)
    corrupted["resolved"]["variants"][0]["cell_ids"] = ["b1"]
    with pytest.raises(SpatialError, match="integrity"):
        h.run_hypothesis_plan(project, frozen_plan=corrupted)
    corrupted["plan_sha256"] = _hash({k: v for k, v in corrupted.items() if k != "plan_sha256"})
    with pytest.raises(SpatialError, match="do not match"):
        h.run_hypothesis_plan(project, frozen_plan=corrupted)
    values = deepcopy(frozen["spec"])
    values["variants"][0]["selector"] = {"cell_ids": ["missing"]}
    draft = h.create_hypothesis_plan(project, values)
    with pytest.raises(SpatialError, match="unknown observations"):
        h.freeze_hypothesis_plan(project, draft["plan_id"], draft["version"])
    assert h.get_hypothesis_plan(project, draft["plan_id"])["status"] == "draft"


def test_embedded_frozen_plan_recomputes_without_ledger_or_project_root(project):
    from spatial_collab.replay import _BundleProject

    frozen = freeze(project)
    run = h.run_hypothesis_plan(project, frozen_plan=frozen)
    folder = project.export_bundle(run["run_id"])["export_path"]
    from pathlib import Path
    files = {p.name: json.loads(p.read_text()) for p in Path(folder).glob("*.json")}
    adapter = _BundleProject(files["source_snapshot.json"],
                             {r["revision_id"]: r for r in files["revisions.json"]},
                             {s["selection_id"]: s for s in files["selections.json"]}, files["provenance.json"])
    replayed = h.run_hypothesis_plan(adapter, frozen_plan=run["parameters"]["frozen_plan"])
    for key in ("parameters", "results", "difference_order", "status_counts", "data_hashes", "selection", "methods"):
        assert replayed[key] == run[key]


def test_context_detects_draft_freeze_revise_without_loading_counts(project, monkeypatch):
    empty = h.hypothesis_context(project)
    assert empty["latest"] is None
    values = spec(project)
    draft = h.create_hypothesis_plan(project, values)
    first = h.hypothesis_context(project)
    frozen = h.freeze_hypothesis_plan(project, draft["plan_id"], draft["version"])
    second = h.hypothesis_context(project)
    revised = h.revise_hypothesis_plan(project, frozen["plan_id"], frozen["version"], values)
    monkeypatch.setattr(project, "summary", lambda: pytest.fail("context loaded source payload"))
    monkeypatch.setattr(project, "cells", lambda *args: pytest.fail("context materialized counts"))
    third = h.hypothesis_context(project)
    assert len({x["state_sha256"] for x in (empty, first, second, third)}) == 4
    assert third["latest"] == {"plan_id": revised["plan_id"], "version": 3,
                               "status": "draft", "plan_sha256": revised["plan_sha256"]}


def test_selected_background_observed_zero_is_not_mislabeled_missing(project):
    values = spec(project)
    values["source_label"] = values["target_label"] = "O"
    selection = project.set_selection(values["revision_id"], cell_ids=["o1"])
    background = project.set_selection(values["revision_id"], cell_ids=["o1", "s1"])
    values["selection_id"] = selection["selection_id"]
    values["background"] = {"kind": "selection", "selection_id": background["selection_id"]}
    result = h.run_hypothesis_plan(project, frozen_plan=freeze(project, values))
    before = result["results"][0]["before"]
    assert before["status"] == "computed"
    assert before["graph"]["target_count"] == 2
    assert before["methods"]["edge_weighted"]["null_fraction"] == 0
    assert before["methods"]["source_equal_weighted"]["null_fraction"] == 0


def test_background_containing_only_source_has_undefined_denominator(project):
    values = spec(project)
    background = project.set_selection(values["revision_id"], cell_ids=["s1"])
    values["background"] = {"kind": "selection", "selection_id": background["selection_id"]}
    result = h.run_hypothesis_plan(project, frozen_plan=freeze(project, values))
    before = result["results"][0]["before"]
    assert before["status"] == "unknown"
    assert before["graph"]["sources_with_unavailable_background"] == 1


@pytest.mark.parametrize("mutation", [
    lambda s: s.update(min_effect=True),
    lambda s: s.update(radii_um=[0]),
    lambda s: s.update(graph_scopes=["whole_slice", "whole_slice"]),
    lambda s: s.update(graph_scopes=[{}]),
    lambda s: s.update(background={"kind": []}),
    lambda s: s["variants"][0].update(changes={"counts": {"G": 9}}),
    lambda s: s["variants"][0].update(selector={"predicate": {"python": "True"}}),
])
def test_invalid_specs_are_rejected_before_ledger_write(project, mutation):
    values = spec(project)
    mutation(values)
    with pytest.raises(SpatialError):
        h.create_hypothesis_plan(project, values)
    assert not (project.root / "hypotheses.sqlite3").exists()
