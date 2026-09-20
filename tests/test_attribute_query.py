"""Attribute query semantics: missing evidence, exact types and pinned cursors."""
from __future__ import annotations

import asyncio
from copy import deepcopy

import pytest

from spatial_collab.exploration import query_observations
from spatial_collab.server import ToolService, create_server
from spatial_collab.store import Project, SpatialError


@pytest.fixture
def project(tmp_path):
    values = [
        ("a", "A", True, 0, {"qc": 0, "method": "boundary", "flag": False}),
        ("b", "A", True, 2, {"qc": 2, "method": "expansion", "flag": True}),
        ("c", "B", False, 5, {"qc": 2.0, "method": "boundary", "flag": False}),
        ("d", "A", True, 5, {"qc": 5, "method": "boundary", "flag": True}),
        ("e", "A", True, 5, {"qc": None, "method": None, "flag": None}),
        ("f", "A", True, 5, {}),
        ("g", "A", True, 5, {"qc": "2", "method": 2, "flag": 1}),
        ("h", "A", True, 5, {"qc": True, "method": False, "flag": "true"}),
    ]
    cells = [{"cell_id": identifier, "x": index, "y": index, "label": label,
              "included": included, "counts": {"G": marker}, "attributes": attributes}
             for index, (identifier, label, included, marker, attributes) in enumerate(values)]
    return Project.create(tmp_path / "project", cells, {
        "name": "Attribute filter test", "slice_id": "synthetic", "source_kind": "synthetic",
        "coordinate_system": "native_xy", "units": "micrometer", "panel_genes": ["G"],
    })


def head(project):
    return project.summary()["head_revision"]


def query(project, field="qc", operator="eq", value=2, **options):
    return query_observations(project, options.pop("revision_id", head(project)),
                              attribute_filter={"field": field, "operator": operator, "value": value}, **options)


def identifiers(result):
    return [row["cell_id"] for row in result["observations"]]


@pytest.mark.parametrize("operator,expected", [
    ("eq", ["b", "c"]), ("ne", ["a", "d"]), ("lt", ["a"]),
    ("le", ["a", "b", "c"]), ("gt", ["d"]), ("ge", ["b", "c", "d"]),
])
def test_numeric_operator_truth_table_and_unknowns(project, operator, expected):
    result = query(project, operator=operator)
    assert identifiers(result) == expected
    assert result["attribute_unavailable_count"] == 4
    assert result["matching_count"] == len(expected)
    assert result["complete"] and result["next_cursor"] is None
    assert all(row["matched_attribute"]["qc"] in (0, 2, 5) for row in result["observations"])


def test_missing_null_and_incompatible_are_not_true_for_ne(project):
    result = query(project, operator="ne", value=2)
    assert identifiers(result) == ["a", "d"]
    assert result["attribute_unavailable_count"] == 4
    assert "unknown" in result["attribute_interpretation"]
    assert query(project, field="not_recorded", operator="ne", value=0)["attribute_unavailable_count"] == 8
    assert query(project, field="not_recorded", operator="ne", value=0)["matching_count"] == 0


def test_numeric_strings_and_booleans_are_not_coerced(project):
    string = query(project, value="2")
    assert identifiers(string) == ["g"] and string["attribute_unavailable_count"] == 7
    boolean = query(project, value=True)
    assert identifiers(boolean) == ["h"] and boolean["attribute_unavailable_count"] == 7
    assert identifiers(query(project, field="flag", value=True)) == ["b", "d"]
    assert identifiers(query(project, field="flag", value=False)) == ["a", "c"]
    assert identifiers(query(project, field="flag", value=1)) == ["g"]
    assert identifiers(query(project, operator="eq", value=2.0)) == ["b", "c"]


def test_string_equality_and_inequality_are_exact(project):
    assert identifiers(query(project, field="method", value="boundary")) == ["a", "c", "d"]
    assert identifiers(query(project, field="method", operator="ne", value="boundary")) == ["b"]
    assert query(project, field="method", value="Boundary")["matching_count"] == 0


@pytest.mark.parametrize("attribute_filter", [
    [], "qc=2", 1, {}, {"field": "qc", "operator": "eq"},
    {"field": "qc", "operator": "eq", "value": 2, "extra": True},
    {"field": "", "operator": "eq", "value": 2},
    {"field": "  ", "operator": "eq", "value": 2},
    {"field": "a" * 129, "operator": "eq", "value": 2},
    {"field": 1, "operator": "eq", "value": 2},
    {"field": "qc", "operator": "contains", "value": 2},
    {"field": "qc", "operator": "EQ", "value": 2},
    {"field": "qc", "operator": [], "value": 2},
    {"field": "qc", "operator": {}, "value": 2},
    {"field": "qc", "operator": None, "value": 2},
    {"field": "qc", "operator": True, "value": 2},
    {"field": "qc", "operator": "eq", "value": None},
    {"field": "qc", "operator": "eq", "value": []},
    {"field": "qc", "operator": "eq", "value": {}},
    {"field": "qc", "operator": "eq", "value": float("nan")},
    {"field": "qc", "operator": "eq", "value": float("inf")},
    {"field": "qc", "operator": "eq", "value": -float("inf")},
    {"field": "qc", "operator": "lt", "value": "2"},
    {"field": "qc", "operator": "ge", "value": True},
])
def test_invalid_attribute_filters_raise_domain_error(project, attribute_filter):
    with pytest.raises(SpatialError):
        query_observations(project, head(project), attribute_filter=attribute_filter)


def test_combined_filter_intersection_and_unknown_denominator(project):
    result = query(project, operator="ge", value=2, labels=["A"], included=True,
                   bounds=[0, 0, 5, 5], gene="G", min_count=2)
    assert identifiers(result) == ["b", "d"]
    # e/f passed all preceding filters, but their QC evidence is unavailable.
    assert result["attribute_unavailable_count"] == 2
    assert [row["marker_count"] for row in result["observations"]] == [2, 5]
    assert query(project, operator="ge", value=2, labels=[])["attribute_unavailable_count"] == 0
    with pytest.raises(SpatialError, match="not measured"):
        query(project, gene="UNKNOWN", min_count=0)


def test_query_does_not_change_selection_revision_context_or_input_filter(project):
    revision = head(project)
    selection = project.set_selection(revision, cell_ids=["e", "f"], name="Keep selection")
    context = project.context()
    specification = {"field": "qc", "operator": "ge", "value": 0}
    saved_specification = deepcopy(specification)
    result = query_observations(project, revision, attribute_filter=specification)
    assert result["selection_changed"] is False
    assert specification == saved_specification
    assert project.get_selection()["selection_id"] == selection["selection_id"]
    assert project.context()["cursor"] == context["cursor"]
    assert head(project) == revision


def test_cursor_binding_covers_entire_attribute_predicate(project):
    revision = head(project)
    first = query(project, operator="ge", value=0, limit=1)
    assert identifiers(first) == ["a"]
    for specification in (
        {"field": "qc", "operator": "gt", "value": 0},
        {"field": "qc", "operator": "ge", "value": 1},
        {"field": "flag", "operator": "ge", "value": 0}, None,
    ):
        with pytest.raises(SpatialError, match="Cursor"):
            query_observations(project, revision, cursor=first["next_cursor"], attribute_filter=specification)
    # Dictionary key ordering does not change the query; page size can change.
    second = query_observations(project, revision, limit=2, cursor=first["next_cursor"],
                                 attribute_filter={"value": 0, "operator": "ge", "field": "qc"})
    third = query(project, revision_id=revision, operator="ge", value=0, limit=2, cursor=second["next_cursor"])
    assert identifiers(first) + identifiers(second) + identifiers(third) == ["a", "b", "c", "d"]
    assert third["complete"]


def test_cursor_remains_pinned_across_a_committed_revision(project):
    revision = head(project)
    first = query(project, operator="ge", value=0, included=True, limit=1)
    selection = project.set_selection(revision, cell_ids=["b"])
    proposal = project.propose_revision(revision, selection["selection_id"], {"included": False}, "Synthetic exclusion")
    updated = project.apply_revision(proposal["proposal_id"], revision, "Synthetic reviewer", True)["revision_id"]
    historic = query(project, revision_id=revision, operator="ge", value=0, included=True, cursor=first["next_cursor"])
    assert identifiers(historic) == ["b", "d"] and historic["historical"]
    assert identifiers(query(project, revision_id=updated, operator="ge", value=0, included=True)) == ["a", "d"]
    with pytest.raises(SpatialError, match="Cursor"):
        query(project, revision_id=updated, operator="ge", value=0, included=True, cursor=first["next_cursor"])


def test_service_dispatch_preserves_attribute_types_and_rejects_bad_contract(project):
    service = ToolService(project.root)
    revision = head(project)
    result = service.call("query_observations", {"revision_id": revision,
        "attribute_filter": {"field": "flag", "operator": "eq", "value": True}})
    assert identifiers(result) == ["b", "d"]
    for bad in ([{"field": "qc", "operator": "eq", "value": 2}], "qc=2", True):
        with pytest.raises(SpatialError):
            service.call("query_observations", {"revision_id": revision, "attribute_filter": bad})
    with pytest.raises(SpatialError):
        service.call("query_observations", {"revision_id": revision, "attribute_filters": []})


def test_mcp_schema_and_errors_expose_same_read_only_attribute_contract(project):
    server = create_server(project.root)
    tool = next(tool for tool in asyncio.run(server.list_tools()) if tool.name == "query_observations")
    assert "attribute_filter" in tool.inputSchema["properties"]
    assert "attribute_filter" not in tool.inputSchema.get("required", [])
    assert tool.annotations.readOnlyHint is True
    result = asyncio.run(server.call_tool("query_observations", {
        "revision_id": head(project), "attribute_filter": {"field": "qc", "operator": "ne", "value": 2}}))
    assert result.isError is False
    assert identifiers(result.structuredContent) == ["a", "d"]
    invalid = asyncio.run(server.call_tool("query_observations", {
        "revision_id": head(project), "attribute_filter": {"field": "qc", "operator": [], "value": 2}}))
    assert invalid.isError is True
    assert "operator" in invalid.structuredContent["error"]
