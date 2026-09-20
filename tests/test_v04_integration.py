"""Exercise identity, unknown scale, shared plans and offline reconstruction together."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest

from spatial_collab.hypotheses import create_hypothesis_plan, freeze_hypothesis_plan, run_hypothesis_plan
from spatial_collab.identity import qualify_observation_id
from spatial_collab.replay import verify_bundle
from spatial_collab.server import ToolService
from spatial_collab.store import Project, SpatialError, _hash, _json, _validate_import


def make_project(tmp_path, units="array_index"):
    sample = "sample-A"
    cells = [{"cell_id": qualify_observation_id(sample, cid), "sample_id": sample, "source_cell_id": cid,
              "x": x, "y": 0, "label": label, "counts": counts}
             for cid, x, label, counts in [("barcode-1", 0, "T", {"EN1": 2, "EN2": 7}),
                                           ("barcode-2", 1, "B", {"EN1": 10}),
                                           ("barcode-3", 2, "O", {})]]
    return Project.create(tmp_path / "project", cells, {
        "name": "Explicit identity fixture", "slice_id": "slice", "sample_id": sample,
        "identity_scope": "sample_qualified", "coordinate_system": "native_xy", "units": units,
        "panel_genes": ["EN1", "EN2"], "features": [{"feature_id": "EN1", "symbol": "Repeated"},
                                                        {"feature_id": "EN2", "symbol": "Repeated"}],
        "source_kind": "synthetic", "observation_unit": "spot", "label_semantics": "spot_annotation"})


def plan_spec(project):
    rev = project.summary()["head_revision"]
    selected = project.set_selection(rev, cell_ids=[qualify_observation_id("sample-A", "barcode-1")])
    return {"name": "All declared alternatives", "rationale": "Compare working labels without approving them.",
            "revision_id": rev, "selection_id": selected["selection_id"], "source_label": "T", "target_label": "B",
            "radii_um": [1.5, 3.0], "graph_scopes": ["whole_slice", "roi_induced"], "min_effect": 0.1,
            "background": {"kind": "neighbor_universe"}, "variants": [
                {"name": "Alternative", "rationale": "Working alternative", "changes": {"label": "B"},
                 "selector": {"cell_ids": [qualify_observation_id("sample-A", "barcode-3")]}},
                {"name": "Missing marker", "rationale": "Retain absence of a measurement as unknown",
                 "changes": {"included": False},
                 "selector": {"predicate": {"all": [{"field": "counts.UNMEASURED", "op": "gt", "value": 0}]}}}]}


def test_feature_ambiguity_survives_marker_query_region_export(tmp_path):
    project = make_project(tmp_path)
    service = ToolService(project.root)
    rev = project.summary()["head_revision"]
    first = qualify_observation_id("sample-A", "barcode-1")
    selection = project.set_selection(rev, cell_ids=[first])
    summary = service.inspect_selection(["Repeated", "feature_id:EN1", "feature_id:EN2"])
    assert summary["expression"]["Repeated"]["status"] == "ambiguous"
    assert summary["expression"]["feature_id:EN1"]["mean"] == 2
    assert summary["expression"]["feature_id:EN2"]["mean"] == 7
    with pytest.raises(SpatialError, match="Ambiguous"):
        service.query_observations(rev, gene="Repeated")
    query = service.query_observations(rev, gene="feature_id:EN2", min_count=1)
    assert query["observations"][0]["cell_id"] == first
    assert query["observations"][0]["source_cell_id"] == "barcode-1"
    run = service.run_region_comparison(rev, rev, selection["selection_id"], genes=["Repeated", "feature_id:EN1"])
    assert run["before"]["foreground"]["markers"]["Repeated"]["status"] == "ambiguous"
    assert run["before"]["marker_contrast"]["feature_id:EN1"]["mean_raw_count_difference"] == -3
    result = verify_bundle(project.export_bundle(run["run_id"], compact=True)["export_path"])
    assert result["recompute_status"] == "matched"
    with pytest.raises(SpatialError, match="Unknown cell"):
        project.set_selection(rev, cell_ids=[qualify_observation_id("sample-B", "barcode-1")])
    catalog = service.get_feature_catalog("Repeated", limit=1)
    assert catalog["matching_count"] == 2 and catalog["next_offset"] == 1
    assert catalog["features"][0]["symbol_ambiguous"] and not catalog["counts_merged"]


@pytest.mark.parametrize("units", ["array_index", "micrometer"])
@pytest.mark.parametrize("compact", [True, False])
def test_plan_roundtrip_retains_all_unknowns_and_shared_context(tmp_path, units, compact):
    project = make_project(tmp_path, units)
    values = plan_spec(project)
    before = project.context()
    draft = create_hypothesis_plan(project, values)
    after = project.context(before["cursor"])
    assert after["changed"] and after["head_revision"] == before["head_revision"]
    assert after["hypotheses"]["latest"]["version"] == 1
    frozen = freeze_hypothesis_plan(project, draft["plan_id"], 1)
    run = run_hypothesis_plan(project, frozen["plan_id"], frozen["version"])
    assert len(run["results"]) == len(run["difference_order"]) == 8
    assert run["status_counts"]["unknown"] >= 4
    folder = Path(project.export_bundle(run["run_id"], compact=compact)["export_path"])
    manifest = json.loads((folder / "manifest.json").read_text())
    assert ("base_cells.json" not in manifest["files"]) is compact
    result = verify_bundle(folder)
    assert result["recompute_status"] == "matched"
    assert "difference_order" in result["compared_fields"]
    assert "status_counts" in result["compared_fields"]
    assert project.context()["head_revision"] == before["head_revision"]


def test_compact_without_result_and_conflicting_layout_rejected(tmp_path):
    project = make_project(tmp_path)
    folder = Path(project.export_bundle(compact=True)["export_path"])
    assert verify_bundle(folder)["recompute_status"] == "no_result_in_bundle"
    manifest = json.loads((folder / "manifest.json").read_text())
    manifest.pop("layout")
    (folder / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(SpatialError, match="incomplete"):
        verify_bundle(folder)


def test_host_metadata_is_bounded_but_internal_feature_mapping_complete(tmp_path):
    project = make_project(tmp_path)
    service = ToolService(project.root)
    state = service.open_project()
    assert "features" not in state["metadata"]
    assert state["metadata"]["feature_catalog"]["duplicate_symbol_count"] == 1
    assert state["metadata"]["panel_gene_count"] == 2
    assert len(project.summary()["metadata"]["features"]) == 2
    assert not state["coordinate_capabilities"]["physical_radius_analysis"]
    assert service.list_assets()["missing_evidence_status"] == "UNKNOWN"


def test_legacy_source_accepts_explicit_catalog_queries_without_rewriting_source(tmp_path):
    source = make_project(tmp_path)
    meta = deepcopy(source.summary()["metadata"])
    meta.pop("features")
    legacy = Project.create(tmp_path / "legacy", source.cells(), meta)
    service = ToolService(legacy.root)
    before = legacy.summary()
    rev = before["head_revision"]
    selection = legacy.set_selection(rev, cell_ids=[qualify_observation_id("sample-A", "barcode-1")])
    summary = service.inspect_selection(["EN1", "feature_id:EN1", "symbol:EN1"])
    assert summary["expression"]["EN1"]["mean"] == 2
    for query in ("feature_id:EN1", "symbol:EN1"):
        assert summary["expression"][query]["status"] == "measured"
        assert summary["expression"][query]["mean"] == 2
        assert len(service.query_observations(rev, gene=query, min_count=1)["observations"]) == 2
    assert summary["cells"][0]["counts"] == {"EN1": 2}
    run = service.run_region_comparison(rev, rev, selection["selection_id"], genes=["feature_id:EN1"])
    assert run["before"]["marker_contrast"]["feature_id:EN1"]["mean_raw_count_difference"] == -3
    assert verify_bundle(legacy.export_bundle(run["run_id"], compact=True)["export_path"])["recompute_matches"]
    assert legacy.summary()["source_sha256"] == before["source_sha256"]
    assert "features" not in legacy.summary()["metadata"]


@pytest.mark.parametrize("case", ["metadata_mismatch", "mixed_namespace", "missing_sample"])
def test_qualified_identity_cannot_disagree_with_declared_sample(tmp_path, case):
    project = make_project(tmp_path)
    meta, cells = deepcopy(project.summary()["metadata"]), project.cells()
    if case == "metadata_mismatch":
        meta["sample_id"] = "sample-B"
    elif case == "missing_sample":
        meta.pop("sample_id")
    else:
        cells[0]["sample_id"] = "sample-B"
        cells[0]["cell_id"] = qualify_observation_id("sample-B", cells[0]["source_cell_id"])
    with pytest.raises(SpatialError, match="sample_id"):
        _validate_import(cells, meta)


@pytest.mark.parametrize("recompute", [True, False])
def test_rehashed_plan_bundle_cannot_rebind_outer_result_to_another_roi(tmp_path, recompute):
    project = make_project(tmp_path, "micrometer")
    values = plan_spec(project)
    draft = create_hypothesis_plan(project, values)
    frozen = freeze_hypothesis_plan(project, draft["plan_id"], 1)
    other = project.set_selection(values["revision_id"], cell_ids=[qualify_observation_id("sample-A", "barcode-2")])
    run = run_hypothesis_plan(project, frozen["plan_id"], frozen["version"])
    folder = Path(project.export_bundle(run["run_id"], compact=True)["export_path"])
    assert verify_bundle(folder)["recompute_matches"]
    result = json.loads((folder / "result.json").read_text())
    result["selection_id"] = other["selection_id"]
    provenance = json.loads((folder / "provenance.json").read_text())
    provenance["run_sha256"] = _hash(result)
    manifest = json.loads((folder / "manifest.json").read_text())
    for name, content in {"result.json": result, "run_selection.json": other, "provenance.json": provenance}.items():
        data = (_json(content) + "\n").encode()
        (folder / name).write_bytes(data)
        manifest["files"][name] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    (folder / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(SpatialError, match="anchors.*frozen plan"):
        verify_bundle(folder, recompute=recompute)
