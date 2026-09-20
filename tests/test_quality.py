from copy import deepcopy
import json
import math

import pytest

from spatial_collab.quality import inspect_quality
from spatial_collab.replay import verify_bundle
from spatial_collab.store import Project, SpatialError, _validate_import


@pytest.fixture
def metadata():
    return {"name": "QC fixture", "slice_id": "s", "coordinate_system": "s-xy",
            "source_kind": "synthetic", "units": "micrometer", "panel_genes": ["G", "H"]}


@pytest.fixture
def cells():
    return [
        {"cell_id": "a", "x": 1, "y": 2, "label": "A", "counts": {"G": 2, "H": 3},
         "attributes": {"transcript_counts": 50, "total_counts": 60, "control_probe_counts": 0,
                        "genomic_control_counts": 2, "cell_area": 10.0, "nucleus_area": 2.0,
                        "nucleus_count": 1, "segmentation_method": "nucleus expansion",
                        "vendor_note": "source only", "flag": False}},
        {"cell_id": "b", "x": 2, "y": 2, "label": "A", "included": False, "counts": {},
         "attributes": {"transcript_counts": None, "control_probe_counts": 4, "cell_area": 30.0,
                        "nucleus_count": 0, "segmentation_method": "cell boundary"}},
        {"cell_id": "c", "x": 3, "y": 2, "label": "B", "counts": {"G": 1}},
    ]


@pytest.fixture
def project(tmp_path, metadata, cells):
    return Project.create(tmp_path / "project", cells, metadata)


def test_source_attributes_roundtrip_overlay_immutable_and_replay(project):
    rev = project.summary()["head_revision"]
    selected = project.set_selection(rev, cell_ids=["a"])
    preview = project.propose_revision(rev, selected["selection_id"], {"included": False}, "fixture exclusion")
    target = project.apply_revision(preview["proposal_id"], rev, "Test driver", True)["revision_id"]
    assert project.cells(rev)[0]["attributes"] == project.cells(target)[0]["attributes"]
    with pytest.raises(SpatialError, match="only label"):
        project.propose_revision(target, selected["selection_id"], {"attributes": {}}, "not an overlay")
    bundle = project.export_bundle()
    receipt = verify_bundle(bundle["export_path"])
    assert receipt["integrity_verified"]


def test_legacy_without_attributes_keeps_source_shape_and_hash(tmp_path, metadata):
    cells = [{"cell_id": "old", "x": 1.0, "y": 2.0, "label": "A", "counts": {"G": 1.0}}]
    first = Project.create(tmp_path / "legacy", cells, metadata)
    assert "attributes" not in first.cells()[0]
    canonical = _validate_import(cells, metadata)
    assert "attributes" not in canonical["cells"][0]
    copied = Project.create(tmp_path / "copy", first.cells(), first.summary()["metadata"])
    assert first.summary()["source_sha256"] == copied.summary()["source_sha256"]
    qc = inspect_quality(first, first.summary()["head_revision"])
    assert qc["recorded_numeric_fields"]["total_counts"]["status"] == "not_recorded"
    assert qc["recorded_numeric_fields"]["total_counts"]["zero_count"] is None


@pytest.mark.parametrize("attributes", [
    None, [], {str(i): i for i in range(33)}, {"": 1}, {"x" * 129: 1},
    {"note": "x" * 1001}, {"value": float("inf")}, {"value": float("nan")},
    {"value": 10 ** 1000}, {"nested": {}}, {"nested": [1]}, {1: "bad key"},
])
def test_attributes_are_bounded_json_scalars(tmp_path, metadata, attributes):
    with pytest.raises(SpatialError):
        Project.create(tmp_path / "invalid", [{"cell_id": "a", "x": 0, "y": 0,
                                              "label": "A", "attributes": attributes}], metadata)


def test_missing_null_zero_and_imported_feature_denominators_are_distinct(project):
    before = project.context()
    result = inspect_quality(project, before["head_revision"])
    transcript = result["recorded_numeric_fields"]["transcript_counts"]
    assert transcript["valid_count"] == 1 and transcript["mean"] == 50
    assert transcript["absent_count"] == 1 and transcript["null_count"] == 1
    assert transcript["missing_count"] == 2 and transcript["zero_count"] == 0
    probes = result["recorded_numeric_fields"]["control_probe_counts"]
    assert probes["valid_count"] == 2 and probes["zero_count"] == 1 and probes["absent_count"] == 1
    assert result["recorded_numeric_fields"]["genomic_control_counts"]["mean"] == 2
    assert result["observation_count"] == 3 and result["excluded_count"] == 1
    assert result["expression"]["library_counts"]["mean"] == 2
    assert result["expression"]["detected_features"]["mean"] == 1
    assert result["expression"]["zero_library_observation_count"] == 1
    assert result["expression"]["zero_library_examples"] == [
        {"cell_id": "b", "included": False, "label": "A", "imported_panel_count": 0.0}]
    assert result["segmentation_method"]["valid_count"] == 2
    assert result["segmentation_method"]["absent_count"] == 1
    assert result["scientific_authorization"] == "NOT_ESTABLISHED"
    assert not result["changes_applied"] and not result["selection_changed"]
    assert project.context()["cursor"] == before["cursor"]
    json.dumps(result, allow_nan=False)


def test_exact_selection_historical_revision_and_empty_selection(project):
    base = project.summary()["head_revision"]
    selection = project.set_selection(base, cell_ids=["a"])
    preview = project.propose_revision(base, selection["selection_id"], {"included": False}, "fixture")
    target = project.apply_revision(preview["proposal_id"], base, "Test driver", True)["revision_id"]
    qc = inspect_quality(project, base, selection["selection_id"])
    assert qc["historical"] and qc["observation_count"] == 1 and qc["included_count"] == 1
    assert qc["expression"]["library_counts"]["mean"] == 5
    with pytest.raises(SpatialError, match="exact inspected revision"):
        inspect_quality(project, target, selection["selection_id"])
    empty = project.set_selection(target, cell_ids=[])
    qc = inspect_quality(project, target, empty["selection_id"])
    assert qc["observation_count"] == 0 and qc["recorded_numeric_fields"]["cell_area"]["status"] == "no_observations"
    assert qc["expression"]["library_counts"]["mean"] is None


def test_source_bad_types_and_negative_qc_not_coerced_to_zero(tmp_path, metadata):
    values = [False, "0", -1, None, 0]
    observations = [{"cell_id": str(i), "x": i, "y": 0, "label": "unknown",
                     "attributes": {"nucleus_count": value}} for i, value in enumerate(values)]
    project = Project.create(tmp_path / "types", observations, metadata)
    metric = inspect_quality(project, project.summary()["head_revision"])["recorded_numeric_fields"]["nucleus_count"]
    assert metric["valid_count"] == 1 and metric["zero_count"] == 1
    assert metric["invalid_type_count"] == 2 and metric["invalid_range_count"] == 1 and metric["null_count"] == 1


def test_no_imported_features_is_unmeasured_not_zero_library(tmp_path, metadata):
    metadata["panel_genes"] = []
    project = Project.create(tmp_path / "empty_panel", [{"cell_id": "a", "x": 0, "y": 0, "label": "A"}], metadata)
    expression = inspect_quality(project, project.summary()["head_revision"])["expression"]
    assert expression["zero_library_observation_count"] is None
    assert expression["zero_library_examples"] == []
    assert expression["library_counts"]["status"] == "no_imported_features"
    assert expression["detected_features"]["mean"] is None


def test_zero_examples_and_category_output_bounded(tmp_path, metadata):
    observations = [{"cell_id": f"{i:03d}", "x": i, "y": 0, "label": "unknown",
                     "attributes": {"segmentation_method": f"source category {i}"}} for i in range(105)]
    project = Project.create(tmp_path / "bounded", observations, metadata)
    result = inspect_quality(project, project.summary()["head_revision"])
    assert result["expression"]["zero_library_observation_count"] == 105
    assert len(result["expression"]["zero_library_examples"]) == 100
    assert result["expression"]["zero_library_examples_omitted"] == 5
    assert result["segmentation_method"]["category_count"] == 105
    assert len(result["segmentation_method"]["categories"]) == 20
    assert result["segmentation_method"]["omitted_observation_count"] == 85


def test_finite_inputs_with_sum_overflow_remain_json_serializable(tmp_path, metadata, cells):
    cells = deepcopy(cells[:1])
    cells[0]["counts"] = {"G": 1e308, "H": 1e308}
    cells[0]["attributes"]["total_counts"] = 1e308
    project = Project.create(tmp_path / "overflow", cells, metadata)
    result = inspect_quality(project, project.summary()["head_revision"])
    assert result["expression"]["library_counts"]["overflow_count"] == 1
    assert result["expression"]["library_counts"]["mean"] is None
    assert math.isfinite(result["recorded_numeric_fields"]["total_counts"]["mean"])
    json.dumps(result, allow_nan=False)
