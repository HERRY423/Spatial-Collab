import hashlib
import json
from pathlib import Path

import pytest

from spatial_collab.analysis import compare
from spatial_collab.replay import verify_bundle
from spatial_collab.store import Project, SpatialError, _hash


@pytest.fixture
def project(tmp_path):
    cells = [{"cell_id": "t", "x": 0, "y": 0, "label": "T cell", "counts": {"CD3D": 4}},
             {"cell_id": "m", "x": 1, "y": 0, "label": "Myeloid", "counts": {}},
             {"cell_id": "b", "x": 100, "y": 0, "label": "B cell", "counts": {}}]
    metadata = {"name": "replay fixture", "slice_id": "s", "coordinate_system": "s_um", "units": "micrometer",
                "panel_genes": ["CD3D"], "source_kind": "synthetic", "source_files": [],
                "biological_replicates": 0, "limitations": ["Synthetic fixture only"]}
    return Project.create(tmp_path / "project", cells, metadata)


@pytest.fixture
def bundle(project):
    base = project.summary()["head_revision"]
    selection = project.set_selection(base, cell_ids=["t", "m", "b"])
    preview = project.propose_revision(base, selection["selection_id"], {"region": "reviewed"}, "Region review")
    target = project.apply_revision(preview["proposal_id"], base, "Researcher", True)["revision_id"]
    run = compare(project, base, target, selection_id=selection["selection_id"], radius_um=2)
    return Path(project.export_bundle(run["run_id"])["export_path"])


def rewrite(folder, filename, data):
    path = folder / filename
    content = (json.dumps(data, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode()
    path.write_bytes(content)
    manifest = json.loads((folder / "manifest.json").read_text())
    manifest["files"][filename] = {"sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)}
    (folder / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")


def test_replay_reconstructs_and_recomputes_without_writing(bundle):
    before = {path.name: path.read_bytes() for path in bundle.iterdir()}
    result = verify_bundle(bundle)
    assert result["integrity"] == "verified"
    assert result["recompute_status"] == "matched"
    assert result["recompute_matches"]
    assert result["scientific_authorization"] == "NOT_ESTABLISHED"
    assert {path.name: path.read_bytes() for path in bundle.iterdir()} == before
    assert not (bundle / "project.sqlite3").exists()


def test_export_without_run_has_integrity_and_no_computation(project):
    result = verify_bundle(project.export_bundle()["export_path"])
    assert result["integrity_verified"]
    assert result["recompute_status"] == "no_result_in_bundle"


def test_integrity_only_is_explicit(bundle):
    result = verify_bundle(bundle, recompute=False)
    assert result["recompute_status"] == "not_requested"
    assert result["recompute_matches"] is None


def test_modified_bytes_are_detected(bundle):
    with (bundle / "result.json").open("ab") as stream:
        stream.write(b" ")
    with pytest.raises(SpatialError, match="size mismatch"):
        verify_bundle(bundle)


def test_missing_file_detected(bundle):
    (bundle / "base_cells.json").unlink()
    with pytest.raises(SpatialError, match="Missing"):
        verify_bundle(bundle)


@pytest.mark.parametrize("name", ["../outside.json", "..\\outside.json", "C:\\outside.json", "/outside.json"])
def test_manifest_cannot_read_arbitrary_files(bundle, name):
    manifest = json.loads((bundle / "manifest.json").read_text())
    manifest["files"].pop("base_cells.json")
    manifest["files"][name] = {"sha256": "0" * 64, "bytes": 0}
    (bundle / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(SpatialError, match="basenames"):
        verify_bundle(bundle)


def test_rehashed_source_tampering_fails_canonical_reference(bundle):
    source = json.loads((bundle / "source_snapshot.json").read_text())
    source["cells"][0]["x"] += 1
    rewrite(bundle, "source_snapshot.json", source)
    with pytest.raises(SpatialError, match="Canonical source hash"):
        verify_bundle(bundle)


def test_exported_cells_are_checked_against_source_and_overlays(bundle):
    target = json.loads((bundle / "target_cells.json").read_text())
    target[0]["included"] = False
    rewrite(bundle, "target_cells.json", target)
    with pytest.raises(SpatialError, match="source plus revision"):
        verify_bundle(bundle)


def test_rehashed_self_consistent_result_still_requires_recomputation(bundle):
    run = json.loads((bundle / "result.json").read_text())
    run["before"]["observed_fraction"] = 0.125
    rewrite(bundle, "result.json", run)
    provenance = json.loads((bundle / "provenance.json").read_text())
    provenance["run_sha256"] = _hash(run)
    rewrite(bundle, "provenance.json", provenance)
    result = verify_bundle(bundle)
    assert result["integrity_verified"]
    assert result["recompute_status"] == "mismatch"
    assert result["mismatched_fields"] == ["before"]


def test_revision_cannot_edit_source_coordinates(bundle):
    revisions = json.loads((bundle / "revisions.json").read_text())
    revisions[-1]["overlays"]["t"]["x"] = 123
    rewrite(bundle, "revisions.json", revisions)
    with pytest.raises(SpatialError, match="immutable source"):
        verify_bundle(bundle)


def test_no_duplicate_json_keys(bundle):
    content = '{"algorithm":"sha256","algorithm":"sha256","files":{}}'
    (bundle / "manifest.json").write_text(content)
    with pytest.raises(SpatialError, match="Duplicate JSON key"):
        verify_bundle(bundle)


def test_required_file_set_cannot_be_trimmed(bundle):
    manifest = json.loads((bundle / "manifest.json").read_text())
    manifest["files"].pop("run_selection.json")
    (bundle / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(SpatialError, match="member set"):
        verify_bundle(bundle)


def test_frozen_selection_has_exact_cell_universe(bundle):
    selections = json.loads((bundle / "selections.json").read_text())
    selections[0]["cell_ids"].append("missing")
    selections[0]["cell_count"] += 1
    rewrite(bundle, "selections.json", selections)
    with pytest.raises(SpatialError, match="known cell IDs"):
        verify_bundle(bundle)


def test_proposal_preview_cannot_disagree_with_committed_overlay(bundle):
    proposals = json.loads((bundle / "proposals.json").read_text())
    proposals[0]["deltas"][0]["changes"]["region"]["after"] = "different"
    rewrite(bundle, "proposals.json", proposals)
    with pytest.raises(SpatialError, match="Proposal preview"):
        verify_bundle(bundle)


def test_revert_bundle_replays_recorded_overlay_target(project):
    base = project.summary()["head_revision"]
    selection = project.set_selection(base, cell_ids=["t"])
    proposal = project.propose_revision(base, selection["selection_id"], {"label": "Other"}, "review")
    edited = project.apply_revision(proposal["proposal_id"], base, "Researcher", True)["revision_id"]
    reverted = project.revert_revision(base, edited, "Researcher", True)["revision_id"]
    run = compare(project, base, reverted, radius_um=2)
    assert verify_bundle(project.export_bundle(run["run_id"])["export_path"])["recompute_matches"]


def test_stale_result_replays_at_original_revision(project):
    base = project.summary()["head_revision"]
    run = compare(project, base, base, radius_um=2)
    selection = project.set_selection(base, cell_ids=["m"])
    proposal = project.propose_revision(base, selection["selection_id"], {"included": False}, "exclude uncertain cell")
    project.apply_revision(proposal["proposal_id"], base, "Researcher", True)
    bundle = project.export_bundle(run["run_id"])
    result = verify_bundle(bundle["export_path"])
    assert result["recompute_status"] == "matched"
