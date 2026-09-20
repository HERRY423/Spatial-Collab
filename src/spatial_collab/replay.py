"""Offline bundle integrity and independent recomputation from exported inputs.

A self-contained manifest establishes internal consistency, not provenance from
a trusted person or biological validity. No project database is opened or made.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any

from .store import SpatialError, _contains, _hash, _polygon, _text, _validate_import


MAX_BUNDLE_BYTES = 512 * 1024 * 1024
MAX_FILE_BYTES = 256 * 1024 * 1024
_BASE_FILES = {"source_snapshot.json", "revisions.json", "selection.json", "provenance.json",
               "proposals.json", "selections.json"}
_RUN_FILES = {"result.json", "result_status_at_export.json", "run_selection.json", "base_cells.json", "target_cells.json"}
_KNOWN_FILES = _BASE_FILES | _RUN_FILES | {"current_cells.json", "protein_assays.json"}


def _decode(content: bytes, name: str) -> Any:
    def pairs(entries):
        result = {}
        for key, value in entries:
            if key in result:
                raise SpatialError(f"Duplicate JSON key in {name}: {key}")
            result[key] = value
        return result
    def invalid_constant(value):
        raise SpatialError(f"Nonfinite JSON number in {name}: {value}")
    try:
        return json.loads(content.decode("utf-8"), object_pairs_hook=pairs, parse_constant=invalid_constant)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SpatialError(f"Invalid UTF-8 JSON in {name}.") from exc


def _safe_file(folder: Path, name: str) -> Path:
    if not isinstance(name, str) or re.fullmatch(r"[A-Za-z0-9_-]+\.json", name) is None:
        raise SpatialError("Manifest file names must be plain JSON basenames without paths.")
    candidate = folder / name
    if candidate.is_symlink() or getattr(candidate, "is_junction", lambda: False)():
        raise SpatialError(f"Bundle members must not be symbolic links or junctions: {name}")
    if candidate.resolve().parent != folder or not candidate.is_file():
        raise SpatialError(f"Missing or unsafe bundle member: {name}")
    return candidate


def _read_verified(folder: Path) -> tuple[dict, dict]:
    manifest_path = _safe_file(folder, "manifest.json")
    if manifest_path.stat().st_size > 64 * 1024:
        raise SpatialError("Bundle manifest exceeds the supported size.")
    manifest = _decode(manifest_path.read_bytes(), "manifest.json")
    if not isinstance(manifest, dict) or manifest.get("algorithm") != "sha256" or not isinstance(manifest.get("files"), dict):
        raise SpatialError("Bundle requires a SHA256 file manifest.")
    members = manifest["files"]
    if len(members) > len(_KNOWN_FILES):
        raise SpatialError("Bundle contains too many manifest members.")
    for name in members:
        _safe_file(folder, name)
    if set(members) - _KNOWN_FILES:
        raise SpatialError("Bundle contains unsupported manifest members.")
    has_run = "result.json" in members
    compact = manifest.get("layout") == "source_plus_overlays"
    if "layout" in manifest and not compact:
        raise SpatialError("Unsupported bundle layout.")
    result_files = _RUN_FILES - {"base_cells.json", "target_cells.json"} if compact else _RUN_FILES
    required = _BASE_FILES | (result_files if has_run else set() if compact else {"current_cells.json"})
    if "protein_assays.json" in members:
        required = required | {"protein_assays.json"}
    if set(members) != required:
        raise SpatialError(f"Bundle manifest member set is incomplete or inconsistent; required: {', '.join(sorted(required))}")
    total, files = 0, {}
    for name, entry in members.items():
        if not isinstance(entry, dict):
            raise SpatialError(f"Invalid manifest entry for {name}.")
        size, digest = entry.get("bytes"), entry.get("sha256")
        if type(size) is not int or size < 0 or size > MAX_FILE_BYTES:
            raise SpatialError(f"Unsupported file size for {name}.")
        total += size
        if total > MAX_BUNDLE_BYTES:
            raise SpatialError("Bundle exceeds the 512 MiB verification resource limit.")
        if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
            raise SpatialError(f"Invalid SHA256 for {name}.")
        path = _safe_file(folder, name)
        if path.stat().st_size != size:
            raise SpatialError(f"File size mismatch: {name}")
        content = path.read_bytes()
        if len(content) != size or hashlib.sha256(content).hexdigest() != digest:
            raise SpatialError(f"SHA256 mismatch: {name}")
        files[name] = _decode(content, name)
    if (files["provenance.json"].get("materialization") == "source_plus_overlays") != compact:
        raise SpatialError("Manifest and provenance materialization layouts differ.")
    return files, manifest


class _BundleProject:
    """Read-only adapter supplies original IDs without creating another project."""

    def __init__(self, source, revisions, selections, provenance):
        self.source = source
        self.revisions = revisions
        self.selections = selections
        self.provenance = provenance

    def summary(self):
        return {"metadata": deepcopy(self.source["metadata"]), "cell_count": len(self.source["cells"]),
                "project_id": self.provenance["project_id"], "source_sha256": _hash(self.source),
                "head_revision": self.provenance["head_revision_at_export"]}

    def get_revision(self, revision_id):
        try:
            revision = deepcopy(self.revisions[revision_id])
        except KeyError:
            raise SpatialError(f"Unknown exported revision: {revision_id}") from None
        return {**revision, "revision_sha256": _hash(revision)}

    def cells(self, revision_id):
        overlays = self.get_revision(revision_id)["overlays"]
        result = deepcopy(self.source["cells"])
        for cell in result:
            cell.update(overlays.get(cell["cell_id"], {}))
        return result

    def get_selection(self, selection_id=None):
        if selection_id is None:
            return None
        try:
            return deepcopy(self.selections[selection_id])
        except KeyError:
            raise SpatialError(f"Unknown exported selection: {selection_id}") from None

    def save_run(self, result):
        return result


def _validate_records(files: dict) -> tuple[_BundleProject, dict | None]:
    source, provenance = files["source_snapshot.json"], files["provenance.json"]
    if not isinstance(source, dict) or set(source) != {"cells", "metadata"}:
        raise SpatialError("Invalid source snapshot.")
    if _validate_import(source["cells"], source["metadata"]) != source:
        raise SpatialError("Source snapshot is not in canonical imported form.")
    if not isinstance(provenance, dict) or provenance.get("format") != "spatial-collab-review-bundle-1":
        raise SpatialError("Unsupported review bundle format.")
    _text(provenance.get("project_id"), "project_id")
    digest = _hash(source)
    if provenance.get("source_sha256") != digest:
        raise SpatialError("Canonical source hash does not match export provenance.")
    revisions_list = files["revisions.json"]
    if not isinstance(revisions_list, list) or not revisions_list:
        raise SpatialError("Bundle requires revision history.")
    cells_by_id = {cell["cell_id"]: cell for cell in source["cells"]}
    revisions, previous = {}, None
    for revision in revisions_list:
        if not isinstance(revision, dict):
            raise SpatialError("Each revision must be an object.")
        rid = _text(revision.get("revision_id"), "revision_id")
        if rid in revisions or revision.get("parent_revision") != previous:
            raise SpatialError("Revision history must be a unique, ordered, append-only chain.")
        if revision.get("source_sha256") != digest:
            raise SpatialError("Revision references a different source snapshot.")
        if revision.get("scientific_authorization") != "NOT_ESTABLISHED":
            raise SpatialError("Bundle exceeds the supported scientific authorization boundary.")
        overlays = revision.get("overlays")
        if not isinstance(overlays, dict) or set(overlays) - set(cells_by_id):
            raise SpatialError("Revision overlays contain unknown cells.")
        for overlay in overlays.values():
            if not isinstance(overlay, dict) or set(overlay) - {"label", "included", "region"}:
                raise SpatialError("Revision attempts to modify immutable source coordinates or counts.")
            if "included" in overlay and type(overlay["included"]) is not bool:
                raise SpatialError("Revision inclusion overlays must be boolean.")
            for field in ("label", "region"):
                if field in overlay:
                    _text(overlay[field], field, empty=field == "region", limit=512)
        if previous is None:
            if revision.get("kind") != "import" or overlays:
                raise SpatialError("First revision must represent the immutable import.")
        else:
            if revision.get("kind") not in {"edit", "revert"} or revision.get("confirmation") is not True:
                raise SpatialError("Revision requires a recorded edit/revert confirmation.")
            _text(revision.get("reviewer"), "reviewer")
            if revision["kind"] == "revert":
                target = revision.get("reverted_to_revision")
                if target not in revisions or overlays != revisions[target]["overlays"]:
                    raise SpatialError("Revert overlay does not match its recorded target.")
        revisions[rid] = revision
        previous = rid
    if previous != provenance.get("head_revision_at_export"):
        raise SpatialError("Exported head does not match the revision history.")
    expected_hashes = {rid: _hash(revision) for rid, revision in revisions.items()}
    if provenance.get("revision_hashes") != expected_hashes:
        raise SpatialError("Canonical revision hashes do not match export provenance.")
    selections_list = files["selections.json"]
    if not isinstance(selections_list, list):
        raise SpatialError("Exported selections must be a list.")
    selections = {}
    for selection in selections_list:
        if not isinstance(selection, dict):
            raise SpatialError("Each exported selection must be an object.")
        sid = _text(selection.get("selection_id"), "selection_id")
        if sid in selections or selection.get("revision_id") not in revisions:
            raise SpatialError("Selection IDs must be unique and refer to an exported revision.")
        for key in ("slice_id", "coordinate_system", "units"):
            if selection.get(key) != source["metadata"].get(key):
                raise SpatialError(f"Selection {key} does not match the source snapshot.")
        if selection.get("source_sha256") != digest:
            raise SpatialError("Selection references a different source snapshot.")
        ids = selection.get("cell_ids")
        if (not isinstance(ids, list) or any(not isinstance(cid, str) for cid in ids)
                or len(set(ids)) != len(ids) or set(ids) - set(cells_by_id)
                or selection.get("cell_count") != len(ids)):
            raise SpatialError("Selection must contain an exact, unique set of known cell IDs.")
        if selection.get("polygon") is not None:
            polygon = _polygon(selection["polygon"])
            contained = {cid for cid, cell in cells_by_id.items() if _contains(cell["x"], cell["y"], polygon)}
            if contained != set(ids):
                raise SpatialError("Exported polygon does not select its recorded centroid IDs.")
        selections[sid] = selection
    def validate_selected(record):
        if record is None:
            return
        if not isinstance(record, dict):
            raise SpatialError("Invalid selected-region record.")
        _text(record.get("selection_id"), "selection_id")
        canonical = {key: value for key, value in record.items() if key != "stale"}
        if selections.get(record.get("selection_id")) != canonical:
            raise SpatialError("Selected-region record differs from immutable selection history.")
        if type(record.get("stale")) is not bool or record["stale"] != (record.get("revision_id") != previous):
            raise SpatialError("Selected-region staleness is inconsistent with the exported head.")
    validate_selected(files["selection.json"])
    project = _BundleProject(source, revisions, selections, provenance)
    project.protein_assays = files.get("protein_assays.json", [])
    if not isinstance(project.protein_assays, list) or len(project.protein_assays) > 16:
        raise SpatialError("Invalid exported protein assays.")
    from .proteomics import _records
    _records(project)
    proposals_list = files["proposals.json"]
    if not isinstance(proposals_list, list):
        raise SpatialError("Exported proposals must be a list.")
    proposals = {}
    for proposal in proposals_list:
        if not isinstance(proposal, dict):
            raise SpatialError("Each exported proposal must be an object.")
        pid = _text(proposal.get("proposal_id"), "proposal_id")
        base = _text(proposal.get("base_revision"), "base_revision")
        sid = _text(proposal.get("selection_id"), "selection_id")
        if pid in proposals or base not in revisions or sid not in selections or selections[sid]["revision_id"] != base:
            raise SpatialError("Proposal must reference its exact revision and frozen selection.")
        changes = proposal.get("changes")
        if not isinstance(changes, dict) or not changes or set(changes) - {"label", "included", "region"}:
            raise SpatialError("Proposal must contain supported overlay changes.")
        if "included" in changes and type(changes["included"]) is not bool:
            raise SpatialError("Proposal inclusion must be boolean.")
        for field in ("label", "region"):
            if field in changes:
                _text(changes[field], field, empty=field == "region", limit=512)
        chosen = set(selections[sid]["cell_ids"])
        deltas = []
        for cell in project.cells(base):
            if cell["cell_id"] in chosen:
                fields = {key: {"before": cell[key], "after": value} for key, value in changes.items() if cell[key] != value}
                if fields:
                    deltas.append({"cell_id": cell["cell_id"], "changes": fields})
        if (not deltas or proposal.get("deltas") != deltas or proposal.get("selected_cell_count") != len(chosen)
                or proposal.get("changed_cell_count") != len(deltas)):
            raise SpatialError("Proposal preview does not match its exact input cells and changes.")
        proposals[pid] = proposal
    for revision in revisions.values():
        if revision["kind"] != "edit":
            continue
        pid = _text(revision.get("proposal_id"), "revision proposal_id")
        proposal = proposals.get(pid)
        if proposal is None or proposal["base_revision"] != revision["parent_revision"]:
            raise SpatialError("Committed edit does not refer to its reviewed proposal.")
        expected_overlays = deepcopy(revisions[revision["parent_revision"]]["overlays"])
        for delta in proposal["deltas"]:
            expected_overlays.setdefault(delta["cell_id"], {}).update({key: change["after"] for key, change in delta["changes"].items()})
        if expected_overlays != revision["overlays"]:
            raise SpatialError("Committed overlays do not match the reviewed proposal.")
    run = files.get("result.json")
    compact = provenance.get("materialization") == "source_plus_overlays"
    if run is None:
        if not compact and files["current_cells.json"] != project.cells(previous):
            raise SpatialError("Exported current cells do not match source plus revision overlays.")
        return project, None
    if not isinstance(run, dict):
        raise SpatialError("Exported result must be an object.")
    if provenance.get("run_sha256") != _hash(run):
        raise SpatialError("Canonical result hash does not match export provenance.")
    base = _text(run.get("base_revision"), "result base_revision")
    target = _text(run.get("target_revision"), "result target_revision")
    if base not in revisions or target not in revisions:
        raise SpatialError("Result refers to missing revisions.")
    expected = {"source_sha256": digest, "base_revision_sha256": expected_hashes[base], "target_revision_sha256": expected_hashes[target]}
    if run.get("data_hashes") != expected:
        raise SpatialError("Result data hashes do not match the reconstructed inputs.")
    run_provenance = run.get("provenance")
    if not isinstance(run_provenance, dict) or any(run_provenance.get(k) != v for k, v in expected.items()):
        raise SpatialError("Result provenance does not match the reconstructed inputs.")
    if run_provenance.get("project_id") != provenance["project_id"]:
        raise SpatialError("Result refers to a different project.")
    for name, revision_id in (("base_cells.json", base), ("target_cells.json", target)):
        if not compact and files[name] != project.cells(revision_id):
            raise SpatialError(f"{name} does not match immutable source plus revision overlays.")
    status = files["result_status_at_export.json"]
    if status != {"run_id": run.get("run_id"), "stale": target != previous, "head_revision_at_export": previous}:
        raise SpatialError("Exported result staleness is inconsistent with the exported head.")
    selected = files["run_selection.json"]
    validate_selected(selected)
    if (selected or {}).get("selection_id") != run.get("selection_id"):
        raise SpatialError("Result and exported selection IDs do not match.")
    if run.get("analysis_schema") == "spatial-collab.hypothesis-plan-run.v2":
        params = run.get("parameters")
        if not isinstance(params, dict) or set(params) != {"frozen_plan"} or base != target:
            raise SpatialError("Hypothesis plan requires a frozen plan and one real revision anchor.")
        frozen = params["frozen_plan"]
        spec = frozen.get("spec") if isinstance(frozen, dict) else None
        if not isinstance(spec, dict) or (run["selection_id"] != spec.get("selection_id")
                                         or base != spec.get("revision_id")):
            raise SpatialError("Result revision and selection anchors must match the embedded frozen plan.")
    if run.get("analysis_schema") == "spatial-collab.protein-region.v1":
        from .proteomics import get_assay
        assay = get_assay(project, run.get("parameters", {}).get("assay_id"))
        if run.get("assay_sha256") != assay["assay_sha256"]:
            raise SpatialError("Protein result does not match its immutable assay.")
    return project, run


def _equal(left: Any, right: Any) -> bool:
    if isinstance(left, bool) or isinstance(right, bool):
        return type(left) is type(right) and left == right
    if isinstance(left, int) or isinstance(right, int):
        return type(left) is type(right) and left == right
    if isinstance(left, float) and isinstance(right, float):
        return math.isfinite(left) and math.isfinite(right) and math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-12)
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(_equal(left[k], right[k]) for k in left)
    if isinstance(left, list):
        return len(left) == len(right) and all(_equal(a, b) for a, b in zip(left, right))
    return left == right


def verify_bundle(bundle_path: str | Path, recompute: bool = True) -> dict:
    """Verify a review bundle, then recompute its metric from source + overlays.

    Invalid/tampered bundles raise SpatialError. A valid bundle with a differing
    recomputation returns ``recompute_status='mismatch'`` and the mismatched
    top-level scientific fields. The original bundle is never modified.
    """
    if type(recompute) is not bool:
        raise SpatialError("recompute must be true or false.")
    folder = Path(bundle_path).resolve()
    if not folder.is_dir():
        raise SpatialError("Bundle path must be an existing exported review directory.")
    files, manifest = _read_verified(folder)
    project, run = _validate_records(files)
    result = {"bundle_path": str(folder), "integrity": "verified", "integrity_verified": True,
              "files_verified": len(manifest["files"]), "source_sha256": _hash(project.source),
              "revision_count": len(project.revisions), "run_id": run.get("run_id") if run else None,
              "recompute_status": "not_requested" if not recompute else "no_result_in_bundle",
              "recompute_matches": None, "mismatched_fields": [],
              "scientific_authorization": "NOT_ESTABLISHED",
              "limitations": ["Manifest verification establishes internal consistency, not external authenticity.",
                              "Recomputation verifies the declared descriptive calculation, not biological validity or the original scientific conclusion."]}
    if not recompute or run is None:
        return result
    params = run.get("parameters")
    if not isinstance(params, dict):
        raise SpatialError("Result parameters must be an object.")
    schema = run.get("analysis_schema")
    if schema == "spatial-collab.descriptive-sensitivity.v1":
        from .analysis import compare
        names = ("radius_um", "graph_scope", "source_label", "target_label", "min_effect")
    elif schema == "spatial-collab.region-contrast.v1":
        from .exploration import compare_regions as compare
        names = ("background_selection_id", "genes")
    elif schema == "spatial-collab.hypothesis-sensitivity.v1":
        from .sensitivity import run_sensitivity
        names = ("variants", "radii_um", "source_label", "target_label", "graph_scope", "min_effect")
        if any(name not in params for name in names) or run["base_revision"] != run["target_revision"]:
            raise SpatialError("Hypothesis result requires one real revision anchor and complete parameters.")
        recomputed = run_sensitivity(project, run["base_revision"], run["selection_id"], **{name: params[name] for name in names})
        fields = ("results", "parameters", "variants_sha256", "data_hashes", "selection", "semantics", "import_scope",
                  "revision_created", "researcher_approval", "evidence_ceiling", "scientific_authorization")
        mismatches = [field for field in fields if not _equal(run.get(field), recomputed.get(field))]
        result.update(recompute_status="mismatch" if mismatches else "matched", recompute_matches=not mismatches,
                      mismatched_fields=mismatches, compared_fields=list(fields), numeric_tolerance={"relative": 1e-12, "absolute": 1e-12})
        return result
    elif schema == "spatial-collab.hypothesis-plan-run.v2":
        from .hypotheses import run_hypothesis_plan
        recomputed = run_hypothesis_plan(project, frozen_plan=params["frozen_plan"])
        fields = ("base_revision", "target_revision", "selection_id", "parameters", "plan_sha256", "results", "difference_order", "difference_order_rule",
                  "status_counts", "data_hashes", "selection", "semantics", "import_scope", "methods",
                  "revision_created", "researcher_approval", "evidence_ceiling", "scientific_authorization")
        mismatches = [field for field in fields if not _equal(run.get(field), recomputed.get(field))]
        result.update(recompute_status="mismatch" if mismatches else "matched", recompute_matches=not mismatches,
                      mismatched_fields=mismatches, compared_fields=list(fields), numeric_tolerance={"relative": 1e-12, "absolute": 1e-12})
        return result
    elif schema == "spatial-collab.protein-region.v1":
        from .proteomics import compare_protein_regions
        names = ("assay_id", "features", "background_selection_id", "scale", "cofactor", "rna_gene")
        if set(params) != set(names):
            raise SpatialError("Protein result requires complete declared parameters.")
        recomputed = compare_protein_regions(project, base_revision=run["base_revision"], target_revision=run["target_revision"],
                                             selection_id=run["selection_id"], **params)
        fields = ("before", "after", "comparison", "parameters", "assay_sha256", "data_hashes", "selection", "evidence_ceiling", "scientific_authorization")
        mismatches = [field for field in fields if not _equal(run.get(field), recomputed.get(field))]
        result.update(recompute_status="mismatch" if mismatches else "matched", recompute_matches=not mismatches,
                      mismatched_fields=mismatches, compared_fields=list(fields))
        return result
    else:
        raise SpatialError("Unsupported analysis schema for recomputation.")
    if any(name not in params for name in names):
        raise SpatialError("Result lacks required recomputation parameters.")
    recomputed = compare(project, run["base_revision"], run["target_revision"],
                         selection_id=run.get("selection_id"), **{name: params[name] for name in names})
    fields = ("before", "after", "comparison", "parameters", "data_hashes", "selection", "evidence_ceiling", "scientific_authorization")
    if "semantics" in run:
        fields += ("semantics",)
    mismatches = [field for field in fields if not _equal(run.get(field), recomputed.get(field))]
    result.update(recompute_status="mismatch" if mismatches else "matched", recompute_matches=not mismatches,
                  mismatched_fields=mismatches, compared_fields=list(fields), numeric_tolerance={"relative": 1e-12, "absolute": 1e-12})
    return result
