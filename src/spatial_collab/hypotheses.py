"""Versioned, frozen what-if plans and two explicitly different estimands.

Plan versions live in a sidecar ledger, never in the scientific revision chain.
Frozen parameters are embedded in runs so offline replay needs no ledger.
"""
from __future__ import annotations

from collections import Counter
from contextlib import contextmanager
from copy import deepcopy
import json
import math
import sqlite3

import numpy as np
from scipy.spatial import cKDTree

from .analysis import (
    MAX_DIRECTED_SOURCE_EDGES, MAX_NEIGHBORS_PER_SOURCE, MAX_PROJECT_CELLS,
    _composition, _number, _runtime_provenance,
)
from .exploration import semantics
from .store import SpatialError, _hash, _id, _json, _now, _text


PLAN_SCHEMA = "spatial-collab.hypothesis-plan.v2"
RUN_SCHEMA = "spatial-collab.hypothesis-plan-run.v2"
METHODS = ("edge_weighted", "source_equal_weighted")
MAX_VARIANTS = 8
MAX_GRID_ROWS = 96
_SPEC_KEYS = {"name", "rationale", "revision_id", "selection_id", "source_label", "target_label",
              "radii_um", "graph_scopes", "min_effect", "background", "variants"}
_OPS = {"eq", "ne", "gt", "ge", "lt", "le", "in"}


def _scalar(value):
    if value is None or type(value) in (str, bool):
        return True
    try:
        return type(value) in (int, float) and math.isfinite(value)
    except OverflowError:
        return False


def _validate_spec(spec):
    if not isinstance(spec, dict) or set(spec) != _SPEC_KEYS:
        raise SpatialError("Plan spec requires exactly: " + ", ".join(sorted(_SPEC_KEYS)))
    spec = deepcopy(spec)
    for key in ("name", "rationale", "revision_id", "selection_id", "source_label", "target_label"):
        _text(spec[key], key, limit=4000 if key == "rationale" else 512)
    if spec["source_label"] == "Unannotated" or spec["target_label"] == "Unannotated":
        raise SpatialError("Declare working groups; Unannotated is not an inferred biological group.")
    radii = spec["radii_um"]
    if not isinstance(radii, list) or not 1 <= len(radii) <= 6:
        raise SpatialError("Declare 1..6 radii.")
    spec["radii_um"] = [_number(r, "radius_um") for r in radii]
    if len(set(spec["radii_um"])) != len(radii):
        raise SpatialError("Declared radii must be distinct.")
    scopes = spec["graph_scopes"]
    if (not isinstance(scopes, list) or not scopes or len(scopes) > 2
            or any(not isinstance(s, str) or s not in {"whole_slice", "roi_induced"} for s in scopes)
            or len(set(scopes)) != len(scopes)):
        raise SpatialError("Declare distinct whole_slice and/or roi_induced graph scopes.")
    spec["min_effect"] = _number(spec["min_effect"], "min_effect", maximum=1.0)
    background = spec["background"]
    if not isinstance(background, dict) or not isinstance(background.get("kind"), str) or background.get("kind") not in {
        "neighbor_universe", "imported_universe", "selection"
    }:
        raise SpatialError("Declare background kind: neighbor_universe, imported_universe or selection.")
    if background["kind"] == "selection":
        if set(background) != {"kind", "selection_id"}:
            raise SpatialError("Selected background requires exactly kind and selection_id.")
        _text(background["selection_id"], "background selection_id")
    elif set(background) != {"kind"}:
        raise SpatialError("Background has unexpected fields.")
    variants = spec["variants"]
    if not isinstance(variants, list) or not 1 <= len(variants) <= MAX_VARIANTS:
        raise SpatialError(f"Declare 1..{MAX_VARIANTS} variants.")
    if len(variants) * len(radii) * len(scopes) > MAX_GRID_ROWS:
        raise SpatialError("Hypothesis grid exceeds the declared row budget.")
    names = set()
    for variant in variants:
        if not isinstance(variant, dict) or set(variant) != {"name", "rationale", "selector", "changes"}:
            raise SpatialError("Each variant requires name, rationale, selector and changes.")
        name = _text(variant["name"], "variant name", limit=128)
        if name in names or name == "baseline":
            raise SpatialError("Variant names must be distinct and must not be baseline.")
        names.add(name)
        _text(variant["rationale"], "variant rationale")
        changes = variant["changes"]
        if not isinstance(changes, dict) or not changes or set(changes) - {"label", "included"}:
            raise SpatialError("Hypotheses may change label and/or included only.")
        if "label" in changes:
            _text(changes["label"], "hypothetical label", limit=512)
        if "included" in changes and type(changes["included"]) is not bool:
            raise SpatialError("included must be a boolean.")
        selector = variant["selector"]
        if not isinstance(selector, dict) or set(selector) not in ({"cell_ids"}, {"predicate"}):
            raise SpatialError("Selector requires either exact cell_ids or predicate.")
        if "cell_ids" in selector:
            ids = selector["cell_ids"]
            if (not isinstance(ids, list) or not ids or len(ids) > MAX_PROJECT_CELLS
                    or any(not isinstance(x, str) or not x for x in ids) or len(set(ids)) != len(ids)):
                raise SpatialError("Explicit IDs must be distinct, nonempty and bounded.")
            selector["cell_ids"] = sorted(ids)
        else:
            predicate = selector["predicate"]
            if not isinstance(predicate, dict) or set(predicate) != {"all"}:
                raise SpatialError("A predicate is an explicit all-list of clauses, never executable code.")
            clauses = predicate["all"]
            if not isinstance(clauses, list) or not 1 <= len(clauses) <= 8:
                raise SpatialError("Predicates require 1..8 clauses.")
            for clause in clauses:
                if not isinstance(clause, dict) or set(clause) != {"field", "op", "value"}:
                    raise SpatialError("Predicate clauses require field, op and value.")
                field = _text(clause["field"], "predicate field", limit=520)
                if field not in {"label", "included", "region"} and not (
                    field.startswith("attributes.") and len(field) > 11
                    or field.startswith("counts.") and len(field) > 7
                ):
                    raise SpatialError("Predicate field must be label/included/region, attributes.KEY or counts.EXACT_FEATURE_ID.")
                if not isinstance(clause["op"], str) or clause["op"] not in _OPS:
                    raise SpatialError("Unsupported predicate operator.")
                value = clause["value"]
                if clause["op"] == "in":
                    if not isinstance(value, list) or not 1 <= len(value) <= 32 or not all(_scalar(x) for x in value):
                        raise SpatialError("in requires 1..32 finite scalar values.")
                elif not _scalar(value):
                    raise SpatialError("Predicate value must be a finite JSON scalar.")
                if clause["op"] in {"gt", "ge", "lt", "le"} and type(value) not in (int, float):
                    raise SpatialError("Ordered comparisons require a numeric value.")
    _json(spec)
    return spec


@contextmanager
def _ledger(project, *, create=False):
    path = project.root / "hypotheses.sqlite3"
    if path.is_symlink() or getattr(path, "is_junction", lambda: False)():
        raise SpatialError("Hypothesis ledger must be a regular file within the project.")
    if not create and not path.exists():
        raise SpatialError("No hypothesis plans have been created in this project.")
    summary = project.summary()
    mode = "rwc" if create else "rw"
    db = sqlite3.connect(path.as_uri() + "?mode=" + mode, uri=True, timeout=30, isolation_level=None)
    try:
        if create:
            db.executescript("""
                CREATE TABLE IF NOT EXISTS identity(key TEXT PRIMARY KEY,value TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS plans(plan_id TEXT PRIMARY KEY,current_version INTEGER NOT NULL);
                CREATE TABLE IF NOT EXISTS versions(plan_id TEXT,version INTEGER,payload TEXT NOT NULL,
                    sha256 TEXT NOT NULL,PRIMARY KEY(plan_id,version));
                CREATE TRIGGER IF NOT EXISTS versions_no_update BEFORE UPDATE ON versions
                    BEGIN SELECT RAISE(ABORT,'Immutable plan version'); END;
                CREATE TRIGGER IF NOT EXISTS versions_no_delete BEFORE DELETE ON versions
                    BEGIN SELECT RAISE(ABORT,'Immutable plan version'); END;
            """)
            db.execute("BEGIN IMMEDIATE")
            db.executemany("INSERT OR IGNORE INTO identity VALUES(?,?)", [
                ("source_sha256", summary["source_sha256"]), ("project_id", summary["project_id"])])
            db.commit()
        identity = dict(db.execute("SELECT key,value FROM identity"))
        if any(identity.get(k) != summary[k] for k in ("source_sha256", "project_id")):
            raise SpatialError("Hypothesis ledger belongs to a different source or project.")
        yield db
    except sqlite3.Error as exc:
        if db.in_transaction:
            db.rollback()
        raise SpatialError("Hypothesis ledger transaction failed: " + str(exc)) from exc
    finally:
        db.close()


def _read_version(db, plan_id, version=None):
    _text(plan_id, "plan_id")
    if version is None:
        row = db.execute("SELECT current_version FROM plans WHERE plan_id=?", (plan_id,)).fetchone()
        if row is None:
            raise SpatialError("Unknown hypothesis plan.")
        version = row[0]
    if type(version) is not int or version < 1:
        raise SpatialError("Plan version must be a positive integer.")
    row = db.execute("SELECT payload,sha256 FROM versions WHERE plan_id=? AND version=?", (plan_id, version)).fetchone()
    if row is None:
        raise SpatialError("Unknown hypothesis plan version.")
    try:
        record = json.loads(row[0])
    except (TypeError, json.JSONDecodeError) as exc:
        raise SpatialError("Plan version payload is not valid JSON.") from exc
    if _hash(record) != row[1] or record.get("plan_id") != plan_id or record.get("version") != version:
        raise SpatialError("Plan version failed its integrity check.")
    return {**record, "plan_sha256": row[1]}


def _append_version(db, plan_id, expected_version, spec, status, resolved=None):
    db.execute("BEGIN IMMEDIATE")
    try:
        if expected_version is None:
            version, parent_hash = 1, None
            db.execute("INSERT INTO plans VALUES(?,?)", (plan_id, version))
        else:
            if type(expected_version) is not int or expected_version < 1:
                raise SpatialError("expected_version must be a positive integer.")
            current = _read_version(db, plan_id)
            if current["version"] != expected_version:
                raise SpatialError("Stale hypothesis plan version; refresh before editing or freezing.")
            version, parent_hash = expected_version + 1, current["plan_sha256"]
            db.execute("UPDATE plans SET current_version=? WHERE plan_id=? AND current_version=?",
                       (version, plan_id, expected_version))
        record = {"schema": PLAN_SCHEMA, "plan_id": plan_id, "version": version, "status": status,
                  "parent_version": expected_version, "parent_sha256": parent_hash, "created_at": _now(),
                  "spec": spec, "resolved": resolved, "scientific_authorization": "NOT_ESTABLISHED",
                  "freeze_semantics": "Explicit computational plan freeze; not human scientific approval."}
        checksum = _hash(record)
        db.execute("INSERT INTO versions VALUES(?,?,?,?)", (plan_id, version, _json(record), checksum))
        db.commit()
        return {**record, "plan_sha256": checksum}
    except Exception:
        db.rollback()
        raise


def create_hypothesis_plan(project, spec):
    spec = _validate_spec(spec)
    with _ledger(project, create=True) as db:
        return _append_version(db, _id("plan"), None, spec, "draft")


def revise_hypothesis_plan(project, plan_id, expected_version, spec):
    """Create a new draft version; an older frozen version remains executable."""
    spec = _validate_spec(spec)
    with _ledger(project) as db:
        return _append_version(db, plan_id, expected_version, spec, "draft")


def get_hypothesis_plan(project, plan_id, version=None):
    with _ledger(project) as db:
        return _read_version(db, plan_id, version)


def list_hypothesis_plans(project):
    if not (project.root / "hypotheses.sqlite3").exists():
        return {"plans": []}
    with _ledger(project) as db:
        records = [_read_version(db, p, v) for p, v in db.execute("SELECT plan_id,current_version FROM plans ORDER BY rowid")]
        plans = []
        for record in records:
            history = [{"version": v, "plan_sha256": checksum, "status": json.loads(payload)["status"]}
                       for v, checksum, payload in db.execute(
                           "SELECT version,sha256,payload FROM versions WHERE plan_id=? ORDER BY version", (record["plan_id"],))]
            plans.append({k: record[k] for k in ("plan_id", "version", "status", "plan_sha256", "created_at")}
                         | {"name": record["spec"]["name"], "versions": history})
        return {"plans": plans}


def hypothesis_context(project):
    """Lightweight shared-state token; reads hashes/status, never count payloads."""
    path = project.root / "hypotheses.sqlite3"
    if not path.exists():
        return {"plan_count": 0, "latest": None, "state_sha256": _hash([])}
    if path.is_symlink() or getattr(path, "is_junction", lambda: False)():
        raise SpatialError("Hypothesis ledger must be a regular project file.")
    with sqlite3.connect(project.db_path.as_uri() + "?mode=ro", uri=True) as db:
        source_hash = db.execute("SELECT sha256 FROM source WHERE id=1").fetchone()[0]
        project_id = db.execute("SELECT value FROM state WHERE key='project_id'").fetchone()[0]
    with sqlite3.connect(path.as_uri() + "?mode=ro", uri=True) as db:
        identity = dict(db.execute("SELECT key,value FROM identity"))
        if identity != {"source_sha256": source_hash, "project_id": project_id}:
            raise SpatialError("Hypothesis ledger belongs to a different source or project.")
        heads = db.execute("""SELECT v.plan_id,v.version,v.sha256,json_extract(v.payload,'$.status')
            FROM plans p JOIN versions v ON p.plan_id=v.plan_id AND p.current_version=v.version
            ORDER BY v.rowid DESC""").fetchall()
    latest = dict(zip(("plan_id", "version", "plan_sha256", "status"), heads[0])) if heads else None
    return {"plan_count": len(heads), "latest": latest, "state_sha256": _hash(heads)}


def _selection(project, sid, revision_id, summary):
    selection = project.get_selection(sid)
    if (not selection or selection["revision_id"] != revision_id
            or selection.get("source_sha256") != summary["source_sha256"]
            or any(selection.get(k) != summary["metadata"].get(k) for k in ("slice_id", "coordinate_system", "units"))):
        raise SpatialError("Plan selection must match the pinned revision, source and coordinate frame.")
    return {k: deepcopy(selection[k]) for k in
            ("selection_id", "revision_id", "source_sha256", "cell_ids", "slice_id", "coordinate_system", "units")}


def _clause(cell, clause, panel):
    field, op, expected = clause["field"], clause["op"], clause["value"]
    if field.startswith("counts."):
        feature = field[7:]
        if feature not in panel:
            return None
        actual = cell["counts"].get(feature, 0)
    elif field.startswith("attributes."):
        key = field[11:]
        if key not in cell.get("attributes", {}):
            return None
        actual = cell["attributes"][key]
    else:
        if field not in cell:
            return None
        actual = cell[field]
    if actual is None:
        return None
    if op in {"gt", "ge", "lt", "le"}:
        if type(actual) not in (int, float):
            return None
        return {"gt": actual > expected, "ge": actual >= expected,
                "lt": actual < expected, "le": actual <= expected}[op]
    def equal(a, b):
        if type(a) in (int, float) and type(b) in (int, float):
            return a == b
        if type(a) is type(b) and type(a) in (bool, str):
            return a == b
        return None
    if op == "in":
        matches = [equal(actual, item) for item in expected]
        return True if True in matches else None if None in matches else False
    same = equal(actual, expected)
    return same if op == "eq" or same is None else not same


def _resolve(project, spec):
    summary = project.summary()
    if summary["cell_count"] > MAX_PROJECT_CELLS:
        raise SpatialError("Project exceeds the bounded hypothesis-plan observation limit.")
    revision = project.get_revision(spec["revision_id"])
    selected = _selection(project, spec["selection_id"], spec["revision_id"], summary)
    cells = project.cells(spec["revision_id"])
    universe = {c["cell_id"] for c in cells}
    if not selected["cell_ids"] or not set(selected["cell_ids"]) <= universe:
        raise SpatialError("Plan requires a nonempty selection of known observations.")
    background = spec["background"]
    background_selection = None
    if background["kind"] == "selection":
        background_selection = _selection(project, background["selection_id"], spec["revision_id"], summary)
        if not set(background_selection["cell_ids"]) <= universe:
            raise SpatialError("Background contains unknown observations.")
    panel = set(summary["metadata"]["panel_genes"])
    resolved_variants = []
    for variant in spec["variants"]:
        selector = variant["selector"]
        chosen, unknown = [], []
        if "cell_ids" in selector:
            if not set(selector["cell_ids"]) <= universe:
                raise SpatialError("Explicit hypothesis IDs include unknown observations.")
            chosen = selector["cell_ids"]
        else:
            for cell in cells:
                truth = [_clause(cell, clause, panel) for clause in selector["predicate"]["all"]]
                if False in truth:
                    continue
                (unknown if None in truth else chosen).append(cell["cell_id"])
        resolved_variants.append({"name": variant["name"], "cell_ids": sorted(chosen),
                                  "unknown_cell_ids": sorted(unknown), "changes": deepcopy(variant["changes"]),
                                  "resolution_status": "unknown" if unknown else "resolved",
                                  "matched_count": len(chosen), "unknown_count": len(unknown)})
    return {"source_sha256": summary["source_sha256"], "project_id": summary["project_id"],
            "revision_sha256": revision["revision_sha256"], "selection": selected,
            "background_selection": background_selection, "variants": resolved_variants,
            "predicate_semantics": "all clauses; false dominates missing; missing/null/type-incompatible values are unknown; only measured exact feature IDs use zero for absent counts"}


def freeze_hypothesis_plan(project, plan_id, expected_version):
    draft = get_hypothesis_plan(project, plan_id)
    if draft["version"] != expected_version:
        raise SpatialError("Stale hypothesis plan version; refresh before freezing.")
    if draft["status"] != "draft":
        raise SpatialError("This version is already frozen; create a revised draft before freezing again.")
    spec = _validate_spec(draft["spec"])
    resolved = _resolve(project, spec)
    with _ledger(project) as db:
        return _append_version(db, plan_id, expected_version, spec, "frozen", resolved)


def _unavailable(status, reason, *, graph=None, composition=None):
    return {"status": status, "reason": reason, "graph": graph,
            "composition": composition, "methods": {name: {
                "status": status, "reason": reason, "observed_fraction": None,
                "null_fraction": None, "excess_over_abundance": None, "classification": "unknown",
            } for name in METHODS}}


def _evaluate_v2(cells, selected, *, radius, scope, source, target, threshold, background, background_ids):
    active = [c for c in cells if c["included"]]
    nodes = active if scope == "whole_slice" else [c for c in active if c["cell_id"] in selected]
    sources = [i for i, c in enumerate(nodes) if c["cell_id"] in selected and c["label"] == source]
    bg = nodes if background == "neighbor_universe" else (
        active if background == "imported_universe" else [c for c in active if c["cell_id"] in background_ids])
    bg_ids = {c["cell_id"] for c in bg}
    bg_targets = sum(c["label"] == target for c in bg)
    graph = {"node_count": len(nodes), "source_count": len(sources),
             "target_count": sum(c["label"] == target for c in nodes), "background_count": len(bg),
             "background_target_count": bg_targets, "directed_source_edges": 0, "target_edges": 0,
             "isolated_sources": len(sources), "sources_with_neighbors": 0, "max_source_degree": 0,
             "edge_definition": "distinct non-self observations within Euclidean radius_um"}
    composition = _composition(cells, selected)
    if not sources:
        return _unavailable("unknown", "source_label_absent_from_included_selection", graph=graph, composition=composition)
    xy = np.array([[c["x"], c["y"]] for c in nodes], dtype=float)
    tree = cKDTree(xy)
    degrees = tree.query_ball_point(xy[sources], radius, return_length=True, workers=1) - 1
    edge_count, max_degree = int(np.sum(degrees)), int(np.max(degrees))
    if edge_count > MAX_DIRECTED_SOURCE_EDGES or max_degree > MAX_NEIGHBORS_PER_SOURCE:
        raise SpatialError("Graph resource budget exceeded; reduce the declared radius or source selection.")
    graph.update(directed_source_edges=edge_count, max_source_degree=max_degree,
                 isolated_sources=int(np.count_nonzero(degrees == 0)),
                 sources_with_neighbors=int(np.count_nonzero(degrees > 0)))
    observed_per_source, null_per_source, edge_weights = [], [], []
    target_flags = np.array([c["label"] == target for c in nodes], dtype=bool)
    undefined_null = 0
    for i, degree in zip(sources, degrees):
        if degree == 0:
            continue
        neighbors = tree.query_ball_point(xy[i], radius, workers=1, return_sorted=True)
        target_edges = sum(bool(target_flags[j]) for j in neighbors if j != i)
        graph["target_edges"] += target_edges
        in_background = nodes[i]["cell_id"] in bg_ids
        denominator = len(bg) - int(in_background)
        numerator = bg_targets - int(in_background and source == target)
        if denominator <= 0:
            undefined_null += 1
            continue
        observed_per_source.append(target_edges / int(degree))
        null_per_source.append(numerator / denominator)
        edge_weights.append(int(degree))
    graph["sources_with_unavailable_background"] = undefined_null
    if not graph["target_count"]:
        reason = "target_label_absent_from_neighbor_universe"
    elif not edge_count:
        reason = "no_eligible_source_neighbor_edges"
    elif source == target and graph["target_count"] <= 1:
        reason = "nonself_target_label_absent_from_neighbor_universe"
    elif undefined_null or not bg:
        reason = "nonself_target_abundance_background_unavailable"
    else:
        reason = None
    if reason:
        return _unavailable("unknown", reason, graph=graph, composition=composition)
    graph["nonisolated_source_background_fraction_range"] = [min(null_per_source), max(null_per_source)]
    graph["background_self_exclusion"] = "Remove the source only if it belongs to the background; remove its target membership only for same-label comparisons."
    methods = {}
    for name in METHODS:
        weights = edge_weights if name == "edge_weighted" else None
        observed = float(np.average(observed_per_source, weights=weights))
        null = float(np.average(null_per_source, weights=weights))
        excess = observed - null
        methods[name] = {"status": "computed", "reason": None, "observed_fraction": observed,
                         "null_fraction": null, "excess_over_abundance": excess,
                         "classification": "descriptive_supported" if excess >= threshold else "descriptive_not_supported",
                         "denominator": edge_count if weights is not None else graph["sources_with_neighbors"],
                         "weighting": "one vote per directed edge" if weights is not None else "one vote per source with at least one neighbor",
                         "isolated_source_policy": "excluded from the estimand; count retained explicitly"}
    return {"status": "computed", "reason": None, "graph": graph, "composition": composition, "methods": methods}


def _calculate(cells, selected, **kwargs):
    try:
        return _evaluate_v2(cells, selected, **kwargs)
    except Exception as exc:
        # One failed grid cell must not erase all other hypotheses. Cancellation
        # and process exits are BaseException and are intentionally not caught.
        result = _unavailable("failed", "resource_limit" if "resource budget" in str(exc) else "calculation_error")
        result["error"] = {"type": type(exc).__name__, "message": str(exc)[:1000]}
        return result


def run_hypothesis_plan(project, plan_id=None, version=None, *, frozen_plan=None):
    """Run every frozen cell; reference methods are alternate estimands, not truth."""
    if frozen_plan is None:
        if plan_id is None or version is None:
            raise SpatialError("Execution requires an explicit frozen plan_id and version.")
        frozen_plan = get_hypothesis_plan(project, plan_id, version)
    elif plan_id is not None or version is not None:
        raise SpatialError("Supply a ledger version or embedded frozen_plan, not both.")
    record = deepcopy(frozen_plan)
    if not isinstance(record, dict) or record.get("schema") != PLAN_SCHEMA or record.get("status") != "frozen":
        raise SpatialError("Only an explicitly frozen plan can execute.")
    checksum = record.pop("plan_sha256", None)
    if _hash(record) != checksum or record.get("scientific_authorization") != "NOT_ESTABLISHED":
        raise SpatialError("Frozen plan integrity or scientific authorization is invalid.")
    record["plan_sha256"] = checksum
    spec = _validate_spec(record["spec"])
    resolved = _resolve(project, spec)
    if resolved != record["resolved"]:
        raise SpatialError("Frozen plan IDs, selectors or source anchors do not match the pinned inputs.")
    summary = project.summary()
    cells = project.cells(spec["revision_id"])
    selected = set(resolved["selection"]["cell_ids"])
    background_ids = set((resolved["background_selection"] or {}).get("cell_ids", []))
    calibrated = summary["metadata"].get("units") == "micrometer"
    baseline = {}
    for scope in spec["graph_scopes"]:
        for radius in spec["radii_um"]:
            baseline[scope, radius] = _calculate(
                cells, selected, radius=radius, scope=scope, source=spec["source_label"], target=spec["target_label"],
                threshold=spec["min_effect"], background=spec["background"]["kind"], background_ids=background_ids,
            ) if calibrated else _unavailable("unknown", "physical_coordinate_calibration_unavailable")
    results = []
    for variant in resolved["variants"]:
        chosen = set(variant["cell_ids"])
        altered = [dict(c, **variant["changes"]) if c["cell_id"] in chosen else c for c in cells]
        changed = [c["cell_id"] for c in cells if c["cell_id"] in chosen
                   and any(c[k] != v for k, v in variant["changes"].items())]
        for scope in spec["graph_scopes"]:
            for radius in spec["radii_um"]:
                before = baseline[scope, radius]
                if not calibrated:
                    after = _unavailable("unknown", "physical_coordinate_calibration_unavailable")
                elif variant["unknown_count"]:
                    after = _unavailable("unknown", "hypothesis_selector_has_unknown_observations")
                else:
                    after = _calculate(altered, selected, radius=radius, scope=scope, source=spec["source_label"],
                                       target=spec["target_label"], threshold=spec["min_effect"],
                                       background=spec["background"]["kind"], background_ids=background_ids)
                comparisons = {}
                for method in METHODS:
                    left, right = before["methods"][method], after["methods"][method]
                    available = left["status"] == right["status"] == "computed"
                    comparisons[method] = {"status": ("stable" if left["classification"] == right["classification"] else "changed") if available else "indeterminate",
                                           "effect_delta": right["excess_over_abundance"] - left["excess_over_abundance"] if available else None,
                                           "threshold_crossed": left["classification"] != right["classification"] if available else None}
                state = "failed" if "failed" in (before["status"], after["status"]) else (
                    "unknown" if "unknown" in (before["status"], after["status"]) else "computed")
                results.append({"row_id": f"row_{len(results) + 1:04d}", "variant": variant["name"],
                                "radius_um": radius, "graph_scope": scope, "background": deepcopy(spec["background"]),
                                "status": state, "matched_count": len(chosen), "selector_unknown_count": variant["unknown_count"],
                                "actual_changed_count": len(changed), "changed_in_analysis_roi": len(selected.intersection(changed)),
                                "before": deepcopy(before), "after": after, "comparison": comparisons,
                                "reference_comparison": {phase: {
                                    "source_equal_minus_edge_weighted": value["methods"][METHODS[1]]["excess_over_abundance"] - value["methods"][METHODS[0]]["excess_over_abundance"]
                                    if value["status"] == "computed" else None,
                                    "interpretation": "Different weighting estimands; disagreement is not an accuracy ranking."}
                                    for phase, value in (("before", before), ("after", after))}})
    def difference(row):
        values = [v["effect_delta"] for v in row["comparison"].values() if v["effect_delta"] is not None]
        return max(map(abs, values)) if values else None
    ordered = sorted(results, key=lambda row: (difference(row) is None, -(difference(row) or 0), row["row_id"]))
    result = {"analysis_schema": RUN_SCHEMA, "base_revision": spec["revision_id"], "target_revision": spec["revision_id"],
              "selection_id": spec["selection_id"], "selection": resolved["selection"],
              "parameters": {"frozen_plan": record}, "plan_sha256": checksum, "results": results,
              "difference_order": [r["row_id"] for r in ordered],
              "difference_order_rule": "All rows, descending max absolute effect delta across both methods; unknown/failed remain at end; stable IDs break ties. Not a winner selection.",
              "status_counts": dict(Counter(r["status"] for r in results)),
              "semantics": semantics(summary["metadata"]), "import_scope": summary["metadata"].get("import_scope", {}),
              "methods": {"primary": METHODS[0], "reference": METHODS[1], "reference_role": "Alternative descriptive source weighting; not independent validation or a significance test."},
              "revision_created": False, "researcher_approval": "NOT_REQUESTED_HYPOTHESIS_ONLY",
              "scientific_authorization": "NOT_ESTABLISHED", "evidence_ceiling": "single_slice_frozen_hypothesis_descriptive_sensitivity",
              "data_hashes": {"source_sha256": resolved["source_sha256"], "base_revision_sha256": resolved["revision_sha256"],
                              "target_revision_sha256": resolved["revision_sha256"]},
              "provenance": {**_runtime_provenance(), "algorithm_version": "frozen-plan-dual-weighting.v2"},
              "limitations": ["Freezing records computational choices, not a named person's scientific approval or independent preregistration.",
                              "Every variant starts from the same pinned baseline. Every declared grid cell is retained, including unknowns and failures.",
                              "Equal-source weighting excludes isolated sources and reports their count; it is a different estimand, not a more accurate method.",
                              "Background abundance is recomputed after each hypothetical inclusion/label change in the explicitly declared universe; self is excluded per source.",
                              "Missing predicate data and unmeasured exact feature IDs remain unknown, never silently false or zero.",
                              "Physical neighborhoods require calibrated coordinates; absent calibration yields unknown cells throughout the grid.",
                              "Threshold and method differences do not establish cell identities, biological interaction, causal effects or population inference."]}
    return project.save_run(result)
