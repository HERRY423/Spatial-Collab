"""Content-addressed extension records in the existing project database."""
from __future__ import annotations

import json

from .store import SpatialError, _hash, _json, _now


def _exists(db):
    return db.execute("SELECT 1 FROM sqlite_master WHERE name='research_objects'").fetchone()


def put(project, kind, payload):
    body = {**payload, "object_kind": kind}
    digest = _hash(body)
    identifier = kind + "_" + digest[:32]
    if hasattr(project, "research_objects"):
        return {**body, "object_id": identifier, "object_sha256": digest}
    with project._db(write=True) as db:
        db.execute("CREATE TABLE IF NOT EXISTS research_objects "
                   "(id TEXT PRIMARY KEY, kind TEXT NOT NULL, payload TEXT NOT NULL, sha256 TEXT NOT NULL, created_at TEXT NOT NULL)")
        for action in ("UPDATE", "DELETE"):
            db.execute(f"CREATE TRIGGER IF NOT EXISTS research_objects_no_{action.lower()} BEFORE {action} "
                       "ON research_objects BEGIN SELECT RAISE(ABORT,'Immutable research object'); END")
        db.execute("INSERT OR IGNORE INTO research_objects VALUES(?,?,?,?,?)",
                   (identifier, kind, _json(body), digest, _now()))
    return {**body, "object_id": identifier, "object_sha256": digest}


def get(project, identifier, kind=None):
    if hasattr(project, "research_objects"):
        result = project.research_objects.get(identifier)
        if result is None or kind and result["object_kind"] != kind:
            raise SpatialError("Unknown exported research object.")
        body = {k: v for k, v in result.items() if k not in {"object_id", "object_sha256"}}
        if _hash(body) != result["object_sha256"]:
            raise SpatialError("Exported research object integrity failed.")
        return json.loads(_json(result))
    with project._db() as db:
        row = db.execute("SELECT * FROM research_objects WHERE id=?", (identifier,)).fetchone() if _exists(db) else None
    if row is None or (kind and row["kind"] != kind):
        raise SpatialError("Unknown research object or wrong object type.")
    body = json.loads(row["payload"])
    if _hash(body) != row["sha256"] or identifier != row["kind"] + "_" + row["sha256"][:32]:
        raise SpatialError("Research object integrity check failed.")
    return {**body, "object_id": identifier, "object_sha256": row["sha256"]}


def catalog(project, kind):
    with project._db() as db:
        return [r[0] for r in db.execute("SELECT id FROM research_objects WHERE kind=? ORDER BY rowid", (kind,))] if _exists(db) else []
