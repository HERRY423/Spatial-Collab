"""Transactional research workspace with immutable imports and revision overlays.

Reviewer names are user attribution, not authentication or scientific approval.
Only a revision's label, inclusion and region overlays can change source cells.
"""
from __future__ import annotations

from collections import Counter
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sqlite3
from typing import Any
from uuid import uuid4


class SpatialError(ValueError):
    """A readable domain/validation error safe to show in the workbench."""


def _json(value: Any) -> str:
    try:
        return json.dumps(value, ensure_ascii=False, sort_keys=True,
                          separators=(",", ":"), allow_nan=False)
    except (ValueError, TypeError) as exc:
        raise SpatialError("Data must contain finite JSON-compatible values.") from exc


def _hash(value: Any) -> str:
    return hashlib.sha256(_json(value).encode("utf-8")).hexdigest()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _id(prefix: str) -> str:
    return prefix + "_" + uuid4().hex


def _text(value: Any, field: str, *, empty: bool = False, limit: int = 4000) -> str:
    if not isinstance(value, str) or len(value) > limit or (not empty and not value.strip()):
        raise SpatialError(f"{field} must be {'a' if empty else 'a nonempty'} string of at most {limit} characters.")
    return value


def _number(value: Any, field: str, *, nonnegative: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SpatialError(f"{field} must be a finite number.")
    try:
        number = float(value)
    except (OverflowError, ValueError):
        raise SpatialError(f"{field} must be a finite number.") from None
    if not math.isfinite(number) or (nonnegative and number < 0):
        raise SpatialError(f"{field} must be finite" + (" and nonnegative." if nonnegative else "."))
    return number


def _attributes(value: Any) -> dict:
    """Preserve bounded source-reported scalars without inventing QC meaning."""
    if not isinstance(value, dict) or len(value) > 32:
        raise SpatialError("attributes must be an object containing at most 32 scalar fields.")
    result = {}
    for key, item in value.items():
        _text(key, "attribute name", limit=128)
        if item is None or type(item) is bool:
            result[key] = item
        elif type(item) is str:
            result[key] = _text(item, f"attribute {key}", empty=True, limit=1000)
        elif type(item) in (int, float):
            _number(item, f"attribute {key}")
            result[key] = item
        else:
            raise SpatialError(f"attribute {key} must be a JSON scalar, not a nested object or list.")
    return result


def _validate_import(cells: Any, metadata: Any) -> dict:
    if not isinstance(metadata, dict):
        raise SpatialError("metadata must be an object.")
    meta = json.loads(_json(metadata))
    for key in ("name", "slice_id", "coordinate_system", "source_kind"):
        _text(meta.get(key), key)
    from .coordinates import UNITS
    if not isinstance(meta.get("units"), str) or meta["units"] not in UNITS:
        raise SpatialError("Declare units as micrometer, pixel, array_index or unknown; never infer physical scale.")
    # Optional for legacy snapshots: defaults belong in presentation, not in the
    # immutable payload, whose canonical hash must remain replayable.
    if "observation_unit" in meta:
        if meta["observation_unit"] not in {"cell", "spot", "bin"}:
            raise SpatialError("observation_unit must be cell, spot or bin.")
        semantics = {"cell": "cell_type", "spot": "spot_annotation", "bin": "bin_annotation"}
        if meta.get("label_semantics") != semantics[meta["observation_unit"]]:
            raise SpatialError("label_semantics must agree with observation_unit.")
    elif "label_semantics" in meta:
        raise SpatialError("label_semantics requires observation_unit.")
    if "platform" in meta:
        _text(meta["platform"], "platform")
    panel = meta.get("panel_genes")
    if not isinstance(panel, list) or any(not isinstance(g, str) or not g.strip() for g in panel):
        raise SpatialError("panel_genes must be a list of nonempty gene names.")
    if len(set(panel)) != len(panel):
        raise SpatialError("panel_genes must not contain duplicates.")
    meta["panel_genes"] = sorted(panel)
    if "features" in meta:
        from .identity import validate_features
        validate_features(meta)
    if meta.get("identity_scope") == "sample_qualified":
        _text(meta.get("sample_id"), "metadata.sample_id", limit=500)
    replicates = meta.get("biological_replicates", 0)
    if type(replicates) is not int or replicates not in (0, 1):
        raise SpatialError("For this single-slice workspace biological_replicates must be 0 (unknown) or 1; cells are not biological replicates.")
    meta["biological_replicates"] = replicates
    limitations = meta.get("limitations", [])
    if not isinstance(limitations, list) or any(not isinstance(x, str) for x in limitations):
        raise SpatialError("limitations must be a list of strings.")
    meta["limitations"] = limitations
    source_files = meta.get("source_files", [])
    if not isinstance(source_files, list):
        raise SpatialError("source_files must be a list.")
    for item in source_files:
        if not isinstance(item, dict):
            raise SpatialError("Each source file must record a path and SHA256.")
        _text(item.get("path"), "source_files.path")
        digest = item.get("sha256")
        if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdefABCDEF" for c in digest):
            raise SpatialError("Each source file must record a valid SHA256.")
    meta["source_files"] = source_files
    empty_workspace = meta.get("source_kind") == "atlas_workspace"
    if empty_workspace and (cells != [] or panel or replicates or meta["units"] != "unknown"):
        raise SpatialError("An atlas_workspace is an empty container, with no primary observations, panel or physical coordinate claim.")
    if not isinstance(cells, list) or (not cells and not empty_workspace):
        raise SpatialError("cells must be a nonempty list.")
    identifiers, clean = set(), []
    genes = set(panel)
    for cell in cells:
        if not isinstance(cell, dict):
            raise SpatialError("Each cell must be an object.")
        cid = _text(cell.get("cell_id"), "cell_id", limit=512)
        if cid in identifiers:
            raise SpatialError(f"Duplicate cell_id: {cid}")
        identifiers.add(cid)
        included = cell.get("included", True)
        if type(included) is not bool:
            raise SpatialError("included must be true or false.")
        counts = cell.get("counts", {})
        if not isinstance(counts, dict):
            raise SpatialError("counts must be a gene-to-count object.")
        if any(g not in genes for g in counts):
            raise SpatialError("Every counts gene must be declared in panel_genes.")
        clean_counts = {g: _number(v, f"count for {g}", nonnegative=True) for g, v in counts.items()}
        clean.append({"cell_id": cid, "x": _number(cell.get("x"), "x"),
                      "y": _number(cell.get("y"), "y"),
                      "label": _text(cell.get("label"), "label", limit=512),
                      "included": included, "region": _text(cell.get("region", ""), "region", empty=True, limit=512),
                      "counts": clean_counts})
        # Absence stays absent: adding empty defaults would change canonical
        # snapshots and source hashes of all pre-attributes projects.
        if "attributes" in cell:
            clean[-1]["attributes"] = _attributes(cell["attributes"])
        if meta.get("identity_scope") == "sample_qualified" or "sample_id" in cell or "source_cell_id" in cell:
            from .identity import qualify_observation_id
            sample = _text(cell.get("sample_id"), "sample_id", limit=500)
            original_id = _text(cell.get("source_cell_id"), "source_cell_id", limit=500)
            if "sample_id" in meta and sample != meta["sample_id"]:
                raise SpatialError("Observation sample_id must match metadata.sample_id in this single-slice project.")
            if qualify_observation_id(sample, original_id) != cid:
                raise SpatialError("Sample-qualified identity does not match sample_id/source_cell_id.")
            clean[-1].update(sample_id=sample, source_cell_id=original_id)
    return {"metadata": meta, "cells": sorted(clean, key=lambda c: c["cell_id"])}


def _on_segment(p: list[float], a: list[float], b: list[float]) -> bool:
    cross = (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
    tolerance = 1e-10 * max(1.0, abs(b[0] - a[0]), abs(b[1] - a[1]))
    return (abs(cross) <= tolerance and min(a[0], b[0]) - tolerance <= p[0] <= max(a[0], b[0]) + tolerance
            and min(a[1], b[1]) - tolerance <= p[1] <= max(a[1], b[1]) + tolerance)


def _polygon(value: Any) -> list[list[float]]:
    if not isinstance(value, list) or not 3 <= len(value) <= 1001:
        raise SpatialError("polygon needs 3 to 1000 vertices (plus optional closing vertex).")
    points = []
    for point in value:
        if not isinstance(point, (list, tuple)) or len(point) != 2:
            raise SpatialError("Each polygon vertex must contain x and y.")
        points.append([_number(point[0], "polygon x"), _number(point[1], "polygon y")])
    if points[0] == points[-1]:
        points.pop()
    if not 3 <= len(points) <= 1000 or len(set(map(tuple, points))) != len(points):
        raise SpatialError("polygon needs distinct vertices.")
    # Translate before multiplying: absolute-coordinate shoelace loses small
    # ROI areas when the coordinate frame has a large origin offset.
    origin = points[0]
    shifted = [[p[0] - origin[0], p[1] - origin[1]] for p in points]
    try:
        area = math.fsum(p[0] * q[1] - q[0] * p[1] for p, q in zip(shifted, shifted[1:] + shifted[:1]))
    except (OverflowError, ValueError):
        raise SpatialError("polygon coordinate range is too large.") from None
    if not math.isfinite(area) or abs(area) <= 1e-12:
        raise SpatialError("polygon must enclose a nonzero finite area.")
    # Self-intersecting polygons have ambiguous biological regions: reject them.
    def turn(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    n = len(points)
    for i in range(n):
        a, b = points[i], points[(i + 1) % n]
        for j in range(i + 1, n):
            if j == i + 1 or (i == 0 and j == n - 1):
                continue
            c, d = points[j], points[(j + 1) % n]
            cross = turn(a, b, c) * turn(a, b, d) < 0 and turn(c, d, a) * turn(c, d, b) < 0
            if cross or any((_on_segment(c, a, b), _on_segment(d, a, b), _on_segment(a, c, d), _on_segment(b, c, d))):
                raise SpatialError("polygon must be simple, without self-intersections.")
    return points


def _contains(x: float, y: float, polygon: list[list[float]]) -> bool:
    inside = False
    for a, b in zip(polygon, polygon[1:] + polygon[:1]):
        if _on_segment([x, y], a, b):
            return True
        if (a[1] > y) != (b[1] > y) and x < (b[0] - a[0]) * (y - a[1]) / (b[1] - a[1]) + a[0]:
            inside = not inside
    return inside


class Project:
    """A single-slice local workspace. Connections are never shared by threads."""

    def __init__(self, root: str | Path):
        self.root = Path(root).resolve()
        self.db_path = self.root / "project.sqlite3"
        if not self.db_path.is_file():
            raise SpatialError(f"No Spatial Collab project at {self.root}")
        try:
            with self._db() as db:
                source = db.execute("SELECT payload,sha256 FROM source WHERE id=1").fetchone()
                if not source or _hash(json.loads(source["payload"])) != source["sha256"]:
                    raise SpatialError("Immutable input snapshot failed its integrity check.")
                version = db.execute("SELECT value FROM state WHERE key='schema_version'").fetchone()
                if not version or version[0] != "1":
                    raise SpatialError("Unsupported project schema version.")
        except sqlite3.Error as exc:
            raise SpatialError("Invalid or unreadable Spatial Collab project.") from exc

    @contextmanager
    def _db(self, *, write: bool = False):
        try:
            db = sqlite3.connect(self.db_path.as_uri() + "?mode=rw", uri=True, timeout=30, isolation_level=None)
        except sqlite3.Error as exc:
            raise SpatialError("Project database is missing or unavailable.") from exc
        db.row_factory = sqlite3.Row
        db.execute("PRAGMA foreign_keys=ON")
        db.execute("PRAGMA busy_timeout=30000")
        try:
            db.execute("BEGIN IMMEDIATE" if write else "BEGIN")
            yield db
            db.commit()
        except sqlite3.Error as exc:
            db.rollback()
            raise SpatialError("Project database is unavailable or inconsistent; retry after other writes complete.") from exc
        except Exception:
            db.rollback()
            raise
        finally:
            db.close()

    @classmethod
    def create_workspace(cls, root: str | Path, name: str = "Spatial atlas workspace") -> "Project":
        """Start directly with disk-backed datasets, without invented primary cells."""
        return cls.create(root, [], {"name": name, "slice_id": "no_primary_slice",
                         "source_kind": "atlas_workspace", "coordinate_system": "no_primary_frame",
                         "units": "unknown", "panel_genes": [], "biological_replicates": 0,
                         "limitations": ["Empty container: select an explicitly identified atlas for observations, coordinates and ROI analyses."]})

    @classmethod
    def create(cls, root: str | Path, cells: list[dict], metadata: dict) -> "Project":
        snapshot = _validate_import(cells, metadata)
        path = Path(root).resolve()
        path.mkdir(parents=True, exist_ok=True)
        db_path = path / "project.sqlite3"
        try:
            with db_path.open("xb"):
                pass
        except FileExistsError as exc:
            raise SpatialError("Project already exists; imports never overwrite a project.") from exc
        db = sqlite3.connect(db_path, isolation_level=None)
        try:
            db.execute("PRAGMA journal_mode=WAL")
            db.execute("PRAGMA synchronous=FULL")
            db.executescript("""
                CREATE TABLE state(key TEXT PRIMARY KEY,value TEXT NOT NULL);
                CREATE TABLE source(id INTEGER PRIMARY KEY CHECK(id=1),payload TEXT NOT NULL,sha256 TEXT NOT NULL);
                CREATE TABLE revisions(id TEXT PRIMARY KEY,payload TEXT NOT NULL,sha256 TEXT NOT NULL,proposal_id TEXT UNIQUE);
                CREATE TABLE selections(id TEXT PRIMARY KEY,payload TEXT NOT NULL);
                CREATE TABLE proposals(id TEXT PRIMARY KEY,payload TEXT NOT NULL);
                CREATE TABLE runs(id TEXT PRIMARY KEY,payload TEXT NOT NULL,sha256 TEXT NOT NULL);
            """)
            for table in ("source", "revisions", "selections", "proposals", "runs"):
                for action in ("UPDATE", "DELETE"):
                    db.execute(f"CREATE TRIGGER {table}_no_{action.lower()} BEFORE {action} ON {table} BEGIN SELECT RAISE(ABORT,'Immutable record'); END")
            project_id, revision_id = _id("project"), _id("rev")
            source_hash = _hash(snapshot)
            revision = {"revision_id": revision_id, "parent_revision": None, "created_at": _now(),
                        "kind": "import", "overlays": {}, "rationale": "Immutable source import",
                        "reviewer": None, "source_sha256": source_hash,
                        "scientific_authorization": "NOT_ESTABLISHED"}
            db.execute("BEGIN IMMEDIATE")
            db.execute("INSERT INTO source VALUES(1,?,?)", (_json(snapshot), source_hash))
            db.execute("INSERT INTO revisions VALUES(?,?,?,NULL)", (revision_id, _json(revision), _hash(revision)))
            db.executemany("INSERT INTO state VALUES(?,?)", [("schema_version", "1"), ("project_id", project_id),
                           ("head_revision", revision_id), ("active_selection", "")])
            db.commit()
        except Exception:
            db.rollback()
            db.close()
            # Only this exclusively created file belongs to this failed import.
            db_path.unlink(missing_ok=True)
            raise
        finally:
            db.close()
        return cls(path)

    @staticmethod
    def _state(db, key: str) -> str:
        row = db.execute("SELECT value FROM state WHERE key=?", (key,)).fetchone()
        if row is None:
            raise SpatialError("Project state is incomplete.")
        return row[0]

    @staticmethod
    def _source(db) -> tuple[dict, str]:
        row = db.execute("SELECT payload,sha256 FROM source WHERE id=1").fetchone()
        source = json.loads(row["payload"])
        if _hash(source) != row["sha256"]:
            raise SpatialError("Immutable input snapshot failed its integrity check.")
        return source, row["sha256"]

    @staticmethod
    def _record(db, table: str, identifier: str) -> dict:
        _text(identifier, f"{table} ID", limit=512)
        row = db.execute(f"SELECT * FROM {table} WHERE id=?", (identifier,)).fetchone()
        if not row:
            raise SpatialError(f"Unknown {table} ID: {identifier}")
        payload = json.loads(row["payload"])
        if "sha256" in row.keys() and _hash(payload) != row["sha256"]:
            raise SpatialError(f"{table} record failed its integrity check.")
        return payload

    @classmethod
    def _head(cls, db, expected: str) -> str:
        _text(expected, "expected_revision", limit=512)
        head = cls._state(db, "head_revision")
        if expected != head:
            raise SpatialError(f"Stale revision: expected {expected}, current head is {head}. Refresh and review again.")
        return head

    @classmethod
    def _materialize(cls, db, revision_id: str) -> list[dict]:
        source, _ = cls._source(db)
        revision = cls._record(db, "revisions", revision_id)
        for cell in source["cells"]:
            cell.update(revision["overlays"].get(cell["cell_id"], {}))
        return source["cells"]

    def cells(self, revision_id: str | None = None) -> list[dict]:
        with self._db() as db:
            return self._materialize(db, revision_id or self._state(db, "head_revision"))

    def get_revision(self, revision_id: str) -> dict:
        with self._db() as db:
            result = self._record(db, "revisions", revision_id)
            result["revision_sha256"] = _hash(result)
            return result

    @classmethod
    def _selection(cls, db, selection_id: str | None = None) -> dict | None:
        if selection_id is not None:
            _text(selection_id, "selection_id", limit=512)
        sid = selection_id if selection_id is not None else cls._state(db, "active_selection")
        if not sid:
            return None
        result = cls._record(db, "selections", sid)
        result["stale"] = result["revision_id"] != cls._state(db, "head_revision")
        return result

    def get_selection(self, selection_id: str | None = None) -> dict | None:
        with self._db() as db:
            return self._selection(db, selection_id)

    def get_proposal(self, proposal_id: str) -> dict:
        """Read a stored proposal so another interface can review the same delta."""
        with self._db() as db:
            proposal = self._record(db, "proposals", proposal_id)
            applied = db.execute("SELECT id FROM revisions WHERE proposal_id=?", (proposal_id,)).fetchone()
            proposal["applied_revision"] = applied[0] if applied else None
            proposal["stale"] = proposal["base_revision"] != self._state(db, "head_revision")
            return proposal

    @classmethod
    def _run_status(cls, db, result: dict) -> dict:
        result["stale"] = result["target_revision"] != cls._state(db, "head_revision")
        result["applicability"] = "historical_revision" if result["stale"] else "current_revision"
        return result

    def summary(self) -> dict:
        with self._db() as db:
            source, digest = self._source(db)
            revisions = []
            for row in db.execute("SELECT payload,sha256 FROM revisions ORDER BY rowid"):
                rev = json.loads(row["payload"])
                revisions.append({k: v for k, v in rev.items() if k != "overlays"})
                revisions[-1]["revision_sha256"] = row["sha256"]
            runs = []
            for row in db.execute("SELECT payload,sha256 FROM runs ORDER BY rowid"):
                run = self._run_status(db, json.loads(row["payload"]))
                runs.append({k: run[k] for k in ("run_id", "base_revision", "target_revision", "created_at", "stale", "applicability")})
            return {"metadata": source["metadata"], "project_id": self._state(db, "project_id"),
                    "head_revision": self._state(db, "head_revision"), "source_sha256": digest,
                    "cell_count": len(source["cells"]), "revision_count": len(revisions),
                    "selection": self._selection(db), "revisions": revisions, "runs": runs,
                    "scientific_authorization": "NOT_ESTABLISHED"}

    def context(self, after_cursor: str | None = None) -> dict:
        """A cheap, coherent state token; no source matrix or full history read.

        This is a polling cursor, not an append-only event stream. A client must
        refresh the actual objects after a change, and still use CAS for edits.
        """
        if after_cursor is not None:
            _text(after_cursor, "after_cursor", limit=128)
        with self._db() as db:
            state = {"project_id": self._state(db, "project_id"),
                     "head_revision": self._state(db, "head_revision"),
                     "selection_id": self._state(db, "active_selection") or None,
                     "source_sha256": db.execute("SELECT sha256 FROM source WHERE id=1").fetchone()[0]}
            for table in ("revisions", "proposals", "runs"):
                row = db.execute(f"SELECT id FROM {table} ORDER BY rowid DESC LIMIT 1").fetchone()
                state["latest_" + table.rstrip("s") + "_id"] = row[0] if row else None
            if db.execute("SELECT 1 FROM sqlite_master WHERE name='research_objects'").fetchone():
                state["research_objects"] = [r[0] for r in db.execute("SELECT id FROM research_objects ORDER BY rowid")]
            if db.execute("SELECT 1 FROM sqlite_master WHERE name='integration_jobs'").fetchone():
                state["integration_jobs"] = [list(r) for r in db.execute("SELECT id,status,result_id FROM integration_jobs ORDER BY rowid")]
            atlas_view = db.execute("SELECT value FROM state WHERE key='atlas_view'").fetchone()
            if atlas_view:
                state["atlas_view"] = json.loads(atlas_view[0])
        from .hypotheses import hypothesis_context
        state["hypotheses"] = hypothesis_context(self)
        from .proteomics import assay_context
        state["assays"] = assay_context(self)
        cursor = _hash(state)
        return {**state, "cursor": cursor, "changed": after_cursor != cursor,
                "cursor_kind": "current_state_token", "scientific_authorization": "NOT_ESTABLISHED"}

    def set_selection(self, expected_revision: str, cell_ids: list[str] | None = None,
                      polygon: list[list[float]] | None = None, name: str = "Selection") -> dict:
        if (cell_ids is None) == (polygon is None):
            raise SpatialError("Provide exactly one of cell_ids or polygon.")
        name = _text(name, "selection name", limit=512)
        geometry = _polygon(polygon) if polygon is not None else None
        with self._db(write=True) as db:
            head = self._head(db, expected_revision)
            source, source_hash = self._source(db)
            available = {c["cell_id"] for c in source["cells"]}
            if cell_ids is not None:
                if not isinstance(cell_ids, list) or any(not isinstance(x, str) for x in cell_ids):
                    raise SpatialError("cell_ids must be a list of cell ID strings.")
                if len(set(cell_ids)) != len(cell_ids):
                    raise SpatialError("cell_ids must not contain duplicates.")
                unknown = set(cell_ids) - available
                if unknown:
                    raise SpatialError(f"Unknown cell IDs: {', '.join(sorted(unknown)[:10])}")
                chosen = sorted(cell_ids)
            else:
                chosen = [c["cell_id"] for c in source["cells"] if _contains(c["x"], c["y"], geometry)]
            meta = source["metadata"]
            result = {"selection_id": _id("sel"), "revision_id": head, "name": name,
                      "slice_id": meta["slice_id"], "coordinate_system": meta["coordinate_system"],
                      "units": meta["units"], "source_sha256": source_hash, "cell_ids": chosen, "cell_count": len(chosen),
                      "polygon": geometry, "created_at": _now(),
                      "method": "centroid_containment_boundary_included" if geometry else "exact_cell_ids",
                      "interpretation": "Selection identifies cell centroids; it does not edit segmentation."}
            if meta.get("observation_unit", "cell") != "cell":
                result["interpretation"] = "Selection identifies spot/bin centroids; it does not imply pure cell identity or edit segmentation."
            db.execute("INSERT INTO selections VALUES(?,?)", (result["selection_id"], _json(result)))
            db.execute("UPDATE state SET value=? WHERE key='active_selection'", (result["selection_id"],))
            return {**result, "stale": False}

    def inspect_selection(self, genes: list[str] | None = None) -> dict:
        with self._db() as db:
            selection = self._selection(db)
            if selection is None:
                raise SpatialError("Select cells first.")
            if selection["stale"]:
                raise SpatialError("Stale selection: select cells again at the current revision.")
            source, _ = self._source(db)
            panel = set(source["metadata"]["panel_genes"])
            if genes is None:
                genes = sorted(panel)[:20]
            if not isinstance(genes, list) or len(genes) > 100 or any(not isinstance(g, str) or not g.strip() for g in genes):
                raise SpatialError("genes must be a list of at most 100 nonempty gene names.")
            if len(set(genes)) != len(genes):
                raise SpatialError("genes must not contain duplicates.")
            selected = set(selection["cell_ids"])
            cells = [c for c in self._materialize(db, selection["revision_id"]) if c["cell_id"] in selected]
            expression = {}
            for gene in genes:
                resolution = None
                if "features" in source["metadata"] or gene.startswith(("feature_id:", "symbol:")):
                    from .identity import describe_feature_query
                    resolution = describe_feature_query(source["metadata"], gene)
                feature_id = resolution["feature_id"] if resolution else gene
                if feature_id not in panel:
                    expression[gene] = {"status": "unmeasured", "measured": False, "mean": None, "sum": None,
                                        "nonzero_cells": None, "zero_cells": None, "cell_count": len(cells)}
                    if resolution:
                        expression[gene].update(status=resolution["status"], feature_resolution=resolution)
                    continue
                # Panel membership is the import contract for sparse measured zeros.
                values = [c["counts"].get(feature_id, 0.0) for c in cells]
                try:
                    total = math.fsum(values)
                except OverflowError:
                    total = None
                expression[gene] = {"status": "measured" if cells else "no_selected_cells", "measured": True,
                                    "mean": math.fsum(v / len(values) for v in values) if values else None,
                                    "sum": total, "sum_status": "overflow" if total is None else "available", "nonzero_cells": sum(v > 0 for v in values),
                                    "zero_cells": sum(v == 0 for v in values), "cell_count": len(cells)}
                if resolution:
                    expression[gene]["feature_resolution"] = resolution
            from .exploration import semantics
            return {"selection": selection, "semantics": semantics(source["metadata"]),
                    "cell_count": len(cells), "label_counts": dict(sorted(Counter(c["label"] for c in cells).items())),
                    "included_count": sum(c["included"] for c in cells), "expression": expression,
                    "expression_summaries": expression, "cells": cells[:100], "returned_cell_count": min(100, len(cells)),
                    "truncated": len(cells) > 100, "expression_scope": "all selected observations, including excluded observations",
                    "scientific_authorization": "NOT_ESTABLISHED"}

    def propose_revision(self, expected_revision: str, selection_id: str, changes: dict,
                         rationale: str, actor: str = "agent") -> dict:
        if not isinstance(changes, dict) or not changes or set(changes) - {"label", "included", "region"}:
            raise SpatialError("changes may contain only label, included and region; provide at least one change.")
        if "included" in changes and type(changes["included"]) is not bool:
            raise SpatialError("included must be true or false.")
        for key in ("label", "region"):
            if key in changes:
                _text(changes[key], key, empty=key == "region", limit=512)
        _text(rationale, "rationale")
        _text(actor, "actor", limit=512)
        with self._db(write=True) as db:
            head = self._head(db, expected_revision)
            selection = self._selection(db, selection_id)
            if selection is None or selection["stale"]:
                raise SpatialError("Stale or absent selection: select cells at the current revision.")
            chosen = set(selection["cell_ids"])
            deltas = []
            for cell in self._materialize(db, head):
                if cell["cell_id"] in chosen:
                    fields = {k: {"before": cell[k], "after": v} for k, v in changes.items() if cell[k] != v}
                    if fields:
                        deltas.append({"cell_id": cell["cell_id"], "changes": fields})
            if not deltas:
                raise SpatialError("The proposal changes no selected cells.")
            affected = self._current_run_ids(db, head)
            result = {"proposal_id": _id("proposal"), "base_revision": head, "selection_id": selection_id,
                      "changes": changes, "rationale": rationale, "actor": actor, "created_at": _now(),
                      "selected_cell_count": len(chosen), "changed_cell_count": len(deltas),
                      "counts": {"selected": len(chosen), "changed": len(deltas)}, "deltas": deltas,
                      "affected_result_ids": affected, "affected_run_ids": affected,
                      "impact_policy": "Existing results at this revision become historical and require recomputation.",
                      "scientific_authorization": "NOT_ESTABLISHED"}
            db.execute("INSERT INTO proposals VALUES(?,?)", (result["proposal_id"], _json(result)))
            return result

    @classmethod
    def _current_run_ids(cls, db, revision: str) -> list[str]:
        return [row["id"] for row in db.execute("SELECT id,payload FROM runs ORDER BY rowid")
                if json.loads(row["payload"])["target_revision"] == revision]

    @staticmethod
    def _approval(reviewer: str, confirmation: bool):
        _text(reviewer, "named reviewer", limit=512)
        if confirmation is not True:
            raise SpatialError("Explicit researcher confirmation=true is required after reviewing the preview.")

    def _commit(self, db, head: str, overlays: dict, reviewer: str, rationale: str,
                kind: str, *, proposal_id: str | None = None, target_revision: str | None = None) -> dict:
        _, digest = self._source(db)
        invalidated = self._current_run_ids(db, head)
        revision = {"revision_id": _id("rev"), "parent_revision": head, "overlays": overlays,
                    "created_at": _now(), "kind": kind, "reviewer": reviewer, "confirmation": True,
                    "reviewer_identity": "user_supplied_attribution_not_authenticated",
                    "rationale": rationale, "source_sha256": digest,
                    "scientific_authorization": "NOT_ESTABLISHED"}
        if proposal_id:
            revision["proposal_id"] = proposal_id
        if target_revision:
            revision["reverted_to_revision"] = target_revision
        db.execute("INSERT INTO revisions VALUES(?,?,?,?)", (revision["revision_id"], _json(revision), _hash(revision), proposal_id))
        updated = db.execute("UPDATE state SET value=? WHERE key='head_revision' AND value=?", (revision["revision_id"], head))
        if updated.rowcount != 1:
            raise SpatialError("Stale revision: a concurrent edit won. Refresh and review again.")
        return {**revision, "revision_sha256": _hash(revision), "invalidated_run_ids": invalidated}

    def apply_revision(self, proposal_id: str, expected_revision: str,
                       reviewer: str, confirmation: bool) -> dict:
        self._approval(reviewer, confirmation)
        with self._db(write=True) as db:
            head = self._head(db, expected_revision)
            proposal = self._record(db, "proposals", proposal_id)
            if proposal["base_revision"] != head:
                raise SpatialError("Stale proposal: create and review a new preview.")
            overlays = self._record(db, "revisions", head)["overlays"]
            for delta in proposal["deltas"]:
                overlay = overlays.setdefault(delta["cell_id"], {})
                overlay.update({k: v["after"] for k, v in delta["changes"].items()})
            return self._commit(db, head, overlays, reviewer, proposal["rationale"], "edit", proposal_id=proposal_id)

    def revert_revision(self, target_revision: str, expected_revision: str,
                        reviewer: str, confirmation: bool) -> dict:
        self._approval(reviewer, confirmation)
        with self._db(write=True) as db:
            head = self._head(db, expected_revision)
            target = self._record(db, "revisions", target_revision)
            return self._commit(db, head, target["overlays"], reviewer,
                                f"Restore overlays from {target_revision}", "revert", target_revision=target_revision)

    def save_run(self, result: dict) -> dict:
        if not isinstance(result, dict):
            raise SpatialError("Analysis result must be an object.")
        payload = json.loads(_json(result))
        region_method = {"spatial-collab.region-contrast.v1": "region_comparison", "spatial-collab.protein-region.v1": "protein_region_comparison"}.get(payload.get("analysis_schema"))
        if region_method:
            from .method_contracts import contract
            payload["method_contract"] = contract(region_method)
        with self._db(write=True) as db:
            base = self._record(db, "revisions", payload.get("base_revision"))
            target = self._record(db, "revisions", payload.get("target_revision"))
            sid = payload.get("selection_id")
            if sid is not None:
                self._record(db, "selections", sid)
            source, digest = self._source(db)
            provenance = payload.get("provenance", {})
            if not isinstance(provenance, dict):
                raise SpatialError("provenance must be an object.")
            payload.pop("stale", None)
            payload.pop("applicability", None)
            payload["run_id"] = _id("run")
            payload["created_at"] = _now()
            payload["scientific_authorization"] = "NOT_ESTABLISHED"
            payload["provenance"] = {**provenance, "project_id": self._state(db, "project_id"),
                    "source_sha256": digest, "base_revision_sha256": _hash(base), "target_revision_sha256": _hash(target),
                    "source_kind": source["metadata"]["source_kind"], "identity_scope": "local content integrity; not independent scientific validation"}
            db.execute("INSERT INTO runs VALUES(?,?,?)", (payload["run_id"], _json(payload), _hash(payload)))
            return {**self._run_status(db, payload), "run_sha256": _hash({k: v for k, v in payload.items() if k not in ("stale", "applicability")})}

    def get_run(self, run_id: str) -> dict:
        with self._db() as db:
            result = self._record(db, "runs", run_id)
            digest = _hash(result)
            return {**self._run_status(db, result), "run_sha256": digest}

    def export_bundle(self, run_id: str | None = None, compact: bool = False) -> dict:
        if type(compact) is not bool:
            raise SpatialError("compact must be a boolean.")
        # Read every artifact in one transaction so concurrent edits cannot mix states.
        with self._db() as db:
            source, digest = self._source(db)
            head = self._state(db, "head_revision")
            revisions = [json.loads(r[0]) for r in db.execute("SELECT payload FROM revisions ORDER BY rowid")]
            selection = self._selection(db)
            run = self._record(db, "runs", run_id) if run_id is not None else None
            run_selection = self._selection(db, run["selection_id"]) if run and run.get("selection_id") else None
            exported_at = _now()
            provenance = {"format": "spatial-collab-review-bundle-1", "exported_at": exported_at,
                          "project_id": self._state(db, "project_id"), "head_revision_at_export": head,
                          "source_sha256": digest, "scientific_authorization": "NOT_ESTABLISHED",
                          "revision_hashes": {revision["revision_id"]: _hash(revision) for revision in revisions},
                          "run_sha256": _hash(run) if run else None,
                          "identity_scope": "Reviewer names are user attribution, not authenticated identity.",
                          "evidence_scope": "Local reproducibility and descriptive analysis; no independent biological validation."}
            if compact:
                provenance["materialization"] = "source_plus_overlays"
            files = {"source_snapshot.json": source, "revisions.json": revisions, "selection.json": selection,
                     "provenance.json": provenance}
            if run:
                files["result.json"] = run
                files["result_status_at_export.json"] = {"run_id": run["run_id"], "stale": run["target_revision"] != head,
                                                        "head_revision_at_export": head}
                files["run_selection.json"] = run_selection
                if not compact:
                    files["base_cells.json"] = self._materialize(db, run["base_revision"])
                    files["target_cells.json"] = self._materialize(db, run["target_revision"])
            elif not compact:
                files["current_cells.json"] = self._materialize(db, head)
            files["proposals.json"] = [json.loads(r[0]) for r in db.execute("SELECT payload FROM proposals ORDER BY rowid")]
            files["selections.json"] = [json.loads(r[0]) for r in db.execute("SELECT payload FROM selections ORDER BY rowid")]
        # Check a replayable size budget before creating an incomplete export.
        from .proteomics import _records
        assays = _records(self)
        if assays:
            from .array_assays import export_record
            files["protein_assays.json"] = [export_record(self, a) if a["schema"].endswith(".v2") else a for a in assays]
        from . import objects
        research = [objects.get(self, oid, kind) for kind in ("integration", "correspondence", "measurementlayer", "molecularrelation", "study", "studyresult", "atlas", "pyramid", "assaydescriptor", "integrationexperiment", "integrationexperimentrun", "analysisinput", "powerdesign", "powerplan", "analysisresult") for oid in objects.catalog(self, kind)]
        if research:
            files["research_objects.json"] = research
        contents = {name: (_json(value) + "\n").encode("utf-8") for name, value in files.items()}
        if any(len(value) > 256 * 1024**2 for value in contents.values()) or sum(map(len, contents.values())) > 512 * 1024**2:
            raise SpatialError("Export exceeds replay budget (256 MiB per file / 512 MiB total). Use compact=True to avoid redundant materializations or an explicitly bounded project.")
        exports = self.root / "exports"
        if exports.is_symlink() or not exports.resolve().is_relative_to(self.root):
            raise SpatialError("Export directory must remain inside the project.")
        exports.mkdir(exist_ok=True)
        destination = exports / _id("review")
        destination.mkdir()
        manifest = {}
        for name, content in contents.items():
            (destination / name).write_bytes(content)
            manifest[name] = {"sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)}
        manifest_content = {"algorithm": "sha256", "files": manifest,
                            "scientific_authorization": "NOT_ESTABLISHED"}
        if compact:
            manifest_content["layout"] = "source_plus_overlays"
        (destination / "manifest.json").write_text(_json(manifest_content) + "\n", encoding="utf-8")
        return {"export_path": str(destination), "manifest_path": str(destination / "manifest.json"),
                "manifest": manifest_content, "run_id": run_id, "source_sha256": digest,
                "scientific_authorization": "NOT_ESTABLISHED"}
