"""Disk-backed spatial atlases: streamed rows, exact queries and declared LOD.

An atlas is an immutable, independently identified dataset. A bounded review
project can attach it without loading its full expression matrix into JSON.
"""
from contextlib import contextmanager
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3

from . import objects
from .identity import qualify_observation_id
from .store import SpatialError, _hash, _json, _number, _polygon, _text, _id


def matrix_rows(path, spec, db=None):
    """Read raw CSR h5ad / 10x CSC h5 in bounded observation chunks."""
    import h5py
    def strings(dataset):
        return dataset.asstr()[:]
    with h5py.File(path, "r") as h:
        if spec["format"] == "10x_h5":
            group = h["matrix"]
            cells = group["barcodes"]
            genes = strings(group["features"][spec.get("feature_field", "id")])
        elif spec["format"] == "h5ad_csr":
            layer = spec["counts_layer"]
            group = h["X" if layer == "X" else "layers/" + layer]
            if group.attrs.get("encoding-type") != "csr_matrix":
                raise SpatialError("Streaming h5ad reader requires producer raw CSR counts.")
            cells = h["obs"][h["obs"].attrs["_index"]]
            key = spec.get("feature_field", "_index")
            genes = strings(h["var"][h["var"].attrs["_index"] if key == "_index" else key])
        else:
            raise SpatialError("Unknown matrix encoding.")
        if len(set(genes)) != len(genes):
            raise SpatialError("Matrix features must have unique producer IDs.")
        if len(group["indptr"]) != len(cells)+1 or int(group["indptr"][0]) != 0 or int(group["indptr"][-1]) != len(group["data"]) or len(group["indices"]) != len(group["data"]):
            raise SpatialError("Malformed sparse matrix pointers or data lengths.")
        if db is not None:
            if set(genes) != {r[0] for r in db.execute("SELECT id FROM features")} or len(cells) != db.execute("SELECT count(*) FROM cells").fetchone()[0]:
                raise SpatialError("Matrix axes must exactly equal the declared observations and complete feature panel.")
            db.execute("CREATE TEMP TABLE matrix_axis(source_id TEXT PRIMARY KEY) WITHOUT ROWID")
        for start in range(0, len(cells), 256):
            pointers = group["indptr"][start:min(start+256, len(cells))+1]
            lo, hi = int(pointers[0]), int(pointers[-1])
            if any(b < a for a, b in zip(pointers, pointers[1:])) or lo < 0 or hi > len(group["data"]):
                raise SpatialError("Sparse matrix pointers are invalid.")
            if hi-lo > 5000000:
                raise SpatialError("256-row raw matrix block exceeds 5M nonzero budget.")
            indices, values = group["indices"][lo:hi], group["data"][lo:hi]
            if len(indices) and (indices.min() < 0 or indices.max() >= len(genes)):
                raise SpatialError("Sparse matrix feature index is outside its declared axis.")
            cell_ids = cells.asstr()[start:min(start+256, len(cells))]
            if db is not None:
                db.executemany("INSERT INTO matrix_axis VALUES(?)", ((cid,) for cid in cell_ids))
            for j in range(len(pointers)-1):
                for k in range(int(pointers[j])-lo, int(pointers[j+1])-lo):
                    yield {"cell_id": cell_ids[j], "feature": genes[indices[k]], "value": values[k]}
        if db is not None:
            if db.execute("SELECT source_id FROM matrix_axis EXCEPT SELECT source_id FROM cells LIMIT 1").fetchone():
                raise SpatialError("Matrix contains unknown observation IDs, including zero-library rows.")
            db.execute("DROP TABLE matrix_axis")


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def rows(path):
    """Bounded CSV/gzip or Parquet batches; no dataframe of all observations."""
    path = Path(path)
    if path.suffix.lower() == ".parquet":
        import pyarrow.parquet as pq
        for batch in pq.ParquetFile(path).iter_batches(batch_size=4096):
            yield from batch.to_pylist()
    else:
        opener = gzip.open if path.suffix.lower() == ".gz" else open
        with opener(path, "rt", encoding="utf-8-sig", newline="") as stream:
            yield from csv.DictReader(stream)


def _float(value, field):
    if isinstance(value, bool):
        raise SpatialError(f"{field} cannot be boolean.")
    try:
        return _number(float(value), field)
    except (TypeError, ValueError, OverflowError) as exc:
        raise SpatialError(f"{field} must be finite numeric data.") from exc


def _bound(value):
    if not isinstance(value, list) or len(value) != 4:
        raise SpatialError("bounds requires [xmin,ymin,xmax,ymax].")
    b = [_number(x, "bounds") for x in value]
    if b[0] >= b[2] or b[1] >= b[3]:
        raise SpatialError("Viewport must have positive width and height.")
    return b


def build(project, spec):
    """Local-only import. Explicit column maps; foreign-key errors abort import."""
    meta = spec.get("metadata", {})
    for k in ("sample_id", "name", "coordinate_system", "platform"):
        _text(meta.get(k), k)
    if meta.get("units") not in {"micrometer", "pixel", "array_index"}:
        raise SpatialError("Atlas requires explicit coordinate units.")
    if meta.get("observation_unit") not in {"cell", "spot", "bin"}:
        raise SpatialError("Declare cell, spot or bin; bins are not pure cells.")
    if meta.get("species") not in {"human", "mouse"}:
        raise SpatialError("Declare human or mouse species.")
    sources = {}
    for kind in ("cells", "features", "counts", "transcripts", "boundaries"):
        if kind in spec:
            path = Path(spec[kind]["path"]).resolve(strict=True)
            sources[kind] = {"path": str(path), "sha256": digest(path)}
    if not {"cells", "features", "counts"} <= set(sources):
        raise SpatialError("Cells, complete feature panel and raw counts are required.")
    folder = project.root / "atlases"
    folder.mkdir(exist_ok=True)
    if not folder.resolve().is_relative_to(project.root):
        raise SpatialError("Atlas storage must remain inside project.")
    path = folder / (_id("atlasdata") + ".sqlite3")
    with path.open("xb"):
        pass
    db = sqlite3.connect(path)
    try:
        db.execute("PRAGMA foreign_keys=ON")
        db.execute("PRAGMA temp_store=FILE")
        db.executescript("""
        CREATE TABLE features(id TEXT PRIMARY KEY,symbol TEXT);
        CREATE TABLE cells(id INTEGER PRIMARY KEY,cell_id TEXT UNIQUE NOT NULL,source_id TEXT UNIQUE NOT NULL,x REAL,y REAL,label TEXT);
        CREATE VIRTUAL TABLE cell_rtree USING rtree(id,x0,x1,y0,y1);
        CREATE TABLE counts(cell INTEGER REFERENCES cells(id),feature TEXT REFERENCES features(id),value REAL,PRIMARY KEY(cell,feature)) WITHOUT ROWID;
        CREATE TABLE transcripts(id INTEGER PRIMARY KEY,source_id TEXT UNIQUE,feature TEXT,x REAL,y REAL,z REAL,qv REAL,cell_id TEXT);
        CREATE VIRTUAL TABLE transcript_rtree USING rtree(id,x0,x1,y0,y1);
        CREATE TABLE vertices(boundary_id TEXT,cell_id TEXT,kind TEXT,vertex INTEGER,x REAL,y REAL,PRIMARY KEY(boundary_id,vertex)) WITHOUT ROWID;
        CREATE TABLE boundaries(id INTEGER PRIMARY KEY,boundary_id TEXT UNIQUE,cell_id TEXT,kind TEXT,polygon TEXT);
        CREATE VIRTUAL TABLE boundary_rtree USING rtree(id,x0,x1,y0,y1);
        CREATE TABLE lod(level INTEGER,gx INTEGER,gy INTEGER,n INTEGER,PRIMARY KEY(level,gx,gy)) WITHOUT ROWID;
        CREATE TABLE transcript_lod(level INTEGER,gx INTEGER,gy INTEGER,n INTEGER,PRIMARY KEY(level,gx,gy)) WITHOUT ROWID;
        """)
        def mapped(kind):
            if kind == "counts" and spec[kind].get("format") in {"10x_h5", "h5ad_csr"}:
                yield from matrix_rows(sources[kind]["path"], spec[kind], db)
                return
            mapping = spec[kind]["columns"]
            for row in rows(sources[kind]["path"]):
                if any(col not in row for col in mapping.values()):
                    raise SpatialError(f"Missing declared {kind} column.")
                yield {**{k: row[col] for k, col in mapping.items()}, **spec[kind].get("constants", {})}
        for row in mapped("features"):
            db.execute("INSERT INTO features VALUES(?,?)", (_text(row["id"], "feature"), row.get("symbol") or None))
        for i, row in enumerate(mapped("cells"), 1):
            cid = _text(str(row["id"]), "cell ID", limit=500)
            x, y = _float(row["x"], "x"), _float(row["y"], "y")
            if abs(x) > 1e12 or abs(y) > 1e12:
                raise SpatialError("Coordinates exceed stable spatial index range.")
            db.execute("INSERT INTO cells VALUES(?,?,?,?,?,?)", (i, qualify_observation_id(meta["sample_id"], cid), cid, x, y, _text(row.get("label", "Unannotated"), "label")))
            db.execute("INSERT INTO cell_rtree VALUES(?,?,?,?,?)", (i, x, x, y, y))
        n = db.execute("SELECT count(*) FROM cells").fetchone()[0]
        nf = db.execute("SELECT count(*) FROM features").fetchone()[0]
        if not n or not nf:
            raise SpatialError("Empty observation or feature axis.")
        for row in mapped("counts"):
            value = _float(row["value"], "raw count")
            if value < 0 or value > 2**53-1 or value != int(value):
                raise SpatialError("Counts must be nonnegative exact integers.")
            cell = db.execute("SELECT id FROM cells WHERE source_id=?", (str(row["cell_id"]),)).fetchone()
            if not cell:
                raise SpatialError("Count references unknown cell ID.")
            db.execute("INSERT INTO counts VALUES(?,?,?)", (cell[0], row["feature"], value))
        db.execute("CREATE INDEX counts_feature ON counts(feature,cell)")
        if "transcripts" in sources:
            if spec["transcripts"].get("projection") != "xy_projection":
                raise SpatialError("Explicit xy_projection required; z is retained but not rendered as 3D.")
            for i, row in enumerate(mapped("transcripts"), 1):
                x, y = _float(row["x"], "x"), _float(row["y"], "y")
                fid = _text(row["feature"], "transcript feature")
                # Controls and unassigned molecules remain visible; never relabel
                # a transcript's assignment as measured cell identity.
                db.execute("INSERT INTO transcripts VALUES(?,?,?,?,?,?,?,?)", (i, _text(str(row["id"]), "transcript ID"), fid, x, y, _float(row["z"], "z") if row.get("z") not in (None, "") else None, _float(row["qv"], "qv") if row.get("qv") not in (None, "") else None, str(row["cell_id"]) if row.get("cell_id") is not None else None))
                db.execute("INSERT INTO transcript_rtree VALUES(?,?,?,?,?)", (i, x, x, y, y))
            db.execute("CREATE INDEX transcript_feature ON transcripts(feature)")
        if "boundaries" in sources:
            for row in mapped("boundaries"):
                if row["kind"] not in {"cell", "nucleus"}:
                    raise SpatialError("Boundary kind must explicitly be cell or nucleus.")
                cid = str(row["cell_id"])
                if not db.execute("SELECT 1 FROM cells WHERE source_id=?", (cid,)).fetchone():
                    raise SpatialError("Boundary references unknown cell.")
                vi = _float(row["vertex"], "vertex index")
                if vi < 0 or vi != int(vi):
                    raise SpatialError("Vertex order requires nonnegative integers.")
                db.execute("INSERT INTO vertices VALUES(?,?,?,?,?,?)", (_text(str(row["id"]), "boundary ID"), cid, row["kind"], int(vi), _float(row["x"], "x"), _float(row["y"], "y")))
            for i, (bid,) in enumerate(db.execute("SELECT DISTINCT boundary_id FROM vertices"), 1):
                vertices = db.execute("SELECT * FROM vertices WHERE boundary_id=? ORDER BY vertex LIMIT 1002", (bid,)).fetchall()
                if len({(v[1], v[2]) for v in vertices}) != 1:
                    raise SpatialError("Boundary ID mixes objects or segmentation kinds.")
                poly = _polygon([[v[4], v[5]] for v in vertices])
                db.execute("INSERT INTO boundaries VALUES(?,?,?,?,?)", (i, bid, vertices[0][1], vertices[0][2], _json(poly)))
                xs, ys = zip(*poly)
                db.execute("INSERT INTO boundary_rtree VALUES(?,?,?,?,?)", (i, min(xs), max(xs), min(ys), max(ys)))
        b = list(db.execute("SELECT min(x),min(y),max(x),max(y) FROM cells").fetchone())
        # A molecule or segmentation may lie beyond its cell centroid. The
        # whole-dataset extent must not silently discard those observations.
        for extra in (db.execute("SELECT min(x),min(y),max(x),max(y) FROM transcripts").fetchone(), db.execute("SELECT min(x0),min(y0),max(x1),max(y1) FROM boundary_rtree").fetchone()):
            if extra[0] is not None:
                b = [min(b[0], extra[0]), min(b[1], extra[1]), max(b[2], extra[2]), max(b[3], extra[3])]
        if b[2] == b[0]:
            b[2] += 1
        if b[3] == b[1]:
            b[3] += 1
        for level in (16, 32, 64, 128, 256):
            db.execute("INSERT INTO lod SELECT ?,min(?,cast((x-?)/?*? AS INTEGER)),min(?,cast((y-?)/?*? AS INTEGER)),count(*) FROM cells GROUP BY 2,3", (level, level-1, b[0], b[2]-b[0], level, level-1, b[1], b[3]-b[1], level))
            db.execute("INSERT INTO transcript_lod SELECT ?,min(?,cast((x-?)/?*? AS INTEGER)),min(?,cast((y-?)/?*? AS INTEGER)),count(*) FROM transcripts GROUP BY 2,3", (level, level-1, b[0], b[2]-b[0], level, level-1, b[1], b[3]-b[1], level))
        totals = {t: db.execute(f"SELECT count(*) FROM {t}").fetchone()[0] for t in ("cells", "features", "counts", "transcripts", "boundaries")}
        for item in sources.values():
            if digest(item["path"]) != item["sha256"]:
                raise SpatialError("Source changed during atlas import.")
        db.commit()
        db.execute("PRAGMA optimize")
    except Exception as exc:
        db.close()
        path.unlink(missing_ok=True)
        if isinstance(exc, sqlite3.IntegrityError):
            raise SpatialError("Duplicate ID/count pair or feature/observation foreign-key mismatch; atlas import aborted.") from exc
        raise
    finally:
        db.close()
    stat = path.stat()
    return objects.put(project, "atlas", {"schema": "spatial-collab.atlas.v1", "metadata": meta, "sources": sources, "column_contract": {k: spec[k] for k in sources}, "bounds": b, "extent_scope": "all_layers", "totals": totals, "database": path.relative_to(project.root).as_posix(), "database_sha256": digest(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns, "lod_levels": [16, 32, 64, 128, 256], "storage": "streamed_sqlite_rtree_sparse_counts", "scientific_authorization": "NOT_ESTABLISHED"})


@contextmanager
def connect(project, atlas_id, verify=False):
    record = objects.get(project, atlas_id, "atlas")
    path = (project.root / record["database"]).resolve(strict=True)
    if not path.is_relative_to(project.root):
        raise SpatialError("Atlas path escaped project.")
    stat = path.stat()
    if (stat.st_size, stat.st_mtime_ns) != (record["size"], record["mtime_ns"]) or verify and digest(path) != record["database_sha256"]:
        raise SpatialError("Atlas storage changed; reimport or restore immutable dataset.")
    db = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)
    db.row_factory = sqlite3.Row
    try:
        yield db, record
    finally:
        db.close()


def catalog(project):
    return {"atlases": [objects.get(project, oid, "atlas") for oid in objects.catalog(project, "atlas")]}


def share_view(project, atlas_id, bounds, layer="cells", feature=None, image_id=None):
    if layer not in {"cells", "transcripts"}:
        raise SpatialError("Choose cells or transcripts for shared viewport context.")
    record = objects.get(project, atlas_id, "atlas")
    if feature is not None:
        _text(feature, "feature")
    if image_id and objects.get(project, image_id, "pyramid")["atlas_id"] != atlas_id:
        raise SpatialError("Image is registered to a different atlas.")
    context = {"atlas_id": atlas_id, "atlas_sha256": record["object_sha256"], "sample_id": record["metadata"]["sample_id"], "coordinate_system": record["metadata"]["coordinate_system"], "units": record["metadata"]["units"], "bounds": _bound(bounds), "layer": layer, "feature": feature, "image_id": image_id, "selection_semantics": "viewport only; not an annotation or individual-cell selection"}
    context["context_id"] = _hash(context)
    with project._db(write=True) as db:
        db.execute("INSERT INTO state(key,value) VALUES('atlas_view',?) ON CONFLICT(key) DO UPDATE SET value=excluded.value", (_json(context),))
    return context


def _query(db, kind, b, feature=None):
    table, tree = ("cells", "cell_rtree") if kind == "cells" else ("transcripts", "transcript_rtree")
    where = "r.x1>=? AND r.x0<=? AND r.y1>=? AND r.y0<=? AND t.x>=? AND t.x<=? AND t.y>=? AND t.y<=?"
    args = [b[0], b[2], b[1], b[3]] * 2
    if feature is not None:
        if kind != "transcripts":
            raise SpatialError("Feature filter applies to transcript molecules.")
        # Selective feature lookup must not scan a whole-slide R-tree first.
        return "FROM transcripts t INDEXED BY transcript_feature WHERE t.feature=? AND t.x>=? AND t.x<=? AND t.y>=? AND t.y<=?", [_text(feature, "feature"), b[0], b[2], b[1], b[3]]
    return f"FROM {tree} r JOIN {table} t ON t.id=r.id WHERE {where}", args


def view(project, atlas_id, bounds=None, layer="cells", feature=None, limit=10000, grid=64):
    if layer not in {"cells", "transcripts", "boundaries"} or type(limit) is not int or not 1 <= limit <= 10000 or grid not in (16, 32, 64, 128):
        raise SpatialError("Invalid atlas layer, point budget or LOD grid.")
    with connect(project, atlas_id) as (db, record):
        b = _bound(bounds or record["bounds"])
        if layer == "boundaries":
            sql = "FROM boundary_rtree r JOIN boundaries t ON r.id=t.id WHERE r.x1>=? AND r.x0<=? AND r.y1>=? AND r.y0<=?"
            args = [b[0], b[2], b[1], b[3]]
            total = db.execute("SELECT count(*) " + sql, args).fetchone()[0]
            limit = min(limit, 2000)
            result = [dict(r) for r in db.execute("SELECT t.* " + sql + " ORDER BY t.id LIMIT ?", [*args, limit])] if total <= limit else []
            vertices = sum(len(json.loads(r["polygon"])) for r in result)
            if vertices > 100000:
                result = []
            for row in result:
                row["polygon"] = json.loads(row["polygon"])
            return {"atlas_id": atlas_id, "layer": layer, "bounds": b, "total": total, "mode": "exact_boundaries" if result or not total else "zoom_required", "records": result, "complete": len(result) == total, "vertex_budget": 100000}
        sql, args = _query(db, layer, b, feature)
        whole = b == record["bounds"] and feature is None and (layer == "cells" or record.get("extent_scope") == "all_layers")
        total = record["totals"][layer] if whole else db.execute("SELECT count(*) " + sql, args).fetchone()[0]
        if total <= limit:
            records = [dict(r) for r in db.execute("SELECT t.* " + sql + " ORDER BY t.id", args)]
            return {"atlas_id": atlas_id, "layer": layer, "bounds": b, "total": total, "mode": "exact_points", "records": records, "complete": True, "feature": feature}
        dx, dy = (b[2]-b[0])/grid, (b[3]-b[1])/grid
        lod_table = "lod" if layer == "cells" else "transcript_lod"
        if whole and db.execute("SELECT 1 FROM sqlite_master WHERE name=?", (lod_table,)).fetchone():
            bins = db.execute(f"SELECT gx,gy,n FROM {lod_table} WHERE level=?", (grid,)).fetchall()
        else:
            bins = db.execute("SELECT min(?,cast((t.x-?)/? AS INTEGER)) gx,min(?,cast((t.y-?)/? AS INTEGER)) gy,count(*) n " + sql + " GROUP BY gx,gy", [grid-1, b[0], dx, grid-1, b[1], dy, *args]).fetchall()
        aggregate = [{"x": b[0]+r[0]*dx, "y": b[1]+r[1]*dy, "width": dx, "height": dy, "count": r[2]} for r in bins]
        return {"atlas_id": atlas_id, "layer": layer, "bounds": b, "total": total, "mode": "aggregate_density", "records": aggregate, "complete": True, "aggregate_total": sum(r["count"] for r in aggregate), "feature": feature, "notice": "All matching objects counted; bins are density aggregates, not cells or a sampled point cloud."}


def extract(project, atlas_id, bounds, max_cells=10000):
    """Freeze exact ROI counts into the executable analysis input contract."""
    from scipy import sparse
    from .workflow_inputs import register_counts
    if type(max_cells) is not int or not 1 <= max_cells <= 100000:
        raise SpatialError("Choose a bounded analysis ROI of at most 100000 objects.")
    b = _bound(bounds)
    with connect(project, atlas_id, verify=True) as (db, record):
        sql, args = _query(db, "cells", b)
        n = db.execute("SELECT count(*) " + sql, args).fetchone()[0]
        if not 0 < n <= max_cells:
            raise SpatialError(f"ROI contains {n} objects; limit is {max_cells}. No sampling performed.")
        cells = db.execute("SELECT t.* " + sql + " ORDER BY t.id", args).fetchall()
        features = db.execute("SELECT * FROM features ORDER BY id").fetchall()
        feature_index = {r[0]: i for i, r in enumerate(features)}
        data, indices, indptr = [], [], [0]
        for c in cells:
            for f, value in db.execute("SELECT feature,value FROM counts WHERE cell=? ORDER BY feature", (c["id"],)):
                data.append(value)
                indices.append(feature_index[f])
                if len(data) > 25000000:
                    raise SpatialError("ROI nonzero budget exceeded; narrow the ROI.")
            indptr.append(len(data))
        matrix = sparse.csr_matrix((data, indices, indptr), shape=(n, len(features)))
        return register_counts(project, matrix, [c["cell_id"] for c in cells], [f[0] for f in features], feature_symbols=[f[1] for f in features], sample_id=record["metadata"]["sample_id"], species=record["metadata"]["species"], kind="spatial", coordinates=[[c["x"], c["y"]] for c in cells], labels=[c["label"] for c in cells], provenance={"atlas_id": atlas_id, "atlas_sha256": record["object_sha256"], "database_sha256": record["database_sha256"], "bounds": b, "scope": "exact_atlas_roi", "units": record["metadata"]["units"], "coordinate_system": record["metadata"]["coordinate_system"]})
