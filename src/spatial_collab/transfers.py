"""Bounded, project-owned file transfers. Never accept client filesystem paths or URLs."""
from __future__ import annotations

import base64
import binascii
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import re
import shutil
import sqlite3
import time
import uuid
import zipfile

from .store import Project, SpatialError

CHUNK_BYTES = 256 * 1024
UPLOAD_BYTES = 128 * 1024**2
STORAGE_BYTES = 512 * 1024**2
TTL_SECONDS = 24 * 3600


class Transfers:
    def __init__(self, project, *, create=True):
        self.project = project
        self.root = project.root / "transfers"
        if create:
            self.root.mkdir(exist_ok=True)
        elif not self.root.is_dir():
            raise SpatialError("No transfers exist in this project.")
        if not self.root.resolve().is_relative_to(project.root):
            raise SpatialError("Transfer storage must remain inside the project.")
        if create:
            with self.db() as db:
                db.execute("CREATE TABLE IF NOT EXISTS files (id TEXT PRIMARY KEY, name TEXT, size INTEGER, sha TEXT, offset INTEGER, state TEXT, expires REAL)")

    @contextmanager
    def db(self):
        db = sqlite3.connect(self.root / "transfers.sqlite3", timeout=30)
        try:
            with db:
                yield db
        finally:
            db.close()

    def path(self, file_id):
        if not isinstance(file_id, str) or not re.fullmatch(r"[a-f0-9]{32}", file_id):
            raise SpatialError("Invalid transfer ID.")
        path = self.root / file_id
        if not path.resolve().is_relative_to(self.root.resolve()) or path.is_symlink():
            raise SpatialError("Invalid transfer storage.")
        return path

    def record(self, db, file_id):
        self.path(file_id)
        row = db.execute("SELECT * FROM files WHERE id=?", (file_id,)).fetchone()
        if not row or row[6] <= time.time():
            raise SpatialError("Transfer absent or expired; start a new transfer.")
        return row

    def cleanup(self):
        with self.db() as db:
            db.execute("BEGIN IMMEDIATE")
            for (file_id,) in db.execute("SELECT id FROM files WHERE expires<=?", (time.time(),)):
                self.path(file_id).unlink(missing_ok=True)
            db.execute("DELETE FROM files WHERE expires<=?", (time.time(),))

    def begin(self, filename, size, sha256):
        self.cleanup()
        if not isinstance(filename, str) or len(filename) > 120 or not re.fullmatch(r"[\w .-]+\.(json|h5ad)", filename):
            raise SpatialError("Upload a plain .json snapshot or .h5ad filename, without directories.")
        if type(size) is not int or not 0 < size <= UPLOAD_BYTES:
            raise SpatialError("Upload must be 1 byte to 128 MiB. Use local import for larger datasets.")
        if not re.fullmatch(r"[a-f0-9]{64}", sha256):
            raise SpatialError("A SHA256 of the complete file is required.")
        file_id = uuid.uuid4().hex
        with self.db() as db:
            db.execute("BEGIN IMMEDIATE")
            used = db.execute("SELECT COALESCE(sum(size),0) FROM files").fetchone()[0]
            if used + size > STORAGE_BYTES:
                raise SpatialError("Transfer quota exceeded. Remove completed transfers or wait for expiry.")
            self.path(file_id).touch(exist_ok=False)
            db.execute("INSERT INTO files VALUES(?,?,?,?,?,?,?)", (file_id, filename, size, sha256, 0, "uploading", time.time()+TTL_SECONDS))
        return {"file_id": file_id, "chunk_bytes": CHUNK_BYTES, "expires_in_seconds": TTL_SECONDS}

    def status(self, file_id):
        with self.db() as db:
            row = self.record(db, file_id)
        return dict(zip(("file_id", "filename", "size", "sha256", "offset", "state", "expires_at"), row))

    def append(self, file_id, offset, data_base64):
        if not isinstance(data_base64, str) or len(data_base64) > (CHUNK_BYTES + 2)//3*4:
            raise SpatialError("Transfer chunk exceeds 256 KiB.")
        try:
            data = base64.b64decode(data_base64, validate=True)
        except (ValueError, binascii.Error):
            raise SpatialError("Invalid base64 chunk.") from None
        if type(offset) is not int or offset < 0 or not data or len(data) > CHUNK_BYTES:
            raise SpatialError("Invalid chunk offset or size.")
        with self.db() as db:
            db.execute("BEGIN IMMEDIATE")
            row = self.record(db, file_id)
            if row[5] != "uploading" or offset + len(data) > row[2]:
                raise SpatialError("Upload state or declared size mismatch.")
            with self.path(file_id).open("r+b") as stream:
                if offset < row[4]:
                    stream.seek(offset)
                    if offset + len(data) > row[4] or stream.read(len(data)) != data:
                        raise SpatialError("Conflicting retry chunk.")
                elif offset == row[4]:
                    stream.seek(offset)
                    stream.write(data)
                    stream.truncate()
                    db.execute("UPDATE files SET offset=? WHERE id=?", (offset+len(data), file_id))
                else:
                    raise SpatialError("Upload chunks must be sequential; resume from the saved offset.")
        return self.status(file_id)

    def complete(self, file_id):
        with self.db() as db:
            db.execute("BEGIN IMMEDIATE")
            row = self.record(db, file_id)
            if row[5] not in {"uploading", "ready"} or row[4] != row[2]:
                raise SpatialError("Upload is incomplete.")
            with self.path(file_id).open("rb") as stream:
                digest = hashlib.file_digest(stream, "sha256").hexdigest()
            if digest != row[3]:
                raise SpatialError("Upload checksum mismatch. No data was imported.")
            db.execute("UPDATE files SET state='ready' WHERE id=?", (file_id,))
        return self.status(file_id)

    def discard(self, file_id):
        with self.db() as db:
            db.execute("BEGIN IMMEDIATE")
            self.record(db, file_id)
            self.path(file_id).unlink(missing_ok=True)
            db.execute("DELETE FROM files WHERE id=?", (file_id,))
        return {"file_id": file_id, "deleted": True, "imported_projects_unchanged": True}

    def export(self, run_id=None, compact=True):
        self.cleanup()
        exported = self.project.export_bundle(run_id, compact=compact)
        folder = Path(exported["export_path"]).resolve()
        if not folder.is_relative_to(self.project.root / "exports"):
            raise SpatialError("Export escaped project storage.")
        file_id = uuid.uuid4().hex
        path = self.path(file_id)
        try:
            with zipfile.ZipFile(path, "x", zipfile.ZIP_DEFLATED) as archive:
                for item in sorted(folder.rglob("*")):
                    if item.is_symlink() or not item.resolve().is_relative_to(folder):
                        raise SpatialError("Export contains an external link.")
                    if item.is_file():
                        archive.write(item, item.relative_to(folder).as_posix())
            size = path.stat().st_size
            with path.open("rb") as stream:
                digest = hashlib.file_digest(stream, "sha256").hexdigest()
            with self.db() as db:
                db.execute("BEGIN IMMEDIATE")
                if size + db.execute("SELECT COALESCE(sum(size),0) FROM files").fetchone()[0] > STORAGE_BYTES:
                    raise SpatialError("Download exceeds transfer quota; original review bundle is retained locally.")
                db.execute("INSERT INTO files VALUES(?,?,?,?,?,?,?)", (file_id, "spatial-review.zip", size, digest, size, "download", time.time()+TTL_SECONDS))
        except Exception:
            path.unlink(missing_ok=True)
            raise
        return self.status(file_id)

    def read(self, file_id, offset=0, length=CHUNK_BYTES):
        with self.db() as db:
            row = self.record(db, file_id)
            if row[5] != "download" or type(offset) is not int or not 0 <= offset <= row[2] or type(length) is not int or not 1 <= length <= CHUNK_BYTES:
                raise SpatialError("Invalid download state or range.")
            with self.path(file_id).open("rb") as stream:
                stream.seek(offset)
                data = stream.read(length)
        return {"file_id": file_id, "offset": offset, "next_offset": offset+len(data), "size": row[2],
                "sha256": row[3], "data_base64": base64.b64encode(data).decode(), "complete": offset+len(data) == row[2]}

    def import_project(self, file_id, options):
        if self.project.root.parent.name == "imported-projects":
            raise SpatialError("Import new files through the primary connection project (project_id=null).")
        record = self.complete(file_id)
        parent = self.project.root / "imported-projects"
        parent.mkdir(exist_ok=True)
        if not parent.resolve().is_relative_to(self.project.root):
            raise SpatialError("Imported projects must remain inside the connection workspace.")
        # The source handle owns the destination. Repeated import never replaces data.
        target = parent / file_id
        if target.exists():
            raise SpatialError("This file already has an import attempt; upload a new copy to retry.")
        source = self.path(file_id)
        if record["filename"].endswith(".json"):
            if options:
                raise SpatialError("JSON snapshot carries its own metadata; options must be empty.")
            if record["size"] > 32 * 1024**2:
                raise SpatialError("JSON snapshot limit is 32 MiB; use H5AD or local import.")
            data = json.loads(source.read_text(encoding="utf-8"))
            if not isinstance(data, dict) or set(data) != {"cells", "metadata"}:
                raise SpatialError("JSON requires exactly cells and metadata, with explicit source and measurement semantics.")
            metadata = dict(data["metadata"])
            metadata["source_files"] = [{"path": str(target / "sources" / record["filename"]), "sha256": record["sha256"]}]
            project = Project.create(target, data["cells"], metadata)
        else:
            allowed = {"label_key", "spatial_key", "counts_layer", "slice_id", "slice_key", "coordinate_system", "units", "name", "biological_replicates", "observation_unit", "platform", "sample_id", "feature_id_key", "feature_symbol_key", "allow_unannotated"}
            if not isinstance(options, dict) or set(options) - allowed:
                raise SpatialError("Unknown H5AD import options; paths and external resources are not allowed.")
            for key in ("counts_layer", "slice_id", "coordinate_system", "units", "observation_unit"):
                if key not in options:
                    raise SpatialError(f"Explicit {key} is required for H5AD import.")
            check_h5ad(source)
            from .importers import import_h5ad
            originals = target / "sources"
            originals.mkdir(parents=True)
            retained = originals / record["filename"]
            shutil.copyfile(source, retained)
            project = import_h5ad(retained, target, **options)
        if record["filename"].endswith(".json"):
            (target / "sources").mkdir()
            shutil.copyfile(source, target / "sources" / record["filename"])
        return {"project_id": project.summary()["project_id"], "summary": project.summary(),
                "source_sha256": record["sha256"], "primary_project_unchanged": True,
                "next": "Call open_project with this project_id. Include project_id on every subsequent tool call."}


def check_h5ad(path):
    """Reject external/soft links, virtual data and oversized decoded arrays before AnnData reads."""
    import h5py
    with h5py.File(path, "r") as handle:
        pending, seen, total = [handle], set(), 0
        while pending:
            group = pending.pop()
            key = hash(group.id)
            if key in seen:
                continue
            seen.add(key)
            if len(seen) > 10000:
                raise SpatialError("H5AD group budget exceeded.")
            for name in group:
                if not isinstance(group.get(name, getlink=True), h5py.HardLink):
                    raise SpatialError("H5AD external and soft links are not accepted.")
                item = group[name]
                if isinstance(item, h5py.Group):
                    pending.append(item)
                else:
                    if item.is_virtual or item.external:
                        raise SpatialError("H5AD external/virtual datasets are not accepted.")
                    total += item.size * max(item.dtype.itemsize, 8)
                    if item.size > 10_000_000 or total > 512 * 1024**2:
                        raise SpatialError("H5AD decoded array budget exceeded; import a declared subset locally.")


def projects(root_project):
    result = [root_project]
    parent = root_project.root / "imported-projects"
    if parent.is_dir() and parent.resolve().is_relative_to(root_project.root):
        for child in sorted(parent.iterdir()):
            if not child.is_symlink() and child.resolve().is_relative_to(parent.resolve()) and (child / "project.sqlite3").is_file():
                try:
                    result.append(Project(child))
                except SpatialError:
                    pass
    return result
