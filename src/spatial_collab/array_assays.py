"""Chunked immutable protein arrays; metadata stays JSON, pages are memory mapped."""
from collections.abc import Mapping
import hashlib

import numpy as np

from .store import SpatialError, _hash, _json

SCHEMA = "spatial-collab.protein-assay.v2"


class ChunkValues(Mapping):
    def __init__(self, folder, chunks):
        self.folder, self.chunks = folder, chunks
        self.index = {cid: (i, j) for i, chunk in enumerate(chunks) for j, cid in enumerate(chunk["ids"])}
        self.loaded = {}

    def __len__(self):
        return len(self.index)

    def __iter__(self):
        return iter(self.index)

    def __getitem__(self, cid):
        i, j = self.index[cid]
        if i not in self.loaded:
            self.loaded = {i: np.load(self.folder / self.chunks[i]["file"], mmap_mode="r", allow_pickle=False)}
        return [None if np.isnan(x) else float(x) for x in self.loaded[i][j]]


def save(project, metadata, ids, matrix):
    folder = project.root / "assays"
    if not folder.resolve().is_relative_to(project.root):
        raise SpatialError("Assay arrays must remain inside the project.")
    folder.mkdir(exist_ok=True)
    if len(list(folder.glob("protein_*.json"))) >= 16:
        raise SpatialError("At most 16 protein assays per project.")
    if matrix.shape != (len(ids), len(metadata["features"])):
        raise SpatialError("Protein matrix does not match feature/observation axes.")
    chunks, missing = [], 0
    for offset in range(0, len(ids), 2048):
        block = matrix[offset:offset + 2048]
        block = block.toarray() if hasattr(block, "toarray") else np.asarray(block)
        block = np.asarray(block, dtype=np.float64)
        if np.isinf(block).any() or (block < 0).any() or metadata["measurement_type"] == "antibody_count" and np.any(block[np.isfinite(block)] % 1):
            raise SpatialError("Protein arrays require nonnegative values or missing NaN, and integer antibody counts.")
        missing += int(np.isnan(block).sum())
        import io
        buf = io.BytesIO()
        np.save(buf, block, allow_pickle=False)
        content = buf.getvalue()
        digest = hashlib.sha256(content).hexdigest()
        filename = "matrix_" + digest + ".npy"
        path = folder / filename
        if path.exists():
            if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
                raise SpatialError("Existing array chunk failed integrity.")
        else:
            with path.open("xb") as stream:
                stream.write(content)
        chunks.append({"file": filename, "sha256": digest, "shape": list(block.shape), "ids": ids[offset:offset + len(block)]})
    record = {**metadata, "schema": SCHEMA, "missing_measurements": missing, "chunks": chunks,
              "storage": "chunked_npy_float64", "observation_count": len(ids)}
    record["assay_id"] = "protein_" + _hash(record)[:24]
    record["assay_sha256"] = _hash(record)
    path = folder / (record["assay_id"] + ".json")
    with path.open("x", encoding="utf-8") as stream:
        stream.write(_json(record) + "\n")
    return record


def validate(project, record):
    # Exported snapshots carry the same chunks as lists, preserving chunk checksums.
    payload = {k: v for k, v in record.items() if k not in {"assay_sha256", "values", "embedded_chunks"}}
    if _hash(payload) != record.get("assay_sha256") or record.get("assay_id") != "protein_" + _hash({k: v for k, v in payload.items() if k != "assay_id"})[:24]:
        raise SpatialError("Protein array metadata integrity failed.")
    context = project.context() if hasattr(project, "context") else project.summary()
    if any(record.get(k) != context[k] for k in ("project_id", "source_sha256")):
        raise SpatialError("Protein array belongs to another project/source.")
    chunks = record["chunks"]
    all_ids = [cid for chunk in chunks for cid in chunk["ids"]]
    if len(all_ids) != len(set(all_ids)) or len(all_ids) != record["observation_count"]:
        raise SpatialError("Protein array ID axis mismatch.")
    embedded = record.get("embedded_chunks")
    missing = 0
    values = {}
    for i, chunk in enumerate(chunks):
        if embedded is not None:
            import io
            array = np.asarray(embedded[i], dtype=float)
            buf = io.BytesIO()
            np.save(buf, array, allow_pickle=False)
            content = buf.getvalue()
        else:
            path = project.root / "assays" / chunk["file"]
            if not path.resolve().is_relative_to(project.root / "assays") or path.stat().st_size > 16 * 1024**2:
                raise SpatialError("Invalid array path/size.")
            content = path.read_bytes()
            array = np.load(path, mmap_mode="r", allow_pickle=False)
        if hashlib.sha256(content).hexdigest() != chunk["sha256"] or list(array.shape) != chunk["shape"] or array.shape != (len(chunk["ids"]), len(record["features"])):
            raise SpatialError("Protein array chunk integrity failed.")
        if np.isinf(array).any() or (array < 0).any():
            raise SpatialError("Invalid protein array measurement.")
        missing += int(np.isnan(array).sum())
        if embedded is not None:
            values.update({cid: [None if np.isnan(x) else float(x) for x in row] for cid, row in zip(chunk["ids"], array)})
    if missing != record["missing_measurements"]:
        raise SpatialError("Protein array missing-value count mismatch.")
    return {**record, "values": values if embedded is not None else ChunkValues(project.root / "assays", chunks)}


def export_record(project, record):
    result = {k: v for k, v in record.items() if k != "values"}
    result["embedded_chunks"] = []
    for chunk in record["chunks"]:
        a = np.load(project.root / "assays" / chunk["file"], allow_pickle=False)
        result["embedded_chunks"].append([[None if np.isnan(x) else float(x) for x in row] for row in a])
    return result
