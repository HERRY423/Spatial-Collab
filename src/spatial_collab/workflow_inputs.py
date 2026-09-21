"""Immutable, portable sparse count inputs for executable scientific workflows."""

import base64
import hashlib
import io
from pathlib import Path
import zipfile

import numpy as np
from scipy import sparse

from . import objects
from .store import SpatialError, _text


def register_counts(
    project,
    matrix,
    observation_ids,
    features,
    *,
    sample_id,
    species,
    kind,
    provenance,
    coordinates=None,
    labels=None,
    batches=None,
    conditions=None,
    feature_symbols=None,
):
    if species not in {"human", "mouse"} or kind not in {"spatial", "reference"}:
        raise SpatialError("Declare human/mouse species and spatial/reference input kind.")
    _text(sample_id, "sample_id")
    if not isinstance(provenance, dict) or not provenance:
        raise SpatialError("Input source provenance is required.")
    ids, genes = list(observation_ids), list(features)
    if feature_symbols is not None:
        feature_symbols = list(feature_symbols)
        if len(feature_symbols) != len(genes):
            raise SpatialError("Declared gene symbols must match the feature ID axis.")
        for symbol in feature_symbols:
            if symbol is not None:
                _text(symbol, "gene symbol", limit=512)
    for axis in (ids, genes):
        if not axis or len(axis) != len(set(axis)):
            raise SpatialError("Count axes must be nonempty and unique; no silent aggregation.")
        for value in axis:
            _text(value, "axis ID", limit=512)
    x = sparse.csr_matrix(matrix, dtype=np.float64)
    x.sum_duplicates()
    x.eliminate_zeros()
    x.sort_indices()
    if x.shape != (len(ids), len(genes)) or max(x.shape) > 100000 or x.nnz > 25000000:
        raise SpatialError(
            "Input axes mismatch or sparse count budget exceeded (100k axes, 25M nonzero entries)."
        )
    if not np.isfinite(x.data).all() or np.any(x.data < 0) or np.any(x.data % 1):
        raise SpatialError("Algorithms require nonnegative integer raw counts, not normalized expression.")
    if np.any(np.asarray(x.sum(axis=1)).ravel() <= 0):
        raise SpatialError("Review zero-library observations before registering an analysis input.")
    xy = None
    if kind == "spatial":
        xy = np.asarray(coordinates, dtype=float)
        if xy.shape != (len(ids), 2) or not np.isfinite(xy).all():
            raise SpatialError("Spatial input requires exact finite N x 2 coordinates.")
        xy = xy.tolist()
    columns = {}
    for name, column in (("labels", labels), ("batches", batches), ("conditions", conditions)):
        if column is not None:
            column = list(column)
            if len(column) != len(ids):
                raise SpatialError(f"{name} does not match the observation axis.")
            for value in column:
                _text(value, name)
        columns[name] = column
    if kind == "reference" and (not columns["labels"] or len(set(columns["labels"])) < 2):
        raise SpatialError("Single-cell references require at least two explicitly annotated cell types.")
    buf = io.BytesIO()
    sparse.save_npz(buf, x, compressed=True)
    content = buf.getvalue()
    if len(content) > 96 * 1024**2:
        raise SpatialError("Compressed input exceeds 96 MiB; choose an explicit scientific subset.")
    return objects.put(
        project,
        "analysisinput",
        {
            "schema": "spatial-collab.analysis-input.v1",
            "kind": kind,
            "sample_id": sample_id,
            "species": species,
            "observation_ids": ids,
            "features": genes,
            "feature_symbols": feature_symbols,
            "coordinates": xy,
            **columns,
            "measurement": "raw_counts",
            "provenance": provenance,
            "shape": list(x.shape),
            "nnz": x.nnz,
            "counts_sha256": hashlib.sha256(content).hexdigest(),
            "counts_npz_base64": base64.b64encode(content).decode("ascii"),
        },
    )


def load_counts(project, input_id):
    record = objects.get(project, input_id, "analysisinput")
    try:
        content = base64.b64decode(record["counts_npz_base64"], validate=True)
        if hashlib.sha256(content).hexdigest() != record["counts_sha256"]:
            raise ValueError("checksum")
        with zipfile.ZipFile(io.BytesIO(content)) as archive:
            if sum(m.file_size for m in archive.infolist()) > 320 * 1024**2:
                raise ValueError("uncompressed size")
        matrix = sparse.load_npz(io.BytesIO(content)).tocsr()
        if list(matrix.shape) != record["shape"] or matrix.nnz != record["nnz"]:
            raise ValueError("shape")
    except (ValueError, KeyError, OSError, zipfile.BadZipFile) as exc:
        raise SpatialError("Analysis count snapshot integrity failed.") from exc
    return record, matrix


def snapshot(project, revision_id, species, observation_ids=None):
    cells = [c for c in project.cells(revision_id) if c["included"]]
    if observation_ids is not None:
        if (
            not observation_ids
            or len(set(observation_ids)) != len(observation_ids)
            or not set(observation_ids) <= {c["cell_id"] for c in cells}
        ):
            raise SpatialError("Explicit ROI must contain unique included IDs from this revision.")
        selected = set(observation_ids)
        cells = [c for c in cells if c["cell_id"] in selected]
    meta = project.summary()["metadata"]
    genes = meta["panel_genes"]
    symbols = {f["feature_id"]: f.get("symbol") for f in meta.get("features", [])}
    lookup = {g: i for i, g in enumerate(genes)}
    rows, cols, values = [], [], []
    for i, cell in enumerate(cells):
        for gene, value in cell["counts"].items():
            if value:
                rows.append(i)
                cols.append(lookup[gene])
                values.append(value)
    x = sparse.csr_matrix((values, (rows, cols)), shape=(len(cells), len(genes)))
    return register_counts(
        project,
        x,
        [c["cell_id"] for c in cells],
        genes,
        feature_symbols=[symbols.get(g) for g in genes] if any(symbols.get(g) for g in genes) else None,
        sample_id=meta.get("sample_id", meta["slice_id"]),
        species=species,
        kind="spatial",
        provenance={
            "source_sha256": project.context()["source_sha256"],
            "revision_id": revision_id,
            "revision_sha256": project.get_revision(revision_id)["revision_sha256"],
            "project_id": project.context()["project_id"],
            "units": meta["units"],
            "coordinate_system": meta["coordinate_system"],
            "scope": "explicit_roi" if observation_ids is not None else "whole_slice",
        },
        coordinates=[[c["x"], c["y"]] for c in cells],
        labels=[c["label"] for c in cells],
    )


def import_h5ad(
    project,
    path,
    *,
    sample_id,
    species,
    kind,
    counts_layer="counts",
    label_key=None,
    batch_key=None,
    condition_key=None,
    spatial_key="spatial",
    units="unknown",
    feature_id_key=None,
    feature_symbol_key=None,
):
    import anndata as ad

    path = Path(path).resolve(strict=True)
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    data = ad.read_h5ad(path)
    if counts_layer != "X" and counts_layer not in data.layers:
        raise SpatialError("Declared raw count layer is missing; X is not silently substituted.")

    def column(key):
        if key is None:
            return None
        if key not in data.obs or data.obs[key].isna().any():
            raise SpatialError(f"Missing declared observation metadata: {key}.")
        return data.obs[key].astype(str).tolist()

    def feature_column(key):
        if key == "_index":
            return data.var_names.tolist()
        if key not in data.var or data.var[key].isna().any():
            raise SpatialError(f"Missing declared feature metadata: {key}.")
        return data.var[key].astype(str).tolist()

    return register_counts(
        project,
        data.X if counts_layer == "X" else data.layers[counts_layer],
        data.obs_names.tolist(),
        feature_column(feature_id_key or "_index"),
        feature_symbols=feature_column(feature_symbol_key) if feature_symbol_key else None,
        sample_id=sample_id,
        species=species,
        kind=kind,
        provenance={
            "file_sha256": digest,
            "filename": path.name,
            "counts_layer": counts_layer,
            "units": units,
            "spatial_key": spatial_key,
            "feature_id_key": feature_id_key or "_index",
            "feature_symbol_key": feature_symbol_key,
        },
        coordinates=data.obsm.get(spatial_key),
        labels=column(label_key),
        batches=column(batch_key),
        conditions=column(condition_key),
    )


def list_inputs(project):
    return [
        {
            k: v
            for k, v in objects.get(project, oid, "analysisinput").items()
            if k
            not in {
                "counts_npz_base64",
                "observation_ids",
                "coordinates",
                "labels",
                "batches",
                "conditions",
                "features",
                "feature_symbols",
            }
        }
        for oid in objects.catalog(project, "analysisinput")
    ]
