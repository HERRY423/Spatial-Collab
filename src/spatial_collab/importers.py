"""Explicit, read-only adapters into a bounded analysis snapshot.

The snapshot is a working index. AnnData/Xenium source files remain unchanged;
it is not a replacement file format for image, segmentation or transcript data.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse

MAX_CELLS = 100_000
MAX_GENES = 20_000
MAX_NNZ = 2_000_000
MAX_H5AD_FEATURES = 100_000
MAX_H5AD_NNZ = 10_000_000
MAX_H5AD_SNAPSHOT_BYTES = 256 * 1024**2
MAX_FILE_BYTES = 2 * 1024**3
MAX_XENIUM_SOURCE_CELLS = 2_000_000
MAX_XENIUM_SOURCE_ID_CHARS = 64 * 1024**2
XENIUM_COUNT_FIELDS = ("transcript_counts", "control_probe_counts", "genomic_control_counts",
                       "control_codeword_counts", "unassigned_codeword_counts",
                       "deprecated_codeword_counts", "total_counts", "nucleus_count")
XENIUM_AREA_FIELDS = ("cell_area", "nucleus_area")


def source_digest(path: str | Path) -> dict:
    path = Path(path).resolve(strict=True)
    if not path.is_file() or path.stat().st_size > MAX_FILE_BYTES:
        raise ValueError("Source must be a file no larger than 2 GiB for this bounded alpha adapter.")
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return {"path": str(path), "sha256": h.hexdigest()}


def _stable_sources(before: list[dict]) -> None:
    for item in before:
        if source_digest(item["path"])["sha256"] != item["sha256"]:
            raise ValueError("Source changed during import. No project was created; retry from a frozen source.")


def _text(value, field: str) -> str:
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if not isinstance(value, str) or not value.strip() or len(value) > 500:
        raise ValueError(f"{field} must be a nonempty string of at most 500 characters.")
    return value


def _unique(values, field: str) -> list[str]:
    values = [_text(v, field) for v in values]
    if len(values) != len(set(values)):
        raise ValueError(f"Duplicate {field}; resolve ambiguous IDs before import.")
    return values


def _number(value, field: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"Non-numeric {field}; booleans are not coordinates.")
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError):
        raise ValueError(f"{field} must be a finite number.") from None
    if not math.isfinite(value):
        raise ValueError(f"Nonfinite {field}.")
    return value


def _metadata(name, slice_id, coordinate_system, panel, sources, source_kind, biological_replicates,
              *, units="micrometer"):
    if type(biological_replicates) is not int or biological_replicates not in (0, 1):
        raise ValueError("This adapter imports one slice: biological_replicates must be 0 (unknown) or 1.")
    return {
        "name": _text(name, "name"), "slice_id": _text(slice_id, "slice_id"),
        "coordinate_system": _text(coordinate_system, "coordinate_system"), "units": units,
        "panel_genes": panel, "source_kind": source_kind, "source_files": sources,
        "biological_replicates": biological_replicates,
        "limitations": [
            "Single-slice descriptive exploration; cells and ROIs are not independent biological replicates.",
            "Centroid review only: no morphology, segmentation-mask or individual-transcript validation.",
            "Input annotations are working hypotheses; editing does not establish a true cell identity.",
            "No patient-level, causal, differential-expression or cell-communication inference.",
        ],
    }


def _matrix_counts(matrix, ids, genes, *, max_nnz=MAX_NNZ):
    if matrix.shape != (len(ids), len(genes)):
        raise ValueError("Matrix dimensions do not match cell and feature IDs.")
    if sparse.issparse(matrix):
        if matrix.nnz > max_nnz:
            raise ValueError(f"Sparse matrix exceeds the declared {max_nnz} nonzero alpha budget; import a defined subset.")
        if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
            raise ValueError("Counts must be finite and nonnegative, including duplicate sparse entries.")
        if (matrix.data > 2**53 - 1).any():
            raise ValueError("Counts exceed the exact integer range supported by this adapter.")
        matrix = matrix.astype(np.float64).tocsr(copy=True)
        matrix.sum_duplicates()
        matrix.eliminate_zeros()
    else:
        if matrix.size > max_nnz:
            raise ValueError(f"Dense matrix exceeds {max_nnz} elements; provide sparse raw counts or a subset.")
        matrix = sparse.csr_matrix(matrix)
    if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
        raise ValueError("Counts must be finite and nonnegative.")
    if (matrix.data > 2**53 - 1).any():
        raise ValueError("Aggregated counts exceed the exact integer range supported by this adapter.")
    if not np.equal(matrix.data, np.floor(matrix.data)).all():
        raise ValueError("Raw integer counts required; choose a raw-counts layer rather than normalized values.")
    output = []
    for i in range(len(ids)):
        lo, hi = matrix.indptr[i:i+2]
        output.append({genes[int(j)]: float(v) for j, v in zip(matrix.indices[lo:hi], matrix.data[lo:hi])})
    return output


def _h5ad_snapshot_budget(cells, metadata):
    """Bound the eventual canonical JSON without making one giant JSON copy.

    The importer already supplies float coordinates/counts and every source-cell
    field, so store canonicalization changes ordering but not the encoded size.
    Reserve the newline written by review-bundle export as well.
    """
    def size(value):
        return len(json.dumps(value, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode("utf-8"))
    total = size({"metadata": metadata, "cells": []}) + max(0, len(cells) - 1) + 1
    for cell in cells:
        total += size(cell)
        if total > MAX_H5AD_SNAPSHOT_BYTES:
            raise ValueError("H5AD snapshot exceeds the 256 MiB replayable-source budget; declare a justified input scope rather than silently dropping features.")
    return total


def import_h5ad(path, destination, *, label_key="cell_type", spatial_key="spatial", counts_layer="counts",
                 slice_id, coordinate_system, units, name="Spatial transcriptomics", biological_replicates=0,
                 slice_key=None, allow_unannotated=False, platform="unknown", observation_unit="cell",
                 sample_id=None, feature_id_key=None, feature_symbol_key=None,
                 max_features=MAX_GENES, max_nnz=MAX_NNZ):
    """Import raw counts and one explicitly declared centroid coordinate frame.

    Use counts_layer='X' only if X really holds raw counts. No coordinate transform
    or unit inference occurs. A var subset is recorded as measured scope, not full panel.
    """
    from .store import Project
    from .identity import qualify_observation_id
    if not isinstance(units, str) or units not in {"micrometer", "pixel", "array_index", "unknown"}:
        raise ValueError("units must be micrometer, pixel, array_index or unknown; coordinates are never inferred or converted.")
    for value, hard_limit, field in ((max_features, MAX_H5AD_FEATURES, "max_features"),
                                     (max_nnz, MAX_H5AD_NNZ, "max_nnz")):
        if type(value) is not int or not 0 < value <= hard_limit:
            raise ValueError(f"{field} must be an integer in 1..{hard_limit}.")
    if sample_id is not None:
        sample_id = _text(sample_id, "sample_id")
    for field, value in (("feature_id_key", feature_id_key), ("feature_symbol_key", feature_symbol_key)):
        if value is not None:
            _text(value, field)
    if type(allow_unannotated) is not bool:
        raise ValueError("allow_unannotated must be a boolean.")
    if not isinstance(observation_unit, str) or observation_unit not in {"cell", "spot", "bin"}:
        raise ValueError("observation_unit must be cell, spot or bin.")
    import h5py
    from anndata.io import read_elem
    path = Path(path).resolve(strict=True)
    sources = [source_digest(path)]
    with h5py.File(path, "r") as handle:
        matrix_key = "X" if counts_layer == "X" else f"layers/{counts_layer}"
        if matrix_key not in handle:
            raise ValueError(f"Raw counts layer {counts_layer!r} missing; no implicit fallback to X.")
        matrix_node = handle[matrix_key]
        if isinstance(matrix_node, h5py.Group):
            shape = tuple(matrix_node.attrs.get("shape", ()))
            if "data" not in matrix_node or len(matrix_node["data"]) > max_nnz:
                raise ValueError(f"Sparse matrix exceeds the declared {max_nnz} nonzero import budget.")
        else:
            shape = matrix_node.shape
            if math.prod(shape) > max_nnz:
                raise ValueError(f"Dense matrix exceeds the declared {max_nnz} element import budget.")
        if len(shape) != 2 or not 0 < shape[0] <= MAX_CELLS or not 0 < shape[1] <= max_features:
            raise ValueError(f"Expected 1..100000 observations and 1..{max_features} features; use a declared bounded budget or subset.")
        obs, var = read_elem(handle["obs"]), read_elem(handle["var"])
        if len(obs) != shape[0] or len(var) != shape[1]:
            raise ValueError("Matrix and annotation dimensions differ.")
        if (label_key not in obs and not allow_unannotated) or f"obsm/{spatial_key}" not in handle:
            raise ValueError(f"Required obs[{label_key!r}] or obsm[{spatial_key!r}] missing.")
        # Common sample columns are checked even when the caller omitted slice_key.
        sample_columns = [k for k in ("slice_id", "sample_id", "library_id") if k in obs]
        if slice_key is not None and slice_key not in sample_columns:
            sample_columns.append(slice_key)
        for key in sample_columns:
            if key not in obs or obs[key].isna().any() or obs[key].nunique() != 1:
                raise ValueError("Mixed or missing slice identifiers; import one slice at a time.")
        if sample_id is not None and "sample_id" in obs and obs["sample_id"].iloc[0] != sample_id:
            raise ValueError("Declared sample_id differs from the source obs['sample_id']; preserve source identity.")
        source_ids = _unique(obs.index.tolist(), "cell ID")
        ids = [qualify_observation_id(sample_id, value) for value in source_ids] if sample_id is not None else source_ids
        for key in (feature_id_key, feature_symbol_key):
            if key is not None and key not in var:
                raise ValueError(f"Required var[{key!r}] missing; feature mapping is never guessed.")
        genes = _unique(var[feature_id_key].tolist() if feature_id_key is not None else var.index.tolist(),
                        "stable feature ID" if feature_id_key is not None else "gene name")
        symbols = [_text(value, "feature symbol") for value in
                   (var[feature_symbol_key].tolist() if feature_symbol_key is not None else var.index.tolist())]
        if label_key not in obs:
            labels = ["Unannotated"] * len(ids)
        else:
            missing_labels = obs[label_key].isna().tolist()
            labels = ["Unannotated" if allow_unannotated and (missing or isinstance(v, str) and not v.strip())
                      else _text(v, "observation label") for v, missing in zip(obs[label_key].tolist(), missing_labels)]
        xy = np.asarray(read_elem(handle[f"obsm/{spatial_key}"]))
        if xy.shape != (len(ids), 2) or not np.isfinite(xy).all():
            raise ValueError("Spatial coordinates must be finite N by 2 centroids in the declared units.")
        matrix = read_elem(matrix_node)
        counts = _matrix_counts(matrix, ids, genes, max_nnz=max_nnz)
        cells = [{"cell_id": cell_id, "x": float(xy[i, 0]), "y": float(xy[i, 1]),
                  "label": labels[i], "included": True, "region": "", "counts": counts[i],
                  **({"sample_id": sample_id, "source_cell_id": source_ids[i]} if sample_id is not None else {})}
                 for i, cell_id in enumerate(ids)]
    _stable_sources(sources)
    metadata = _metadata(name, slice_id, coordinate_system, genes, sources, "anndata_h5ad", biological_replicates,
                         units=units)
    metadata["features"] = [{"feature_id": feature_id, "symbol": symbol} for feature_id, symbol in zip(genes, symbols)]
    if sample_id is not None:
        metadata.update(sample_id=sample_id, identity_scope="sample_qualified")
    metadata["import_parameters"] = {"label_key": label_key, "spatial_key": spatial_key,
                                     "counts_layer": counts_layer, "slice_key": slice_key,
                                     "allow_unannotated": allow_unannotated, "sample_id": sample_id,
                                     "feature_id_key": feature_id_key, "feature_symbol_key": feature_symbol_key,
                                     "max_features": max_features, "max_nnz": max_nnz}
    metadata.update(platform=_text(platform, "platform"), observation_unit=observation_unit,
                    label_semantics={"cell": "cell_type", "spot": "spot_annotation", "bin": "bin_annotation"}[observation_unit],
                    annotation_status="unannotated" if all(v == "Unannotated" for v in labels) else "partially_unannotated" if "Unannotated" in labels else "provided_working_labels",
                    unannotated_observation_count_at_import=labels.count("Unannotated"))
    if "Unannotated" in labels:
        metadata["limitations"].append("Unannotated is an exploration placeholder, not an inferred biological group.")
    if observation_unit != "cell":
        metadata["limitations"].append(f"Observations are {observation_unit}s; labels are regional annotations, not pure cell identities or deconvolution.")
    metadata["limitations"].append("Only genes retained in this AnnData object are treated as measured.")
    if units != "micrometer":
        metadata["limitations"].append("Coordinate scale is unverified for physical distances; expression and annotation exploration remain available, but micrometre-radius analysis is unavailable.")
    _h5ad_snapshot_budget(cells, metadata)
    return Project.create(destination, cells, metadata)


def _read_csv(path: Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Missing or duplicate CSV headers.")
        rows = []
        for row in reader:
            if len(rows) >= MAX_CELLS:
                raise ValueError("CSV exceeds 100000-cell alpha budget; import a defined subset.")
            if None in row or any(value is None for value in row.values()):
                raise ValueError("CSV row width does not match its header.")
            rows.append(row)
        return rows


def _xenium_attributes(row):
    """Preserve provided instrument QC scalars, without creating confidence scores."""
    attributes = {}
    for field in XENIUM_COUNT_FIELDS + XENIUM_AREA_FIELDS:
        if field not in row:
            continue
        value = row[field]
        if value is None or str(value).strip().casefold() in {"", "nan", "na", "null", "none"}:
            attributes[field] = None
            continue
        number = _number(value, field)
        if number < 0 or (field in XENIUM_COUNT_FIELDS and (number != math.floor(number) or number > 2**53 - 1)):
            raise ValueError(f"Xenium {field} must be nonnegative" + (" exact integer count." if field in XENIUM_COUNT_FIELDS else "."))
        attributes[field] = int(number) if field in XENIUM_COUNT_FIELDS else number
    if "segmentation_method" in row:
        value = row["segmentation_method"]
        attributes["segmentation_method"] = _text(value, "segmentation_method") if value is not None and str(value).strip() else None
    return attributes


def _validated_bounds(bounds):
    if not isinstance(bounds, list) or len(bounds) != 4:
        raise ValueError("bounds must be [xmin, ymin, xmax, ymax] in native micrometers.")
    if any(isinstance(v, bool) or not isinstance(v, (int, float)) for v in bounds):
        raise ValueError("bounds must contain finite numeric coordinates, not booleans or strings.")
    values = [_number(value, "bounds") for value in bounds]
    if values[0] >= values[2] or values[1] >= values[3]:
        raise ValueError("bounds must have strictly positive width and height.")
    return values


def _import_xenium_h5_window(outs, annotations, destination, *, bounds, slice_id, name,
                            biological_replicates, allow_unannotated):
    """Stream full centroid IDs, then read only selected HDF5 CSC columns.

    The explicit window defines the complete imported observation universe. Its
    exterior is unavailable to downstream whole-slice metrics in this project.
    """
    import h5py
    from .format_readers import _annotations, _csv_rows, _one_file
    from .store import Project
    limits = _validated_bounds(bounds)
    csv_path = _one_file(outs, ["cells.csv", "cells.csv.gz"], "Xenium cell summary")
    matrix_path = outs / "cell_feature_matrix.h5"
    if not matrix_path.is_file():
        raise ValueError("Bounded Xenium spatial-window import currently requires cell_feature_matrix.h5.")
    sources = [source_digest(path) for path in (csv_path, matrix_path)]
    selected, remaining_ids, source_count, source_id_chars = {}, set(), 0, 0
    for row in _csv_rows(csv_path, max_rows=MAX_XENIUM_SOURCE_CELLS):
        if not {"cell_id", "x_centroid", "y_centroid"} <= set(row):
            raise ValueError("Xenium cell summary requires cell_id,x_centroid,y_centroid.")
        identifier = _text(row["cell_id"], "cell_id")
        if identifier in remaining_ids:
            raise ValueError(f"Duplicate source cell ID: {identifier}.")
        source_id_chars += len(identifier)
        if source_id_chars > MAX_XENIUM_SOURCE_ID_CHARS:
            raise ValueError("Source identity index exceeds the bounded memory budget.")
        remaining_ids.add(identifier)
        source_count += 1
        x, y = _number(row["x_centroid"], "x_centroid"), _number(row["y_centroid"], "y_centroid")
        if limits[0] <= x <= limits[2] and limits[1] <= y <= limits[3]:
            if len(selected) >= MAX_CELLS:
                raise ValueError("Spatial window exceeds 100000 selected cells; choose a smaller declared window.")
            attributes = _xenium_attributes(row)
            selected[identifier] = {"cell_id": identifier, "x": x, "y": y,
                                    **({"attributes": attributes} if attributes else {})}
    if not selected:
        raise ValueError("The declared spatial window contains no cell centroids.")

    with h5py.File(matrix_path, "r") as handle:
        if "matrix" not in handle:
            raise ValueError("Xenium H5 requires a matrix group.")
        group = handle["matrix"]
        shape = tuple(int(v) for v in group["shape"][:])
        if len(shape) != 2 or not 0 < shape[0] <= MAX_GENES or shape[1] != source_count:
            raise ValueError("Xenium full matrix dimensions must match the source centroid count and feature budget.")
        if any(key not in group for key in ("barcodes", "features", "data", "indices", "indptr")):
            raise ValueError("Xenium H5 is missing required CSC or feature metadata.")
        if len(group["barcodes"]) != shape[1] or len(group["indptr"]) != shape[1] + 1:
            raise ValueError("Xenium barcode/column metadata dimensions do not match the full matrix.")
        names = [_text(value, "feature name") for value in group["features"]["name"][:]]
        types = [_text(value, "feature type") for value in group["features"]["feature_type"][:]]
        if len(names) != shape[0] or len(types) != shape[0]:
            raise ValueError("Xenium feature metadata dimensions do not match the full matrix.")
        keep = [i for i, value in enumerate(types) if value == "Gene Expression"]
        genes = _unique([names[i] for i in keep], "gene name")
        if not genes:
            raise ValueError("No Gene Expression features found.")
        ids, columns = [], []
        # Chunked barcode scan validates the entire source identity join without
        # loading the source count arrays or keeping a second full string index.
        for start in range(0, shape[1], 8192):
            for offset, value in enumerate(group["barcodes"][start:start + 8192]):
                identifier = _text(value, "matrix cell ID")
                if identifier not in remaining_ids:
                    raise ValueError("Matrix has duplicate or unmatched source barcodes.")
                remaining_ids.remove(identifier)
                if identifier in selected:
                    ids.append(identifier)
                    columns.append(start + offset)
        if remaining_ids:
            raise ValueError("Source centroid IDs and matrix barcodes do not match exactly.")
        pointers = group["indptr"][:]
        if (pointers.dtype.kind not in "iu" or pointers.ndim != 1 or pointers[0] != 0
                or (pointers[1:] < pointers[:-1]).any()
                or int(pointers[-1]) != len(group["data"]) or len(group["data"]) != len(group["indices"])):
            raise ValueError("Malformed full CSC column pointers or stored-entry dimensions.")
        lengths = [int(pointers[column + 1]) - int(pointers[column]) for column in columns]
        selected_nnz = sum(lengths)
        if selected_nnz > MAX_NNZ:
            raise ValueError("Spatial window exceeds 2 million stored matrix entries; choose a smaller declared window.")
        data = np.empty(selected_nnz, dtype=np.float64)
        indices = np.empty(selected_nnz, dtype=np.int32)
        subset_pointers = np.zeros(len(columns) + 1, dtype=np.int64)
        cursor = 0
        for i, (column, length) in enumerate(zip(columns, lengths)):
            lo, hi = int(pointers[column]), int(pointers[column + 1])
            block_indices = group["indices"][lo:hi]
            if block_indices.dtype.kind not in "iu" or (block_indices < 0).any() or (block_indices >= shape[0]).any():
                raise ValueError("Selected CSC feature indices are outside the declared matrix.")
            block_data = group["data"][lo:hi]
            if (not np.isfinite(block_data).all() or (block_data < 0).any()
                    or (block_data > 2**53 - 1).any() or not np.equal(block_data, np.floor(block_data)).all()):
                raise ValueError("Selected raw counts must be finite nonnegative exact integers.")
            data[cursor:cursor + length] = block_data
            indices[cursor:cursor + length] = block_indices
            cursor += length
            subset_pointers[i + 1] = cursor
        matrix = sparse.csc_matrix((data, indices, subset_pointers), shape=(shape[0], len(columns)))
        matrix.check_format(full_check=True)
        counts = _matrix_counts(matrix[keep, :].T.tocsr(), ids, genes)
        source_nnz = int(pointers[-1])
    labels = _annotations(annotations, ids, allow_unannotated, sources)
    cells = [{**selected[identifier], "label": labels[identifier], "counts": count}
             for identifier, count in zip(ids, counts)]
    _stable_sources(sources)
    metadata = _metadata(name, slice_id, "xenium_native_xy", genes, sources, "xenium_h5_csv", biological_replicates)
    unannotated = sum(label == "Unannotated" for label in labels.values())
    metadata.update(platform="Xenium", observation_unit="cell", label_semantics="cell_type",
                    annotation_status="unannotated" if unannotated == len(labels) else "partially_unannotated" if unannotated else "provided_working_labels",
                    unannotated_observation_count_at_import=unannotated,
                    import_scope={"kind": "spatial_window", "bounds": limits, "units": "micrometer",
                                  "source_observation_count": source_count, "selected_observation_count": len(cells),
                                  "selection_rule": "all cell centroids with xmin <= x <= xmax and ymin <= y <= ymax",
                                  "sampling": "none", "gene_selection": "all measured Gene Expression features",
                                  "source_stored_entries": source_nnz, "selected_stored_entries_including_controls": selected_nnz},
                    import_parameters={"feature_type": "Gene Expression", "matrix_format": "10x_h5_selected_csc_columns",
                                       "bounds": limits, "allow_unannotated": allow_unannotated},
                    preserved_attribute_fields=sorted({key for cell in cells for key in cell.get("attributes", {})}))
    metadata["limitations"].append("Imported observations are a declared spatial window, not the complete source slice; whole_slice graph scope means only this imported universe. Neighbors outside the window are unavailable.")
    metadata["limitations"].append("Instrument QC attributes are immutable reported values; no automatic quality threshold or confidence score was inferred.")
    if unannotated:
        metadata["limitations"].append("Unannotated is an exploration placeholder, not an inferred biological group.")
    return Project.create(destination, cells, metadata)


def import_xenium(outs, annotations, destination, *, slice_id, name="Xenium review", biological_replicates=0,
                  allow_unannotated=False, bounds=None):
    """Read official cells.csv[.gz] + 10x cell_feature_matrix.h5 and explicit labels.

    The annotation CSV must contain exactly one row per cell with cell_id,label.
    Controls and protein features are excluded using the HDF5 feature_type field.
    """
    from .store import Project
    if type(allow_unannotated) is not bool:
        raise ValueError("allow_unannotated must be a boolean.")
    outs = Path(outs).resolve(strict=True)
    if bounds is not None:
        return _import_xenium_h5_window(outs, annotations, destination, bounds=bounds, slice_id=slice_id,
                                       name=name, biological_replicates=biological_replicates,
                                       allow_unannotated=allow_unannotated)
    if not (outs / "cell_feature_matrix.h5").is_file():
        from .format_readers import import_xenium_mtx
        return import_xenium_mtx(outs, annotations, destination, slice_id=slice_id, name=name,
                                 biological_replicates=biological_replicates, allow_unannotated=allow_unannotated)
    import h5py
    csv_path = outs / "cells.csv.gz"
    if not csv_path.exists():
        csv_path = outs / "cells.csv"
    matrix_path = outs / "cell_feature_matrix.h5"
    if annotations is None and not allow_unannotated:
        raise ValueError("Annotations are required unless allow_unannotated=true is explicit.")
    annotations = Path(annotations).resolve(strict=True) if annotations is not None else None
    sources = [source_digest(p) for p in (csv_path, matrix_path) + ((annotations,) if annotations else ())]
    rows = _read_csv(csv_path)
    annotations_rows = _read_csv(annotations) if annotations else [{"cell_id": row["cell_id"], "label": "Unannotated"} for row in rows]
    if not rows or not annotations_rows:
        raise ValueError("Cells and annotations must be nonempty.")
    ids = _unique([row.get("cell_id") for row in rows], "cell ID")
    ann_ids = _unique([row.get("cell_id") for row in annotations_rows], "annotation cell ID")
    if set(ids) != set(ann_ids):
        raise ValueError("Annotation cell IDs must match Xenium cells exactly; no implicit dropping or relabeling.")
    labels = {row["cell_id"]: "Unannotated" if allow_unannotated and not str(row.get("label", "")).strip()
              else _text(row.get("label"), "label") for row in annotations_rows}
    with h5py.File(matrix_path, "r") as handle:
        group = handle["matrix"]
        shape = tuple(int(i) for i in group["shape"][:])
        if len(shape) != 2 or shape[0] > MAX_GENES or shape[1] > MAX_CELLS:
            raise ValueError("Xenium matrix exceeds bounded alpha size.")
        if len(group["data"]) > MAX_NNZ:
            raise ValueError("Xenium matrix exceeds 2 million nonzero alpha budget.")
        barcodes = _unique(group["barcodes"][:], "matrix cell ID")
        if set(ids) != set(barcodes):
            raise ValueError("Matrix barcodes must match cells exactly.")
        feature_names = [_text(v, "feature name") for v in group["features"]["name"][:]]
        feature_types = [_text(v, "feature type") for v in group["features"]["feature_type"][:]]
        if len(barcodes) != shape[1] or len(feature_names) != shape[0] or len(feature_types) != shape[0]:
            raise ValueError("Xenium feature/barcode metadata dimensions must match the entire matrix.")
        keep = [i for i, kind in enumerate(feature_types) if kind == "Gene Expression"]
        genes = _unique([feature_names[i] for i in keep], "gene name")
        if not genes:
            raise ValueError("No Gene Expression features found.")
        matrix = sparse.csc_matrix((group["data"][:], group["indices"][:], group["indptr"][:]), shape=shape)
        matrix.check_format(full_check=True)
        barcode_index = {v: i for i, v in enumerate(barcodes)}
        matrix = matrix[keep, :][:, [barcode_index[v] for v in ids]].T.tocsr()
        counts = _matrix_counts(matrix, ids, genes)
    cells = [{"cell_id": cell_id, "x": _number(row["x_centroid"], "x_centroid"),
              "y": _number(row["y_centroid"], "y_centroid"), "label": labels[cell_id],
              "included": True, "region": "", "counts": counts[i]}
             for i, (cell_id, row) in enumerate(zip(ids, rows))]
    for cell, row in zip(cells, rows):
        attributes = _xenium_attributes(row)
        if attributes:
            cell["attributes"] = attributes
    _stable_sources(sources)
    metadata = _metadata(name, slice_id, "xenium_native_xy", genes, sources, "xenium_h5_csv",
                         biological_replicates)
    metadata["import_parameters"] = {"feature_type": "Gene Expression", "annotation_columns": ["cell_id", "label"],
                                     "allow_unannotated": allow_unannotated, "matrix_format": "10x_h5"}
    unannotated = sum(label == "Unannotated" for label in labels.values())
    metadata.update(platform="Xenium", observation_unit="cell", label_semantics="cell_type",
                    annotation_status="unannotated" if unannotated == len(labels) else "partially_unannotated" if unannotated else "provided_working_labels",
                    unannotated_observation_count_at_import=unannotated)
    metadata["preserved_attribute_fields"] = sorted({key for cell in cells for key in cell.get("attributes", {})})
    if unannotated:
        metadata["limitations"].append("Unannotated is an exploration placeholder, not an inferred biological group.")
    return Project.create(destination, cells, metadata)
