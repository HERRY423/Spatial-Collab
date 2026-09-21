"""Bounded local, read-only image and segmentation registrations.

Pixels are indexed (row=y, column=x); a required affine maps pixel centers
(x, y) to the explicitly declared project frame. Registration records a
declaration of matching origin, not evidence of biological registration quality.
No array position is ever interpreted as a cell identifier.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any

import numpy as np

from .store import Project, SpatialError

MAX_PIXELS = 32_000_000
MAX_ARRAY_BYTES = 256 * 1024**2
MAX_FILE_BYTES = 2 * 1024**3
MAX_RECORD_BYTES = 32 * 1024**2
MAX_ASSETS = 32
MAX_LABELS = 250_000
KEY = "spatial_collab_asset"
SCHEMA = "spatial-collab.asset.v1"


def _project(project):
    return project if isinstance(project, Project) else Project(project)


def _json(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _hash(value):
    return hashlib.sha256(_json(value).encode()).hexdigest()


def _text(value, field, limit=1000):
    if not isinstance(value, str) or not value.strip() or len(value) > limit:
        raise SpatialError(f"{field} must be nonempty text of at most {limit} characters.")
    return value


def _file(path):
    try:
        path = Path(path).resolve(strict=True)
    except (OSError, TypeError) as exc:
        raise SpatialError("Asset source file is unavailable.") from exc
    if not path.is_file() or path.stat().st_size > MAX_FILE_BYTES:
        raise SpatialError("Asset source must be a local file no larger than 2 GiB.")
    return path


def file_digest(path):
    path = _file(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024**2), b""):
            digest.update(chunk)
    return {"path": str(path), "sha256": digest.hexdigest()}


def validate_array_budget(shape, dtype):
    dtype = np.dtype(dtype)
    if len(shape) != 2 or any(type(n) is not int or n <= 0 for n in shape):
        raise SpatialError("Asset must be an explicitly selected, nonempty two-dimensional plane.")
    if dtype.kind not in "uif" or dtype.hasobject:
        raise SpatialError("Asset pixels must be numeric; object, boolean and complex arrays are unsupported.")
    count = math.prod(shape)
    if count > MAX_PIXELS or count * dtype.itemsize > MAX_ARRAY_BYTES:
        raise SpatialError("Asset exceeds the bounded pixel/memory budget; create an explicitly located crop first.")


def array_digest(data):
    if not isinstance(data, np.ndarray):
        raise SpatialError("Only materialized NumPy arrays are accepted; lazy arrays are never implicitly loaded.")
    validate_array_budget(data.shape, data.dtype)
    digest = hashlib.sha256()
    digest.update(_json({"shape": list(data.shape), "dtype": data.dtype.str}).encode())
    # Iteration avoids a second full contiguous copy of noncontiguous arrays.
    for row in data:
        digest.update(np.ascontiguousarray(row).tobytes())
    return digest.hexdigest()


def validate_affine(value):
    try:
        matrix = np.asarray(value, dtype=float)
    except (ValueError, TypeError, OverflowError) as exc:
        raise SpatialError("pixel_to_world must be an explicit numeric 3x3 affine.") from exc
    if matrix.shape != (3, 3) or not np.isfinite(matrix).all() or not np.array_equal(matrix[2], [0, 0, 1]):
        raise SpatialError("pixel_to_world must be finite, invertible and two-dimensional affine.")
    with np.errstate(over="ignore", invalid="ignore"):
        determinant = np.linalg.det(matrix[:2, :2])
    if not np.isfinite(determinant) or abs(determinant) < 1e-12:
        raise SpatialError("pixel_to_world must be finite, invertible and two-dimensional affine.")
    return matrix


def _load(path, tiff_page):
    path = _file(path)
    if path.suffix.lower() == ".npy":
        if tiff_page is not None:
            raise SpatialError("tiff_page is only valid for TIFF sources.")
        try:
            source = np.load(path, mmap_mode="r", allow_pickle=False)
        except (ValueError, OSError) as exc:
            raise SpatialError("NPY source is not a bounded numeric 2D plane.") from exc
        validate_array_budget(source.shape, source.dtype)
        data = np.array(source, copy=True)
    elif path.suffix.lower() in {".tif", ".tiff"}:
        try:
            import tifffile
        except ImportError as exc:
            raise SpatialError("TIFF reading requires the optional tifffile package.") from exc
        if tiff_page is not None and (type(tiff_page) is not int or tiff_page < 0):
            raise SpatialError("tiff_page must be an explicit nonnegative integer.")
        with tifffile.TiffFile(path) as source:
            if tiff_page is None and len(source.pages) != 1:
                raise SpatialError("Multiple TIFF pages require an explicit tiff_page; no projection is inferred.")
            index = tiff_page if tiff_page is not None else 0
            if index >= len(source.pages):
                raise SpatialError("Requested TIFF page is unavailable.")
            page = source.pages[index]
            validate_array_budget(page.shape, page.dtype)
            data = page.asarray()
    else:
        raise SpatialError("Only local NPY and TIFF planes are supported by this bounded asset reader.")
    if not np.isfinite(data).all():
        raise SpatialError("Asset contains nonfinite pixels.")
    data.flags.writeable = False
    return data


def _binding(project, supplied):
    summary = project.summary()
    metadata = summary["metadata"]
    expected = {"source_sha256": summary["source_sha256"], "slice_id": metadata["slice_id"],
                "coordinate_system": metadata["coordinate_system"], "units": metadata["units"]}
    if metadata.get("sample_id") is not None:
        expected["sample_id"] = metadata["sample_id"]
    if supplied != expected:
        raise SpatialError("Asset source, slice, coordinate frame or units do not exactly match the project.")
    return expected


def _mapping(project, data, kind, supplied):
    if kind != "segmentation":
        if supplied is not None:
            raise SpatialError("Only segmentation assets accept label_to_cell_id.")
        return {}, 0, 0
    if data.dtype.kind not in "ui" or (data < 0).any() or (data > 2**53 - 1).any():
        raise SpatialError("Segmentation pixels must be nonnegative exact integers, with background 0.")
    present = {int(label) for label in np.unique(data) if label}
    if len(present) > MAX_LABELS:
        raise SpatialError("Segmentation label budget exceeded.")
    if supplied is None:
        return {}, len(present), len(present)
    if not isinstance(supplied, dict) or len(supplied) > MAX_LABELS:
        raise SpatialError("label_to_cell_id must be a bounded explicit mapping.")
    ids = {cell["cell_id"] for cell in project.cells()}
    mapping = {}
    for label, cell_id in supplied.items():
        if type(label) is int and label > 0:
            label = str(label)
        if not isinstance(label, str) or not re.fullmatch(r"[1-9][0-9]{0,15}", label):
            raise SpatialError("Mask labels must be canonical positive integers; background 0 cannot map to a cell.")
        if label in mapping or int(label) not in present:
            raise SpatialError("Duplicate or absent mask label in explicit cell mapping.")
        if not isinstance(cell_id, str) or cell_id not in ids:
            raise SpatialError("Mask mapping refers to an unknown exact project cell_id.")
        mapping[label] = cell_id
    return mapping, len(present), len(present) - len(mapping)


def register_asset(project, path, *, kind, source_sha256, slice_id, coordinate_system, units,
                   pixel_to_world, registration_note, label_to_cell_id=None, tiff_page=None,
                   name=None, origin_sources=None, sample_id=None):
    """Register a local source after explicit provenance/frame declaration.

    The only writes are immutable JSON sidecars in project/assets. TIFF/NPY,
    observations, shared selection and annotation revisions are unchanged.
    """
    project = _project(project)
    if kind not in {"image", "segmentation"}:
        raise SpatialError("Asset kind must be image or segmentation.")
    supplied = dict(source_sha256=source_sha256, slice_id=slice_id, coordinate_system=coordinate_system, units=units)
    if sample_id is not None:
        supplied["sample_id"] = sample_id
    binding = _binding(project, supplied)
    matrix = validate_affine(pixel_to_world)
    note = _text(registration_note, "registration_note", 4000)
    source = file_digest(path)
    data = _load(source["path"], tiff_page)
    mapping, label_count, missing = _mapping(project, data, kind, label_to_cell_id)
    if file_digest(source["path"]) != source:
        raise SpatialError("Asset source changed during registration; no record was written.")
    origins = []
    if origin_sources is not None:
        if not isinstance(origin_sources, list) or len(origin_sources) > 8:
            raise SpatialError("origin_sources must contain at most eight path/sha256 records.")
        for item in origin_sources:
            if not isinstance(item, dict) or set(item) != {"path", "sha256"}:
                raise SpatialError("Origin sources require exact path and sha256 fields.")
            actual = file_digest(item["path"])
            if actual["sha256"] != item["sha256"]:
                raise SpatialError("Origin source digest differs from its explicit declaration.")
            origins.append(actual)
    record = {"asset_schema": SCHEMA, **binding, "kind": kind,
              "name": _text(name or Path(path).name, "name", 500), "source": source,
              "origin_sources": origins, "registration_note": note,
              "registration_status": "DECLARED_NOT_INDEPENDENTLY_VALIDATED",
              "pixel_to_world": matrix.tolist(), "pixel_axis_order": ["y", "x"],
              "transform_axis_order": ["x", "y"], "pixel_reference": "integer coordinates are pixel centers",
              "tiff_page": tiff_page, "shape": list(data.shape), "dtype": data.dtype.str,
              "array_sha256": array_digest(data), "label_to_cell_id": mapping,
              "label_count": label_count, "unmapped_label_count": missing,
              "scientific_authorization": "NOT_ESTABLISHED",
              "limitations": ["Source registration is a declared binding, not independent tissue identification.",
                              "Missing image, mask, mapping or calibrated frame leaves the corresponding check unknown.",
                              "Pixel/cell correspondence does not establish a merged cell, doublet or true coexpression.",
                              "This is a 2D plane; no 3D segmentation validity is inferred."]}
    record["asset_id"] = "asset_" + _hash(record)[:24]
    record["record_sha256"] = _hash(record)
    directory = project.root / "assets"
    directory.mkdir(exist_ok=True)
    destination = directory / (record["asset_id"] + ".json")
    encoded = _json(record)
    if len(encoded.encode()) > MAX_RECORD_BYTES:
        raise SpatialError("Asset registration record exceeds the bounded record budget.")
    if destination.exists():
        if destination.read_text(encoding="utf-8") != encoded:
            raise SpatialError("An existing immutable asset record has inconsistent contents.")
        return _public(record)
    if len(list(directory.glob("asset_*.json"))) >= MAX_ASSETS:
        raise SpatialError("Project asset budget exceeded.")
    with destination.open("x", encoding="utf-8") as stream:
        stream.write(encoded)
    return _public(record)


def _read(project, asset_id):
    if not isinstance(asset_id, str) or not re.fullmatch(r"asset_[a-f0-9]{24}", asset_id):
        raise SpatialError("Asset ID is invalid; file paths are not accepted here.")
    path = project.root / "assets" / (asset_id + ".json")
    if not path.is_file() or path.stat().st_size > MAX_RECORD_BYTES:
        raise SpatialError("Asset registration is unavailable or exceeds the record budget.")
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
        unhashed = {key: value for key, value in record.items() if key != "record_sha256"}
        original = {key: value for key, value in unhashed.items() if key != "asset_id"}
        if (record["record_sha256"] != _hash(unhashed) or record["asset_id"] != asset_id
                or record["asset_schema"] != SCHEMA or "asset_" + _hash(original)[:24] != asset_id):
            raise ValueError("inconsistent record")
        supplied = {key: record[key] for key in ("source_sha256", "slice_id", "coordinate_system", "units")}
        if "sample_id" in record:
            supplied["sample_id"] = record["sample_id"]
        _binding(project, supplied)
        validate_affine(record["pixel_to_world"])
    except (ValueError, KeyError, TypeError) as exc:
        raise SpatialError("Asset registration integrity check failed.") from exc
    return record


def _public(record):
    result = deepcopy(record)
    mapping = result.pop("label_to_cell_id")
    result["mapped_label_count"] = len(mapping)
    result["mapping_status"] = ("not_applicable" if record["kind"] == "image" else
                                "missing" if not mapping else
                                "partial" if record["unmapped_label_count"] else "complete")
    return result


def list_assets(project):
    project = _project(project)
    directory = project.root / "assets"
    paths = sorted(directory.glob("asset_*.json")) if directory.exists() else []
    if len(paths) > MAX_ASSETS:
        raise SpatialError("Project asset budget exceeded.")
    assets = [_public(_read(project, path.stem)) for path in paths]
    return {"assets_schema": "spatial-collab.assets.v1", "assets": assets,
            "content_integrity": "not_checked_until_inspection_or_loading",
            "image_available": any(a["kind"] == "image" for a in assets),
            "segmentation_available": any(a["kind"] == "segmentation" for a in assets),
            "missing_evidence_status": "UNKNOWN", "scientific_authorization": "NOT_ESTABLISHED"}


def _verified(project, asset_id):
    record = _read(project, asset_id)
    for source in [record["source"], *record["origin_sources"]]:
        if file_digest(source["path"]) != source:
            raise SpatialError("Registered asset or origin file changed; refresh cannot repair its declared identity.")
    data = _load(record["source"]["path"], record["tiff_page"])
    if array_digest(data) != record["array_sha256"] or file_digest(record["source"]["path"]) != record["source"]:
        raise SpatialError("Asset pixel content differs from the registered snapshot.")
    return record, data


def inspect_asset(project, asset_id, cell_ids=None):
    project = _project(project)
    record, data = _verified(project, asset_id)
    if cell_ids is None:
        cell_ids = []
    if (not isinstance(cell_ids, list) or len(cell_ids) > 100 or any(not isinstance(c, str) for c in cell_ids)
            or len(cell_ids) != len(set(cell_ids))):
        raise SpatialError("Inspect at most 100 distinct exact cell IDs per request.")
    summary = project.summary()
    cells = {c["cell_id"]: c for c in project.cells(summary["head_revision"])}
    if set(cell_ids) - cells.keys():
        raise SpatialError("Unknown exact project cell_id in asset inspection.")
    inverse = np.linalg.inv(record["pixel_to_world"])
    mapping = record["label_to_cell_id"]
    observations = []
    for cid in cell_ids:
        cell = cells[cid]
        position = inverse @ [cell["x"], cell["y"], 1]
        if not np.isfinite(position).all():
            raise SpatialError("Transformed cell coordinates exceed finite numeric range.")
        x, y = [math.floor(float(v) + .5) for v in position[:2]]
        inside = 0 <= y < data.shape[0] and 0 <= x < data.shape[1]
        item = {"cell_id": cid, "world_xy": [cell["x"], cell["y"]],
                "pixel_center_xy": position[:2].tolist(), "nearest_pixel_xy": [x, y],
                "inside_plane": inside, "correspondence": "outside_image" if not inside else "image_only"}
        if inside:
            item["pixel_value"] = data[y, x].item()
            if record["kind"] == "segmentation":
                label = int(data[y, x])
                mapped = mapping.get(str(label))
                item.update(label_id=label, mapped_cell_id=mapped,
                            correspondence="background_at_centroid" if not label else
                            "mapping_unknown" if mapped is None else
                            "same_declared_cell" if mapped == cid else "different_declared_cell")
        observations.append(item)
    return {"inspection_schema": "spatial-collab.asset-inspection.v1", "asset": _public(record),
            "revision_id": summary["head_revision"], "content_integrity": "verified",
            "observations": observations, "changes_applied": False,
            "merge_or_true_coexpression": "NOT_DETERMINED", "scientific_authorization": "NOT_ESTABLISHED"}


def load_asset_layer(project, asset_id):
    project = _project(project)
    record, data = _verified(project, asset_id)
    swap = np.array([[0., 1, 0], [1, 0, 0], [0, 0, 1]])
    affine_yx = swap @ np.asarray(record["pixel_to_world"]) @ swap
    kwargs = {"name": record["name"], "metadata": {KEY: _public(record)},
              "affine": affine_yx, "opacity": .5 if record["kind"] == "segmentation" else 1.0}
    return data, kwargs, "labels" if record["kind"] == "segmentation" else "image"


def validate_layer_capture(project, asset_id, capture):
    record = _read(_project(project), asset_id)
    if capture["metadata"].get(KEY) != _public(record):
        raise SpatialError("Visible asset metadata changed; reload the registered layer.")
    swap = np.array([[0., 1, 0], [1, 0, 0], [0, 0, 1]])
    expected = swap @ np.asarray(record["pixel_to_world"]) @ swap
    if not np.allclose(capture["transform"], expected, rtol=1e-12, atol=1e-12):
        raise SpatialError("Visible image or segmentation was transformed; reload before comparing cell locations.")
    if array_digest(capture["data"]) != record["array_sha256"]:
        raise SpatialError("Visible image or segmentation pixels changed; reload the registered layer.")


def register_array_snapshot(project, data, *, origin_sources, **declaration):
    """Freeze an existing materialized layer; no implicit loading of lazy data."""
    project = _project(project)
    supplied = {key: declaration[key] for key in ("source_sha256", "slice_id", "coordinate_system", "units")}
    if declaration.get("sample_id") is not None:
        supplied["sample_id"] = declaration["sample_id"]
    _binding(project, supplied)
    validate_affine(declaration["pixel_to_world"])
    digest = array_digest(data)
    directory = project.root / "assets" / "snapshots"
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / (digest + ".npy")
    if not destination.exists():
        # Never overwrite a registered snapshot, even after a viewer edit.
        with destination.open("xb") as stream:
            np.save(stream, data, allow_pickle=False)
    elif array_digest(_load(destination, None)) != digest:
        raise SpatialError("Existing frozen layer snapshot differs from its content identity.")
    return register_asset(project, destination, origin_sources=origin_sources, **declaration)


def _encode_png_rgba(rgba: np.ndarray) -> bytes:
    """Encode an H x W x 4 uint8 numpy array into a valid PNG without third-party dependencies."""
    import struct
    import zlib
    h, w, _ = rgba.shape
    raw_scanlines = bytearray()
    for row in rgba:
        raw_scanlines.append(0)  # Filter type 0 (None)
        raw_scanlines.extend(row.tobytes())
    compressed = zlib.compress(bytes(raw_scanlines), level=6)

    def chunk(tag: bytes, data: bytes) -> bytes:
        length = struct.pack(">I", len(data))
        crc = struct.pack(">I", zlib.crc32(tag + data) & 0xffffffff)
        return length + tag + data + crc

    sig = b"\x89PNG\r\n\x1a\n"
    ihdr_data = struct.pack(">IIBBBBB", w, h, 8, 6, 0, 0, 0)
    ihdr = chunk(b"IHDR", ihdr_data)
    idat = chunk(b"IDAT", compressed)
    iend = chunk(b"IEND", b"")
    return sig + ihdr + idat + iend


def _hsv_to_rgb(h: float, s: float, v: float) -> tuple[float, float, float]:
    i = int(h * 6.0)
    f = (h * 6.0) - i
    p = v * (1.0 - s)
    q = v * (1.0 - s * f)
    t = v * (1.0 - s * (1.0 - f))
    i = i % 6
    if i == 0:
        return v, t, p
    if i == 1:
        return q, v, p
    if i == 2:
        return p, v, t
    if i == 3:
        return p, q, v
    if i == 4:
        return t, p, v
    return v, p, q


def get_asset_raster(project, asset_id: str, max_dimension: int = 1024) -> dict[str, Any]:
    """Render a downsampled, verified 2D asset plane to a lightweight base64 PNG data URL."""
    import base64
    project = _project(project)
    record, data = _verified(project, asset_id)
    h, w = data.shape
    step = max(1, math.ceil(max(h, w) / max_dimension))
    downsampled = data[::step, ::step]
    dh, dw = downsampled.shape

    rgba = np.zeros((dh, dw, 4), dtype=np.uint8)
    if record["kind"] == "image":
        if np.issubdtype(downsampled.dtype, np.floating) or np.issubdtype(downsampled.dtype, np.integer):
            p1, p99 = np.percentile(downsampled, (1, 99))
            if p99 > p1:
                norm = np.clip((downsampled.astype(float) - p1) / (p99 - p1), 0.0, 1.0)
            else:
                norm = np.zeros_like(downsampled, dtype=float)
            gray = (norm * 255.0).astype(np.uint8)
            rgba[..., 0] = gray
            rgba[..., 1] = gray
            rgba[..., 2] = gray
            rgba[..., 3] = 255
    else:  # segmentation
        labels = downsampled.astype(np.int64)
        unique_labels = np.unique(labels[labels > 0])
        for lbl in unique_labels:
            hue = (int(lbl) * 0.618033988749895) % 1.0
            r, g, b = _hsv_to_rgb(hue, 0.75, 0.9)
            mask = labels == lbl
            rgba[mask, 0] = int(r * 255)
            rgba[mask, 1] = int(g * 255)
            rgba[mask, 2] = int(b * 255)
            rgba[mask, 3] = 160

    png_bytes = _encode_png_rgba(rgba)
    data_url = f"data:image/png;base64,{base64.b64encode(png_bytes).decode('ascii')}"

    # Downsampled affine: pixel (x, y) in downsampled plane corresponds to (x*step, y*step) in original
    scale_mat = np.array([[step, 0.0, 0.0], [0.0, step, 0.0], [0.0, 0.0, 1.0]])
    downsampled_affine = np.asarray(record["pixel_to_world"]) @ scale_mat

    # Calculate corners in world space
    corners_ds = np.array([[0.0, 0.0, 1.0], [dw, 0.0, 1.0], [dw, dh, 1.0], [0.0, dh, 1.0]])
    corners_world = (downsampled_affine @ corners_ds.T)[:2].T
    bounds = [float(corners_world[:, 0].min()), float(corners_world[:, 1].min()),
              float(corners_world[:, 0].max()), float(corners_world[:, 1].max())]

    return {"asset_id": asset_id, "name": record["name"], "kind": record["kind"],
            "data_url": data_url, "shape": [h, w], "downsampled_shape": [dh, dw],
            "width": dw, "height": dh,
            "step": step, "pixel_to_world": downsampled_affine.tolist(),
            "original_pixel_to_world": record["pixel_to_world"], "bounds": bounds,
            "world_bounds": bounds}

