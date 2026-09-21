"""Immutable paired protein assays on the project's exact observation universe."""
from __future__ import annotations

from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from .store import SpatialError, _hash, _json, _now, _text

SCHEMA = "spatial-collab.protein-assay.v1"
RUN_SCHEMA = "spatial-collab.protein-region.v1"
MAX_BYTES = 32 * 1024**2


def _validate(record, project):
    if isinstance(record, dict) and record.get("schema") == "spatial-collab.protein-assay.v2":
        from .array_assays import validate
        return validate(project, record)
    if not isinstance(record, dict) or record.get("schema") != SCHEMA:
        raise SpatialError("Unsupported protein assay record.")
    payload = {k: v for k, v in record.items() if k != "assay_sha256"}
    if record.get("assay_sha256") != _hash(payload):
        raise SpatialError("Protein assay integrity check failed.")
    if record.get("assay_id") != "protein_" + _hash({k: v for k, v in payload.items() if k != "assay_id"})[:24]:
        raise SpatialError("Protein assay identity integrity check failed.")
    features, values = record.get("features"), record.get("values")
    if not isinstance(features, list) or not 1 <= len(features) <= 512 or not isinstance(values, dict):
        raise SpatialError("Invalid protein feature/matrix structure.")
    from .identity import validate_features
    validate_features({"features": features, "panel_genes": [f.get("feature_id") for f in features]})
    if len(values) != record.get("observation_count") or len(values) > 100_000 or len(values) * len(features) > 5_000_000:
        raise SpatialError("Invalid protein matrix dimensions.")
    if record.get("measurement_type") not in {"antibody_count", "intensity"}:
        raise SpatialError("Invalid protein measurement type.")
    missing = 0
    for cid, row in values.items():
        if not isinstance(cid, str) or not isinstance(row, list) or len(row) != len(features):
            raise SpatialError("Invalid protein observation row.")
        for value in row:
            if value is None:
                missing += 1
            elif type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                raise SpatialError("Invalid protein measurement.")
            elif record["measurement_type"] == "antibody_count" and value % 1:
                raise SpatialError("Protein counts must be integers.")
    if missing != record.get("missing_measurements"):
        raise SpatialError("Protein missing-measurement count differs from matrix.")
    context = project.context() if hasattr(project, "context") else project.summary()
    if any(record.get(k) != context[k] for k in ("project_id", "source_sha256")):
        raise SpatialError("Protein assay belongs to another project/source.")
    return record


def _records(project):
    if hasattr(project, "protein_assays"):
        return [_validate(r, project) for r in project.protein_assays]
    directory = project.root / "assays"
    paths = sorted(directory.glob("protein_*.json")) if directory.exists() else []
    if len(paths) > 16:
        raise SpatialError("At most 16 protein assays per project.")
    records = []
    for path in paths:
        if not path.resolve().is_relative_to(project.root.resolve()):
            raise SpatialError("Protein assay must be inside the project.")
        if path.stat().st_size > MAX_BYTES:
            raise SpatialError("Protein assay exceeds 32 MiB.")
        record = _validate(json.loads(path.read_text(encoding="utf-8")), project)
        if path.stem != record["assay_id"]:
            raise SpatialError("Protein filename differs from assay identity.")
        records.append(record)
    return records


def assay_context(project):
    directory = project.root / "assays"
    paths = sorted(directory.glob("protein_*.json")) if directory.exists() else []
    if len(paths) > 16:
        raise SpatialError("Too many protein assays.")
    # Content-addressed names are a discovery cursor; reads verify actual bytes.
    return {"count": len(paths), "cursor": _hash([p.name for p in paths])}


def list_assays(project):
    return {"assays": [{k: v for k, v in r.items() if k != "values"} for r in _records(project)],
            "shared_objects": "Same observation IDs, coordinates, selections and annotation revisions as RNA.",
            "scientific_authorization": "NOT_ESTABLISHED"}


def get_assay(project, assay_id):
    _text(assay_id, "assay_id", limit=128)
    if hasattr(project, "protein_assays"):
        candidates = [r for r in project.protein_assays if r.get("assay_id") == assay_id]
        return _validate(candidates[0], project) if candidates else _unknown_assay()
    import re
    if not re.fullmatch(r"protein_[0-9a-f]{24}", assay_id):
        raise SpatialError("Invalid protein assay ID.")
    path = project.root / "assays" / (assay_id + ".json")
    if path.is_file():
        if not path.resolve().is_relative_to(project.root.resolve()) or path.stat().st_size > MAX_BYTES:
            raise SpatialError("Invalid protein record path/size.")
        # Read and hash just the requested immutable record, never all assays.
        raw = path.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        cache = getattr(project, "_protein_cache", {})
        if cache.get(assay_id, (None,))[0] == digest and json.loads(raw).get("schema") == SCHEMA:
            return cache[assay_id][1]
        record = _validate(json.loads(raw), project)
        if record["assay_id"] == assay_id:
            project._protein_cache = {assay_id: (digest, record)}
            return record
    return _unknown_assay()


def _unknown_assay():
    raise SpatialError("Unknown protein assay; register a matching local source first.")


def register_protein(project, path, *, sample_id, source_sha256, name, measurement_type,
                     matrix="X", spatial_key="spatial", coordinate_system, units,
                     feature_id_key=None, feature_symbol_key=None, allow_partial=False,
                     registration_note, format_id="h5ad", id_column="cell_id", x_column="x", y_column="y",
                     coordinate_transform=None, coordinate_tolerance=0.0, storage="auto"):
    """CLI-only registration; no guessed cross-sample or nearest-neighbor pairing."""
    path = Path(path).resolve(strict=True)
    if storage not in {"auto", "json", "chunked"}:
        raise SpatialError("Use auto, json or chunked protein storage.")
    if path.stat().st_size > 256 * 1024**2:
        raise SpatialError("Protein input exceeds 256 MiB.")
    if measurement_type not in {"antibody_count", "intensity"} or type(allow_partial) is not bool:
        raise SpatialError("Declare antibody_count or intensity and boolean allow_partial.")
    for key, value in (("sample_id", sample_id), ("name", name), ("registration_note", registration_note)):
        _text(value, key, limit=4000)
    source_digest = hashlib.sha256(path.read_bytes()).hexdigest()
    summary = project.summary()
    meta = summary["metadata"]
    if source_sha256 != summary["source_sha256"] or sample_id != meta.get("sample_id"):
        raise SpatialError("Protein registration requires exact source_sha256 and recorded sample_id.")
    if units != meta["units"] or coordinate_system != meta["coordinate_system"]:
        raise SpatialError("Protein coordinates must use the declared project frame and units.")
    cells = project.cells()
    raw_ids = {c.get("source_cell_id", c["cell_id"]): c for c in cells}
    if len(raw_ids) != len(cells):
        raise SpatialError("Original observation IDs are ambiguous.")
    if format_id == "h5ad":
        import anndata as ad
        data = ad.read_h5ad(path, backed="r")
        try:
            if data.n_obs > 100_000 or data.n_vars > 512 or data.n_obs * data.n_vars > 5_000_000:
                raise SpatialError("Protein assay budget: 100000 observations, 512 features, 5M entries.")
            for key in ("sample_id", "slice_id", "library_id"):
                if key in data.obs and (data.obs[key].isna().any() or data.obs[key].nunique() != 1):
                    raise SpatialError("Protein input mixes or lacks sample identities.")
            if "sample_id" in data.obs and str(data.obs["sample_id"].iloc[0]) != sample_id:
                raise SpatialError("Protein source sample_id disagrees with declaration.")
            ids = list(data.obs_names)
            features = [{"feature_id": str(data.var.iloc[i][feature_id_key]) if feature_id_key else str(g),
                         "symbol": str(data.var.iloc[i][feature_symbol_key]) if feature_symbol_key else str(g)}
                        for i, g in enumerate(data.var_names)]
            if spatial_key not in data.obsm:
                raise SpatialError("Protein spatial coordinates required; no identity pairing by row order.")
            xy = np.asarray(data.obsm[spatial_key])
            raw = data.X if matrix == "X" else data.layers[matrix]
            raw = raw.to_memory() if hasattr(raw, "to_memory") else raw
            values = raw if hasattr(raw, "toarray") else np.asarray(raw)
        finally:
            data.file.close()
    elif format_id == "csv":
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            header = reader.fieldnames or []
            if len(header) != len(set(header)) or not {id_column, x_column, y_column} <= set(header):
                raise SpatialError("CSV needs unique headers and declared ID, x, y columns.")
            names = [key for key in header if key not in {id_column, x_column, y_column}]
            if not 1 <= len(names) <= 512:
                raise SpatialError("CSV requires 1..512 protein feature columns.")
            rows = []
            for row in reader:
                if len(rows) >= min(100_000, 5_000_000 // len(names)):
                    raise SpatialError("Protein CSV exceeds the observation/entry budget.")
                rows.append(row)
        ids = [r[id_column] for r in rows]
        xy = np.asarray([[float(r[x_column]), float(r[y_column])] for r in rows])
        features = [{"feature_id": g, "symbol": g} for g in names]
        values = np.asarray([[np.nan if r[g] in ("", None) else float(r[g]) for g in names] for r in rows])
    else:
        raise SpatialError("Protein reader must be h5ad or csv.")
    feature_ids = [r["feature_id"] for r in features]
    if not ids or len(ids) != len(set(ids)) or len(feature_ids) != len(set(feature_ids)):
        raise SpatialError("Protein observation and stable feature IDs must be unique and nonempty.")
    for fid in feature_ids:
        _text(fid, "protein feature ID", limit=512)
    if set(ids) - set(raw_ids):
        raise SpatialError("Protein observations include foreign IDs; no silent intersection.")
    if not allow_partial and set(ids) != set(raw_ids):
        raise SpatialError("Protein observations differ from RNA; declare allow_partial to retain missing coverage.")
    if xy.shape != (len(ids), 2) or not np.isfinite(xy).all():
        raise SpatialError("Protein coordinates must be finite N x 2.")
    expected = np.asarray([[raw_ids[i]["x"], raw_ids[i]["y"]] for i in ids])
    original_xy_sha256 = _hash(xy.tolist())
    if type(coordinate_tolerance) not in (int, float) or not math.isfinite(coordinate_tolerance) or coordinate_tolerance < 0:
        raise SpatialError("Coordinate tolerance must be nonnegative and finite in declared project units.")
    if coordinate_transform is not None:
        affine = np.asarray(coordinate_transform, dtype=float)
        if affine.shape != (3, 3) or not np.isfinite(affine).all() or not np.array_equal(affine[2], [0, 0, 1]) or abs(np.linalg.det(affine)) < 1e-15:
            raise SpatialError("Declare an invertible finite 3x3 affine coordinate transform.")
        xy = (np.column_stack([xy, np.ones(len(xy))]) @ affine.T)[:, :2]
    if not np.all(np.abs(xy - expected) <= coordinate_tolerance):
        raise SpatialError("Protein coordinates do not exactly match the corresponding original IDs.")
    use_chunks = storage == "chunked" or storage == "auto" and len(ids) * len(features) > 100_000
    if not use_chunks:
        values = values.toarray() if hasattr(values, "toarray") else np.asarray(values)
    checked_values = values.data if hasattr(values, "tocsr") else np.asarray(values)
    if np.isinf(checked_values).any() or (checked_values < 0).any():
        raise SpatialError("Protein measurements must be nonnegative or explicit missing NaN.")
    if measurement_type == "antibody_count" and np.any(checked_values[np.isfinite(checked_values)] % 1):
        raise SpatialError("Antibody counts must be integers.")
    if use_chunks:
        from .array_assays import save
        if hashlib.sha256(path.read_bytes()).hexdigest() != source_digest:
            raise SpatialError("Protein input changed while importing.")
        return save(project, {"modality": "protein", "project_id": summary["project_id"], "source_sha256": source_sha256,
             "sample_id": sample_id, "name": name, "measurement_type": measurement_type, "coordinate_system": coordinate_system,
             "units": units, "features": features, "missing_observation_count": len(cells) - len(ids), "allow_partial": allow_partial,
             "source_file": {"path": str(path), "sha256": source_digest}, "reader": format_id, "matrix": matrix,
             "registration_note": registration_note, "created_at": _now(),
             "pairing": "exact original ID plus explicit coordinate conversion/tolerance in declared frame",
             "coordinate_conversion": {"original_xy_sha256": original_xy_sha256, "transformed_xy_sha256": _hash(xy.tolist()),
                   "affine": coordinate_transform, "tolerance": coordinate_tolerance, "tolerance_units": units,
                   "maximum_absolute_residual": float(np.max(np.abs(xy - expected)))},
             "scientific_authorization": "NOT_ESTABLISHED"}, [raw_ids[cid]["cell_id"] for cid in ids], values)
    values = np.asarray(values, dtype=float)
    if values.shape != (len(ids), len(features)) or np.isinf(values).any() or (values < 0).any():
        raise SpatialError("Protein measurements must be finite nonnegative values or explicit missing NaN.")
    if measurement_type == "antibody_count" and np.any(values[np.isfinite(values)] % 1):
        raise SpatialError("Antibody counts must be integers; do not label normalized values as raw counts.")
    if hashlib.sha256(path.read_bytes()).hexdigest() != source_digest:
        raise SpatialError("Protein input changed while importing.")
    record = {"schema": SCHEMA, "modality": "protein", "project_id": summary["project_id"],
              "source_sha256": source_sha256, "sample_id": sample_id, "name": name,
              "measurement_type": measurement_type, "coordinate_system": coordinate_system, "units": units,
              "features": features, "observation_count": len(ids), "missing_observation_count": len(cells) - len(ids),
              "missing_measurements": int(np.isnan(values).sum()), "allow_partial": allow_partial,
              "source_file": {"path": str(path), "sha256": source_digest}, "reader": format_id,
              "matrix": matrix if format_id == "h5ad" else "wide_csv",
              "registration_note": registration_note, "created_at": _now(),
              "pairing": "exact original ID and exact recorded XY; declared same sample/frame",
              "coordinate_conversion": {"original_xy_sha256": original_xy_sha256, "transformed_xy_sha256": _hash(xy.tolist()),
                   "affine": coordinate_transform, "tolerance": coordinate_tolerance, "tolerance_units": units,
                   "maximum_absolute_residual": float(np.max(np.abs(xy - expected)))},
              "values": {raw_ids[cid]["cell_id"]: [None if np.isnan(v) else float(v) for v in row]
                         for cid, row in zip(ids, values)}, "scientific_authorization": "NOT_ESTABLISHED"}
    record["assay_id"] = "protein_" + _hash(record)[:24]
    record["assay_sha256"] = _hash(record)
    content = (_json(record) + "\n").encode()
    if len(content) > MAX_BYTES:
        raise SpatialError("Serialized protein assay exceeds 32 MiB.")
    folder = project.root / "assays"
    if not folder.resolve().is_relative_to(project.root.resolve()):
        raise SpatialError("Protein registration directory must be inside the project.")
    folder.mkdir(exist_ok=True)
    if len(list(folder.glob("protein_*.json"))) >= 16:
        raise SpatialError("At most 16 protein assays per project.")
    # No visible half-written record: publish a fully flushed sibling.
    temporary = folder / (record["assay_id"] + ".pending")
    with temporary.open("xb") as stream:
        stream.write(content)
    temporary.rename(folder / (record["assay_id"] + ".json"))
    return {k: v for k, v in record.items() if k != "values"}


def _feature(assay, query):
    from .identity import resolve_feature
    try:
        fid = resolve_feature({"panel_genes": [f["feature_id"] for f in assay["features"]], "features": assay["features"]}, query)
    except ValueError as exc:
        raise SpatialError(str(exc)) from exc
    return next(i for i, f in enumerate(assay["features"]) if f["feature_id"] == fid)


def _transform(values, scale, cofactor):
    if scale not in {"raw", "log1p", "asinh"} or type(cofactor) not in (float, int) or not math.isfinite(cofactor) or cofactor <= 0:
        raise SpatialError("Use raw/log1p/asinh and a positive finite cofactor.")
    if scale == "log1p":
        return np.log1p(values)
    if scale == "asinh":
        return np.arcsinh(values / cofactor)
    return values


def _summary(cells, ids, assay, feature, scale, cofactor):
    active = [c for c in cells if c["cell_id"] in ids and c["included"]]
    values = [assay["values"].get(c["cell_id"], [None] * len(assay["features"]))[feature] for c in active]
    measured = np.array([v for v in values if v is not None], dtype=float)
    transformed = _transform(measured, scale, cofactor)
    return {"included_count": len(active), "measured_count": len(measured), "missing_count": len(active) - len(measured),
            "zero_count": int(np.sum(measured == 0)), "mean_raw": float(measured.mean()) if len(measured) else None,
            "mean_transformed": float(transformed.mean()) if len(measured) else None,
            "median_raw": float(np.median(measured)) if len(measured) else None,
            "label_counts": dict(sorted(Counter(c["label"] for c in active).items()))}


def inspect_protein(project, assay_id, revision_id, feature, selection_id=None, scale="raw", cofactor=5.0,
                    bounds=None, limit=10000, offset=0):
    assay = get_assay(project, assay_id)
    index = _feature(assay, feature)
    cells = project.cells(revision_id)
    selection = project.get_selection(selection_id) if selection_id else None
    if selection_id and (not selection or selection["revision_id"] != revision_id):
        raise SpatialError("Protein selection must belong to the requested revision.")
    ids = set(selection["cell_ids"]) if selection else {c["cell_id"] for c in cells}
    stats = _summary(cells, ids, assay, index, scale, cofactor)
    if type(limit) is not int or not 1 <= limit <= 10000 or type(offset) is not int or offset < 0:
        raise SpatialError("Use limit 1..10000 and nonnegative offset.")
    if bounds is not None and (not isinstance(bounds, list) or len(bounds) != 4 or not np.isfinite(bounds).all() or bounds[0] > bounds[2] or bounds[1] > bounds[3]):
        raise SpatialError("Bounds require finite xmin,ymin,xmax,ymax.")
    visible = [c for c in cells if bounds is None or bounds[0] <= c["x"] <= bounds[2] and bounds[1] <= c["y"] <= bounds[3]]
    points = []
    for cell in visible[offset:offset + limit]:
            value = assay["values"].get(cell["cell_id"], [None] * len(assay["features"]))[index]
            points.append({"cell_id": cell["cell_id"], "x": cell["x"], "y": cell["y"], "included": cell["included"],
                           "value": value, "display_value": None if value is None else float(_transform(np.array([value]), scale, cofactor)[0])})
    return {"assay_id": assay_id, "assay_sha256": assay["assay_sha256"], "feature": assay["features"][index],
            "revision_id": revision_id, "selection_id": selection_id, "scale": scale, "cofactor": cofactor,
            "measurement_type": assay["measurement_type"], "summary": stats, "points": points,
            "view_complete": offset == 0 and len(visible) <= limit, "view_notice": "Exact viewport page; statistics use the full selected ID set, never the displayed page.",
            "bounds": bounds, "offset": offset, "total_in_view": len(visible),
            "next_offset": offset + limit if offset + limit < len(visible) else None, "display_sampling": False,
            "scientific_authorization": "NOT_ESTABLISHED"}


def compare_protein_regions(project, assay_id, base_revision, target_revision, selection_id, features,
                            background_selection_id=None, scale="raw", cofactor=5.0, rna_gene=None):
    from .analysis import _runtime_provenance
    from .identity import resolve_feature
    if not isinstance(features, list) or not 1 <= len(features) <= 32 or any(not isinstance(f, str) for f in features) or len(set(features)) != len(features):
        raise SpatialError("Choose 1..32 distinct protein features.")
    assay = get_assay(project, assay_id)
    indices = [_feature(assay, q) for q in features]
    summary = project.summary()
    before, after = project.cells(base_revision), project.cells(target_revision)
    universe = {c["cell_id"] for c in before}
    def selected(sid):
        item = project.get_selection(sid)
        if not item or item["revision_id"] not in {base_revision, target_revision}:
            raise SpatialError("Region selection must belong to a comparison revision.")
        ids = set(item["cell_ids"])
        if not ids or not ids <= universe:
            raise SpatialError("Region must have nonempty known IDs.")
        return ids
    fg = selected(selection_id)
    bg = selected(background_selection_id) if background_selection_id else universe - fg
    if fg & bg or not bg:
        raise SpatialError("Protein foreground and background must be nonempty and disjoint.")
    rna_id = resolve_feature(summary["metadata"], rna_gene) if rna_gene else None
    def evaluate(cells):
        result = {}
        rna_measurements = {}
        if rna_id:
            for cell in cells:
                if cell["cell_id"] in fg and cell["included"]:
                    total = math.fsum(cell["counts"].values())
                    count = cell["counts"].get(rna_id, 0)
                    rna_measurements[cell["cell_id"]] = (count, count / total * 10000 if total > 0 else None)
        for feature, index in zip(features, indices):
            left, right = (_summary(cells, ids, assay, index, scale, cofactor) for ids in (fg, bg))
            delta = None if None in (left["mean_transformed"], right["mean_transformed"]) else left["mean_transformed"] - right["mean_transformed"]
            paired, normalized = [], []
            if rna_id:
                for cell in cells:
                    value = assay["values"].get(cell["cell_id"], [None] * len(assay["features"]))[index]
                    if cell["cell_id"] in fg and cell["included"] and value is not None:
                        count, panel_count = rna_measurements[cell["cell_id"]]
                        paired.append((count, value))
                        if panel_count is not None:
                            normalized.append((panel_count, value))
            def correlation(pairs):
                if len(pairs) < 3 or len({x[0] for x in pairs}) < 2 or len({x[1] for x in pairs}) < 2:
                    return {"n": len(pairs), "rho": None, "status": "unknown", "reason": "insufficient_pairs_or_constant_measurement"}
                return {"n": len(pairs), "rho": float(spearmanr(np.array(pairs)[:, 0], np.array(pairs)[:, 1]).statistic), "status": "descriptive"}
            result[feature] = {"foreground": left, "background": right, "contrast": delta,
                               "paired_rna": {"rna_feature_id": rna_id, "raw": correlation(paired),
                                              "rna_panel_per_10000": correlation(normalized)} if rna_id else None}
        return result
    a, b = evaluate(before), evaluate(after)
    result = {"analysis_schema": RUN_SCHEMA, "base_revision": base_revision, "target_revision": target_revision,
              "selection_id": selection_id, "selection": {"foreground": sorted(fg), "background": sorted(bg)},
              "parameters": {"assay_id": assay_id, "features": features, "background_selection_id": background_selection_id,
                             "scale": scale, "cofactor": cofactor, "rna_gene": rna_gene}, "assay_sha256": assay["assay_sha256"],
              "before": a, "after": b,
              "comparison": {f: {"contrast_delta": None if None in (a[f]["contrast"], b[f]["contrast"]) else b[f]["contrast"] - a[f]["contrast"]} for f in features},
              "data_hashes": {"source_sha256": summary["source_sha256"], "base_revision_sha256": project.get_revision(base_revision)["revision_sha256"],
                              "target_revision_sha256": project.get_revision(target_revision)["revision_sha256"]},
              "scientific_authorization": "NOT_ESTABLISHED", "evidence_ceiling": "paired_single_slice_descriptive_protein_review",
              "provenance": {**_runtime_provenance(), "algorithm_version": "protein-region.v1"},
              "limitations": ["RNA/protein pairing is exact ID/XY plus declared sample identity, not assay specificity validation.",
                              "Missing proteins are not zero. Antibody counts and intensities are not interchangeable.",
                              "Transform is declared visualization/summary scaling, not background correction or batch integration.",
                              "Spearman correlations are descriptive paired spot/object concordance, not causation, expression equivalence or cell fractions.",
                              "RNA panel normalization is a sensitivity view; protein library/background effects remain uncontrolled.",
                              "Shared annotation/inclusion revisions affect denominators; frozen geometric membership never changes silently.",
                              "No p-values, biological replicate inference, automatic annotation or segmentation decision."]}
    return project.save_run(result)
