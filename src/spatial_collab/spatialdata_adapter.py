"""Bounded SpatialData Zarr reader with an explicit table→geometry→frame join.

Only points and shapes are supported. Images/labels are not loaded or presented
as segmentation evidence. Sources are never transformed or written in place.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .importers import MAX_CELLS, MAX_GENES, MAX_NNZ, _matrix_counts, _metadata, _text, _unique, source_digest
from .store import Project, SpatialError


def _members(path):
    result, total = [], 0
    for item in path.rglob("*"):
        if item.is_symlink() or getattr(item, "is_junction", lambda: False)():
            raise SpatialError("SpatialData stores must not contain symlinks or junctions.")
        if item.is_file():
            total += item.stat().st_size
            result.append(item)
            if total > 2 * 1024**3 or len(result) > 10000:
                raise SpatialError("SpatialData source exceeds the bounded 2 GiB / 10000-file import limit; export a declared subset.")
    return sorted(result)


def _preflight(root, category, path):
    # SpatialData's public reader selects categories, not individual elements.
    # Bound all arrays/parquet in the categories it will materialize, including
    # unused expression layers and geometry attributes.
    total = 0
    def walk(group):
        nonlocal total
        for _, array in group.arrays():
            size = int(np.prod(array.shape, dtype=object))
            if size > max(MAX_NNZ, MAX_CELLS, MAX_GENES):
                raise SpatialError("SpatialData array exceeds the 2 million element import budget.")
            total += max(8, array.dtype.itemsize) * size
            if total > 256 * 1024**2:
                raise SpatialError("SpatialData table metadata exceeds the 256 MiB decoded budget.")
        for _, child in group.groups():
            walk(child)
    walk(root["tables"])
    import pyarrow.parquet as pq
    rows, decoded = 0, 0
    for file in (path / category).rglob("*.parquet"):
        if not file.is_file():
            continue
        meta = pq.read_metadata(file)
        rows += meta.num_rows
        decoded += sum(meta.row_group(i).total_byte_size for i in range(meta.num_row_groups))
        if rows > MAX_CELLS or decoded > 256 * 1024**2:
            raise SpatialError("SpatialData geometry category exceeds 100000 rows / 256 MiB decoded budget.")


def import_spatialdata(path, destination, *, table_name, element_name, coordinate_system, units,
                       slice_id, observation_unit, label_key="cell_type", counts_layer="counts",
                       allow_unannotated=False, platform="SpatialData", name="SpatialData review",
                       biological_replicates=0):
    """Import one explicitly chosen linked region into calibrated 2D centroids.

    Units are a researcher declaration: a coordinate-system name such as
    'global' does not itself establish a physical scale. The stored transform is
    applied, but the declared target frame must already be in micrometres.
    """
    if units != "micrometer":
        raise SpatialError("Explicit micrometer target-frame calibration is required.")
    if observation_unit not in {"cell", "spot", "bin"}:
        raise SpatialError("Declare observation_unit as cell, spot or bin.")
    if type(allow_unannotated) is not bool:
        raise SpatialError("allow_unannotated must be a boolean.")
    for key, value in (("table_name", table_name), ("element_name", element_name), ("counts_layer", counts_layer)):
        _text(value, key)
        if "/" in value or "\\" in value or value in {".", ".."}:
            raise SpatialError(f"{key} must name a single element, not a path.")
    source = Path(path).resolve(strict=True)
    if not source.is_dir():
        raise SpatialError("SpatialData reader requires a local Zarr directory.")
    if Path(destination).resolve().is_relative_to(source):
        raise SpatialError("Destination must be outside the immutable source store.")
    import zarr
    from spatialdata import get_centroids, read_zarr
    from spatialdata.transformations import get_transformation
    paths = _members(source)
    root = zarr.open_group(str(source), mode="r")
    if "tables" not in root or table_name not in root["tables"]:
        raise SpatialError("Named SpatialData table not found.")
    categories = [kind for kind in ("points", "shapes") if kind in root and element_name in root[kind]]
    if len(categories) != 1:
        raise SpatialError("Select one Points or Shapes element; raster labels and transcript-to-cell aggregation are unsupported.")
    category = categories[0]
    _preflight(root, category, source)
    sources = [source_digest(p) for p in paths]
    sdata = read_zarr(str(source), selection=("tables", category), on_bad_files="error")
    table, element = sdata.tables[table_name], getattr(sdata, category)[element_name]
    attrs = table.uns.get("spatialdata_attrs", {})
    region_key, instance_key = attrs.get("region_key"), attrs.get("instance_key")
    if region_key not in table.obs or instance_key not in table.obs:
        raise SpatialError("Table must declare region_key and instance_key for an explicit geometry join.")
    mask = table.obs[region_key].astype(str) == element_name
    if not mask.any():
        raise SpatialError("Selected table contains no observations linked to the selected geometry element.")
    table = table[mask].copy()
    if not 0 < table.n_obs <= MAX_CELLS or not 0 < table.n_vars <= MAX_GENES:
        raise SpatialError("Selected region requires 1..100000 observations and 1..20000 features.")
    for key in ("slice_id", "sample_id", "library_id"):
        if key in table.obs and (table.obs[key].isna().any() or table.obs[key].nunique() != 1):
            raise SpatialError("Selected region contains mixed or missing sample identifiers; subset one slice first.")
    ids = _unique(table.obs_names.tolist(), "observation ID")
    genes = _unique(table.var_names.tolist(), "gene name")
    if table.obs[instance_key].isna().any():
        raise SpatialError("Missing geometry instance IDs.")
    instance_ids = table.obs[instance_key].tolist()
    if len(set(instance_ids)) != len(instance_ids):
        raise SpatialError("Duplicate geometry instance IDs in selected table region.")
    if label_key not in table.obs:
        if not allow_unannotated:
            raise SpatialError("Annotation column missing; explicitly allow_unannotated to explore without assigned identities.")
        labels = ["Unannotated"] * len(ids)
    else:
        labels = [_text(v, "annotation") for v in table.obs[label_key].tolist()]
    transform = get_transformation(element, to_coordinate_system=coordinate_system)
    matrix = np.asarray(transform.to_affine_matrix(input_axes=("x", "y"), output_axes=("x", "y")), dtype=float)
    if matrix.shape != (3, 3) or not np.isfinite(matrix).all() or np.linalg.matrix_rank(matrix[:2, :2]) != 2:
        raise SpatialError("Geometry requires a finite invertible 2D transform into the selected frame.")
    if category == "points" and "z" in element.columns:
        raise SpatialError("3D coordinates cannot be silently projected to 2D.")
    if category == "shapes" and (element.geometry.has_z.any() or element.geometry.is_empty.any()
                                  or not element.geometry.is_valid.all()):
        raise SpatialError("Shapes must be valid, nonempty 2D geometries.")
    centroids = get_centroids(element, coordinate_system=coordinate_system).compute()
    if len(centroids) > MAX_CELLS or not centroids.index.is_unique:
        raise SpatialError("Geometry has too many instances or duplicate instance IDs.")
    # Exact set equality avoids positional assignment, lossy joins and invisible
    # table/geometry subsets. The user can export an explicit subset beforehand.
    if set(instance_ids) != set(centroids.index.tolist()):
        raise SpatialError("Table instance IDs and geometry IDs differ; no implicit dropping or positional join.")
    xy = np.asarray(centroids.loc[instance_ids, ["x", "y"]], dtype=float)
    if xy.shape != (len(ids), 2) or not np.isfinite(xy).all():
        raise SpatialError("Transformed coordinates must be finite 2D centroids.")
    if counts_layer != "X" and counts_layer not in table.layers:
        raise SpatialError("Explicit raw-count layer missing; no implicit fallback to X.")
    counts = _matrix_counts(table.X if counts_layer == "X" else table.layers[counts_layer], ids, genes)
    metadata = _metadata(name, slice_id, coordinate_system, genes, sources, "spatialdata_zarr", biological_replicates)
    metadata.update(platform=platform, observation_unit=observation_unit,
                    label_semantics={"cell": "cell_type", "spot": "spot_annotation", "bin": "bin_annotation"}[observation_unit],
                    spatialdata={"table_name": table_name, "element_name": element_name, "element_type": category,
                                 "region_key": region_key, "instance_key": instance_key, "counts_layer": counts_layer,
                                 "transform_xy": matrix.tolist(), "coordinate_units_declaration": units,
                                 "input_table_observations": int(mask.shape[0]), "selected_region_observations": len(ids),
                                 "instance_ids": [str(v) for v in instance_ids]})
    metadata["limitations"].append("SpatialData geometry is reduced to centroids; images, masks, transcript coordinates and morphology were not reviewed.")
    if observation_unit != "cell":
        metadata["limitations"].append("Spot/bin annotations are regional working labels, not pure cell types or deconvolution results.")
    if paths != _members(source) or sources != [source_digest(p) for p in paths]:
        raise SpatialError("Source changed during import; retry from a frozen Zarr store.")
    cells = [{"cell_id": cid, "x": float(p[0]), "y": float(p[1]), "label": label, "counts": count}
             for cid, p, label, count in zip(ids, xy, labels, counts)]
    return Project.create(destination, cells, metadata)
