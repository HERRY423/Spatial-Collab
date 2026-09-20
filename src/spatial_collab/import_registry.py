"""Declarative, lazy reader discovery; source probing never imports a reader."""
from __future__ import annotations

from copy import deepcopy
import csv
import gzip
from importlib import import_module
import json
from pathlib import Path
from typing import Callable


_READERS: dict[str, tuple[dict, Callable | str]] = {}


def register_reader(metadata: dict, reader: Callable | str, *, replace: bool = False) -> None:
    """Register trusted local Python code, never code named by an input dataset.

    ``reader`` is a callable or a lazy ``module:function`` import target. Reader
    call signatures are ``reader(path, destination, **options) -> Project``.
    """
    if not isinstance(metadata, dict):
        raise ValueError("Reader metadata must be an object.")
    clean = json.loads(json.dumps(metadata, allow_nan=False))
    if not isinstance(clean.get("id"), str) or not clean["id"].strip():
        raise ValueError("Reader metadata requires a nonempty id.")
    for key in ("description", "platform", "observation_units", "required_inputs", "required_options", "dependencies", "ceilings"):
        if key not in clean:
            raise ValueError(f"Reader metadata requires {key}.")
    if not callable(reader) and not (isinstance(reader, str) and reader.count(":") == 1):
        raise ValueError("reader must be a trusted callable or module:function target.")
    if type(replace) is not bool:
        raise ValueError("replace must be a boolean.")
    if clean["id"] in _READERS and not replace:
        raise ValueError(f"Reader {clean['id']} is already registered.")
    _READERS[clean["id"]] = (clean, reader)


def list_formats() -> list[dict]:
    """Return capabilities without importing optional scientific dependencies."""
    return [deepcopy(_READERS[key][0]) for key in sorted(_READERS)]


def import_source(format_id: str, path, destination, **options):
    """Import only the explicitly chosen reader; never auto-select on a probe."""
    if not isinstance(format_id, str) or format_id not in _READERS:
        raise ValueError(f"Unknown format {format_id!r}; inspect list_formats().")
    _, reader = _READERS[format_id]
    if isinstance(reader, str):
        module, function = reader.split(":")
        try:
            reader = getattr(import_module(module), function)
        except ModuleNotFoundError as exc:
            raise ValueError(f"Reader {format_id} needs an optional dependency: {exc.name}.") from exc
    return reader(path, destination, **options)


def _header(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8-sig", newline="") as stream:
            line = stream.readline(1_000_001)
        return next(csv.reader([line])) if len(line) <= 1_000_000 else []
    except (OSError, UnicodeError, csv.Error, StopIteration):
        return []


def probe_source(path) -> dict:
    """Inspect filenames and at most one CSV header; no count matrix loading."""
    source = Path(path).resolve(strict=True)
    found: dict[str, dict] = {}

    def add(reader_id: str, evidence: list[str], missing: list[str] | None = None):
        found[reader_id] = {"format_id": reader_id, "evidence": evidence,
                            "missing_inputs": missing or [],
                            "required_options": deepcopy(_READERS[reader_id][0]["required_options"])}

    if source.is_file():
        if source.suffix == ".h5ad":
            add("anndata_h5ad", [".h5ad extension; dataset contents not loaded"])
        if source.name.endswith((".csv", ".csv.gz")):
            header = _header(source)
            if {"CenterX_global_px", "CenterY_global_px", "fov", "cell_ID"} <= set(header):
                add("cosmx_csv", ["CosMx global-coordinate metadata header"], ["count CSV", "annotations or allow_unannotated"])
            if {"center_x", "center_y"} <= set(header):
                add("merscope_csv", ["MERSCOPE-like centroid header; explicitly map ID fields"], ["cell_by_gene.csv", "annotations or allow_unannotated"])
    else:
        names = {entry.name for entry in source.iterdir()}
        if "cells.csv" in names or "cells.csv.gz" in names:
            add("xenium", ["cells.csv[.gz] present"], [] if "cell_feature_matrix.h5" in names or "cell_feature_matrix" in names else ["cell_feature_matrix.h5 or unpacked cell_feature_matrix directory"])
        if "matrix.mtx" in names or "matrix.mtx.gz" in names:
            add("tenx_mtx", ["Matrix Market filename present"], ["positions CSV", "annotations or allow_unannotated"])
        spatial = source / "spatial"
        if spatial.is_dir() and any((spatial / n).is_file() for n in ("tissue_positions.csv", "tissue_positions_list.csv")):
            add("visium", ["legacy Visium tissue positions file present"])
        if {"cell_metadata.csv", "cell_by_gene.csv"} <= names:
            add("merscope_csv", ["MERSCOPE cell_metadata.csv and cell_by_gene.csv present"])
        if any("metadata_file.csv" in n for n in names) or {"metadata.csv", "exprMat.csv"} <= names:
            add("cosmx_csv", ["CosMx-style flat-file filenames; headers and frame still require validation"])
        if source.suffix == ".zarr" or {".zgroup", ".zattrs"} <= names or "zarr.json" in names:
            add("spatialdata_zarr", ["Zarr container marker; SpatialData schema not yet validated"])
    return {"path": str(source), "candidates": [found[key] for key in sorted(found)],
            "ambiguous": len(found) > 1, "automatic_selection": False,
            "inspection_scope": "filenames and bounded CSV headers only; no scientific validity or import success implied"}


def _register_builtin(identifier, description, platform, units, inputs, options, dependencies, target, ceilings):
    register_reader({"id": identifier, "description": description, "platform": platform,
                     "observation_units": units, "required_inputs": inputs, "required_options": options,
                     "dependencies": dependencies, "ceilings": ceilings}, target)


_COMMON = ["Single slice and explicit coordinate frame; no segmentation or image import.",
           "Raw integer RNA counts only; strict ID alignment; no automatic annotation.",
           "Bounded alpha imports: 100000 observations, 20000 features, 2000000 stored nonzero values."]
_register_builtin("anndata_h5ad", "AnnData raw counts plus N x 2 spatial coordinates", "declared_by_user", ["cell", "spot", "bin"],
                  [".h5ad"], ["slice_id", "coordinate_system", "units=micrometer|pixel|array_index|unknown"], ["anndata", "h5py"],
                  "spatial_collab.importers:import_h5ad", _COMMON + [
                      "Explicit sample_id qualifies observation IDs and preserves source_cell_id; optional feature_id_key selects unique stable IDs, feature_symbol_key selects symbols. Duplicate symbols stay separate and queries report ambiguity.",
                      "Defaults remain 20000 features/2000000 entries; explicit max_features and max_nnz may raise bounded limits to at most 100000 features/10000000 entries. Full expression is retained within the declared budget.",
                      "Uncalibrated pixel/array_index/unknown coordinates permit expression and annotation exploration; physical-radius analyses remain unavailable."])
_register_builtin("xenium", "Xenium cells.csv plus H5 or unpacked MTX; optional declared H5 spatial window", "Xenium", ["cell"],
                  ["cells.csv[.gz]", "cell_feature_matrix.h5 or cell_feature_matrix/", "annotations unless allow_unannotated"], ["slice_id"], ["h5py for H5"],
                  "spatial_collab.format_readers:import_xenium_source", _COMMON + ["Optional bounds=[xmin,ymin,xmax,ymax] reads all H5 cells in a native-micrometer window; no sampling or RNA-gene dropping. Full source identity scan limited to 2000000 cells; selected-cell/count budgets remain unchanged.", "Window projects have no observations outside their import bounds. Instrument QC attributes are reported values, not inferred confidence."])
_register_builtin("tenx_mtx", "10x gene-by-observation MTX with explicit centroid and annotation sidecars", "declared_by_user", ["cell", "spot", "bin"],
                  ["matrix.mtx[.gz]", "features.tsv[.gz] or genes.tsv[.gz]", "barcodes.tsv[.gz]", "positions", "annotations unless allow_unannotated"],
                  ["positions", "slice_id", "coordinate_system", "units=micrometer"], ["scipy", "numpy"],
                  "spatial_collab.format_readers:import_tenx_mtx", _COMMON)
_register_builtin("visium", "Legacy standard Visium spots, never HD or segmented cells", "Visium", ["spot"],
                  ["10x MTX directory", "spatial/tissue_positions.csv or tissue_positions_list.csv", "annotations unless allow_unannotated"],
                  ["slice_id", "microns_per_pixel", "in_tissue_only"], ["scipy", "numpy"],
                  "spatial_collab.format_readers:import_visium", _COMMON + ["Explicit full-resolution pixel calibration; spot labels are regional annotations, not pure cell identities."])
_register_builtin("cosmx_csv", "CosMx RNA flat files with global pixel centroids and composite FOV/cell identity", "CosMx", ["cell"],
                  ["metadata CSV", "wide count CSV", "annotations unless allow_unannotated"],
                  ["slice_id", "microns_per_pixel", "gene_columns"], ["numpy", "scipy"],
                  "spatial_collab.format_readers:import_cosmx", _COMMON + ["Global CenterX/Y_global_px only; local FOV coordinates rejected; RNA gene columns explicitly declared."])
_register_builtin("merscope_csv", "MERSCOPE cell_metadata.csv plus cell_by_gene.csv in native global micrometres", "MERSCOPE", ["cell"],
                  ["cell_metadata.csv", "cell_by_gene.csv", "annotations unless allow_unannotated"], ["slice_id"], ["numpy", "scipy"],
                  "spatial_collab.format_readers:import_merscope", _COMMON)
_register_builtin("spatialdata_zarr", "SpatialData Zarr with explicit table and coordinate transform", "declared_by_user", ["cell", "spot", "bin"],
                  ["SpatialData .zarr store"], ["slice_id", "table_name", "element_name", "coordinate_system", "units=micrometer", "observation_unit"], ["spatialdata"],
                  "spatial_collab.spatialdata_adapter:import_spatialdata", _COMMON + ["Explicit region/table linkage and supported calibrated transforms only."])
