"""Explicit spatial-transcriptomics profiles, with bounded strict identity joins.

These adapters copy counts and centroids to a review snapshot. Images, transcript
locations, segmentation and deconvolution are outside this module's contract.
"""
from __future__ import annotations

import csv
from decimal import Decimal, InvalidOperation
import gzip
from pathlib import Path
from urllib.parse import quote

import numpy as np
from scipy import sparse

from .importers import (MAX_CELLS, MAX_GENES, MAX_NNZ, _matrix_counts, _metadata,
                        _number, _stable_sources, _text, _unique, _xenium_attributes, source_digest)

MAX_TEXT_BYTES = 512 * 1024**2
MAX_LINE_CHARS = 1_000_000
MAX_CSV_COUNT_VALUES = 20_000_000
UNANNOTATED = "Unannotated"

FORMAT_DOCUMENTATION = {
    "tenx_mtx": "https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-mex-matrices",
    "xenium": "https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/tutorials/outputs/xoa-output-understanding-outputs",
    "visium": "https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/spatial-outputs",
    "cosmx_csv": "https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/flat-file-exports/flat-files-compare.html",
    "merscope_csv": "https://vizgen.com/wp-content/uploads/2025/09/91600001_MERSCOPE_Instrument_User_Guide_RevK.pdf",
}


def _flag(value, name):
    if type(value) is not bool:
        raise ValueError(f"{name} must be an explicit boolean.")
    return value


def _scale(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("microns_per_pixel must be an explicit finite positive number.")
    factor = _number(value, "microns_per_pixel")
    if factor <= 0:
        raise ValueError("microns_per_pixel must be greater than zero.")
    return factor


def _lines(path):
    opener = gzip.open if Path(path).suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as stream:
        total = 0
        while True:
            line = stream.readline(MAX_LINE_CHARS + 1)
            if not line:
                break
            total += len(line.encode("utf-8"))
            if len(line) > MAX_LINE_CHARS or total > MAX_TEXT_BYTES:
                raise ValueError("Decompressed text exceeds the bounded reader budget; import a declared subset.")
            yield line


def _csv_rows(path, *, max_rows=MAX_CELLS):
    reader = csv.DictReader(_lines(path))
    names = reader.fieldnames
    if not names or len(names) != len(set(names)):
        raise ValueError("CSV requires unique headers.")
    for index, row in enumerate(reader):
        if index >= max_rows:
            raise ValueError("CSV observation count exceeds the alpha resource budget.")
        if None in row or any(value is None for value in row.values()):
            raise ValueError("CSV row width does not match its header.")
        yield row


def _table(path):
    rows = list(_csv_rows(path))
    if not rows:
        raise ValueError("Input CSV must contain observations.")
    return rows


def _required(rows, columns, kind):
    absent = set(columns) - set(rows[0])
    if absent:
        raise ValueError(f"{kind} missing required columns: {', '.join(sorted(absent))}.")


def _one_file(directory, names, description):
    found = [Path(directory) / name for name in names if (Path(directory) / name).is_file()]
    if len(found) != 1:
        raise ValueError(f"Expected exactly one {description}; found {len(found)}. Resolve missing or ambiguous files explicitly.")
    return found[0]


def _raw_count(value):
    try:
        number = Decimal(str(value))
    except InvalidOperation:
        raise ValueError("Raw count values must be finite nonnegative integers.") from None
    if not number.is_finite() or number < 0 or number != number.to_integral_value() or number > 2**53 - 1:
        raise ValueError("Raw count values must be nonnegative exact integers below 2**53.")
    return float(number)


def _annotations(path, ids, allow_unannotated, sources, *, key=None):
    _flag(allow_unannotated, "allow_unannotated")
    if path is None:
        if not allow_unannotated:
            raise ValueError("Annotations are required; explicitly set allow_unannotated=true to explore without labels.")
        return {identifier: UNANNOTATED for identifier in ids}
    source = Path(path).resolve(strict=True)
    sources.append(source_digest(source))
    rows = _table(source)
    _required(rows, ["label"], "Annotation CSV")
    ann_ids = _unique([key(row) if key else row.get("cell_id") for row in rows], "annotation observation ID")
    if set(ann_ids) != set(ids):
        raise ValueError("Annotation IDs must match imported observations exactly; no silent dropping or partial join.")
    output = {}
    for identifier, row in zip(ann_ids, rows):
        value = row["label"]
        output[identifier] = UNANNOTATED if allow_unannotated and not value.strip() else _text(value, "label")
    return output


def _single_slice(rows):
    for key in ("slice_id", "slide_ID", "slide_id", "sample_id", "Run_Tissue_name"):
        if key in rows[0]:
            values = {_text(row[key], key) for row in rows}
            if len(values) != 1:
                raise ValueError(f"Mixed {key} values: import one slice in one global coordinate frame.")


def _finish(destination, cells, *, name, slice_id, coordinate_system, genes, sources,
            source_kind, platform, observation_unit="cell", biological_replicates=0, parameters=None):
    from .store import Project
    if not isinstance(observation_unit, str) or observation_unit not in {"cell", "spot", "bin"}:
        raise ValueError("observation_unit must be cell, spot or bin.")
    _stable_sources(sources)
    metadata = _metadata(name, slice_id, coordinate_system, genes, sources, source_kind, biological_replicates)
    metadata.update(platform=_text(platform, "platform"), observation_unit=observation_unit,
                    label_semantics={"cell": "cell_type", "spot": "spot_annotation", "bin": "bin_annotation"}[observation_unit],
                    import_parameters=parameters or {})
    unannotated = sum(cell["label"] == UNANNOTATED for cell in cells)
    metadata["annotation_status"] = "unannotated" if unannotated == len(cells) else "partially_unannotated" if unannotated else "provided_working_labels"
    metadata["unannotated_observation_count_at_import"] = unannotated
    if unannotated:
        metadata["limitations"].append("Unannotated is a placeholder for exploration, not a biological group or inferred cell identity.")
    if observation_unit != "cell":
        metadata["limitations"].append(f"Observations are {observation_unit}s. Their labels are regional working annotations, not pure cell identities or deconvolution.")
    if source_kind in FORMAT_DOCUMENTATION:
        metadata["format_documentation"] = FORMAT_DOCUMENTATION[source_kind]
    return Project.create(destination, cells, metadata)


def _matrix_market(directory, *, require_feature_types=False):
    directory = Path(directory).resolve(strict=True)
    matrix_path = _one_file(directory, ["matrix.mtx", "matrix.mtx.gz"], "Matrix Market file")
    features_path = _one_file(directory, ["features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"], "feature TSV")
    barcode_path = _one_file(directory, ["barcodes.tsv", "barcodes.tsv.gz"], "barcode TSV")
    sources = [source_digest(path) for path in (matrix_path, features_path, barcode_path)]
    feature_rows = []
    for row in csv.reader(_lines(features_path), delimiter="\t"):
        if len(feature_rows) >= MAX_GENES or len(row) not in (2, 3):
            raise ValueError("Feature TSV must have 2 legacy gene columns or 3 typed feature columns within the gene budget.")
        feature_rows.append(row)
    if not feature_rows or len({len(row) for row in feature_rows}) != 1:
        raise ValueError("Feature TSV must be nonempty with consistent columns.")
    if require_feature_types and len(feature_rows[0]) != 3:
        raise ValueError("Xenium requires typed three-column features to distinguish RNA genes from controls.")
    for row in feature_rows:
        _text(row[1], "feature name")
        if len(row) == 3:
            _text(row[2], "feature type")
    _unique([row[0] for row in feature_rows], "feature ID")
    barcode_rows = []
    for row in csv.reader(_lines(barcode_path), delimiter="\t"):
        if len(barcode_rows) >= MAX_CELLS or len(row) != 1:
            raise ValueError("Barcode TSV must have one column within the observation budget.")
        barcode_rows.append(row[0])
    ids = _unique(barcode_rows, "matrix observation ID")
    if not ids:
        raise ValueError("Matrix requires observations.")
    keep = [i for i, row in enumerate(feature_rows) if len(row) == 2 or row[2] == "Gene Expression"]
    genes = _unique([feature_rows[i][1] for i in keep], "gene name")
    if not genes:
        raise ValueError("Matrix has no Gene Expression features.")
    lines = iter(_lines(matrix_path))
    header = next(lines, "").strip().lower().split()
    if header not in (["%%matrixmarket", "matrix", "coordinate", "integer", "general"],
                      ["%%matrixmarket", "matrix", "coordinate", "real", "general"]):
        raise ValueError("Only general coordinate Matrix Market raw integer counts are supported.")
    dimensions = next((line.split() for line in lines if line.strip() and not line.startswith("%")), [])
    try:
        if len(dimensions) != 3:
            raise ValueError
        n_genes, n_cells, nnz = (int(value) for value in dimensions)
    except ValueError:
        raise ValueError("Malformed Matrix Market dimensions.") from None
    if n_genes != len(feature_rows) or n_cells != len(ids):
        raise ValueError("Matrix dimensions must match all feature and barcode IDs exactly.")
    if not 0 <= nnz <= MAX_NNZ:
        raise ValueError("Matrix stored entries exceed the alpha resource budget.")
    rr, cc, vv = np.empty(nnz, dtype=np.int32), np.empty(nnz, dtype=np.int32), np.empty(nnz, dtype=np.float64)
    position = 0
    for line in lines:
        if not line.strip() or line.startswith("%"):
            continue
        fields = line.split()
        if len(fields) != 3 or position >= nnz:
            raise ValueError("Matrix entries do not match declared stored-entry count.")
        try:
            row, column = int(fields[0]), int(fields[1])
        except ValueError:
            raise ValueError("Matrix indices must be integers.") from None
        if not 1 <= row <= n_genes or not 1 <= column <= n_cells:
            raise ValueError("Matrix index is outside declared dimensions.")
        rr[position], cc[position], vv[position] = row - 1, column - 1, _raw_count(fields[2])
        position += 1
    if position != nnz:
        raise ValueError("Matrix is truncated relative to declared stored-entry count.")
    matrix = sparse.coo_matrix((vv, (rr, cc)), shape=(n_genes, n_cells)).tocsr()
    # Values are validated before duplicate aggregation; float64 avoids integer wrap.
    counts = _matrix_counts(matrix[keep, :].T.tocsr(), ids, genes)
    return ids, genes, counts, sources


def import_tenx_mtx(path, destination, *, positions, annotations=None, slice_id, coordinate_system,
                    units, observation_unit="cell", platform="10x_mtx", allow_unannotated=False,
                    name="10x spatial review", biological_replicates=0):
    if units != "micrometer":
        raise ValueError("Explicit units=micrometer is required; coordinates are never inferred or rescaled here.")
    ids, genes, counts, sources = _matrix_market(path)
    positions = Path(positions).resolve(strict=True)
    sources.append(source_digest(positions))
    rows = _table(positions)
    _required(rows, ["cell_id", "x", "y"], "Position CSV")
    _single_slice(rows)
    position_ids = _unique([row["cell_id"] for row in rows], "position observation ID")
    if set(position_ids) != set(ids):
        raise ValueError("Position IDs must match matrix barcodes exactly.")
    xy = {row["cell_id"]: (_number(row["x"], "x"), _number(row["y"], "y")) for row in rows}
    labels = _annotations(annotations, ids, allow_unannotated, sources)
    cells = [{"cell_id": identifier, "x": xy[identifier][0], "y": xy[identifier][1],
              "label": labels[identifier], "counts": counts[i]} for i, identifier in enumerate(ids)]
    return _finish(destination, cells, name=name, slice_id=slice_id, coordinate_system=coordinate_system,
                   genes=genes, sources=sources, source_kind="tenx_mtx", platform=platform,
                   observation_unit=observation_unit, biological_replicates=biological_replicates,
                   parameters={"feature_filter": "Gene Expression for typed features; all genes for legacy two-column format",
                               "allow_unannotated": allow_unannotated, "positions": str(positions)})


def import_xenium_source(path, destination, *, annotations=None, **options):
    from .importers import import_xenium
    return import_xenium(path, annotations, destination, **options)


def import_xenium_mtx(path, annotations, destination, *, slice_id, name="Xenium review",
                      biological_replicates=0, allow_unannotated=False):
    root = Path(path).resolve(strict=True)
    ids, genes, counts, sources = _matrix_market(root / "cell_feature_matrix", require_feature_types=True)
    positions = _one_file(root, ["cells.csv", "cells.csv.gz"], "Xenium cell summary")
    sources.append(source_digest(positions))
    rows = _table(positions)
    _required(rows, ["cell_id", "x_centroid", "y_centroid"], "Xenium cell summary")
    position_ids = _unique([row["cell_id"] for row in rows], "centroid cell ID")
    if set(ids) != set(position_ids):
        raise ValueError("Xenium centroid IDs must match matrix barcodes exactly.")
    xy = {row["cell_id"]: (_number(row["x_centroid"], "x_centroid"), _number(row["y_centroid"], "y_centroid")) for row in rows}
    labels = _annotations(annotations, ids, allow_unannotated, sources)
    cells = [{"cell_id": identifier, "x": xy[identifier][0], "y": xy[identifier][1],
              "label": labels[identifier], "counts": counts[i]} for i, identifier in enumerate(ids)]
    rows_by_id = {row["cell_id"]: row for row in rows}
    for cell in cells:
        attributes = _xenium_attributes(rows_by_id[cell["cell_id"]])
        if attributes:
            cell["attributes"] = attributes
    return _finish(destination, cells, name=name, slice_id=slice_id, coordinate_system="xenium_native_xy",
                   genes=genes, sources=sources, source_kind="xenium", platform="Xenium",
                   biological_replicates=biological_replicates,
                   parameters={"matrix_format": "10x_mtx", "allow_unannotated": allow_unannotated})


def _visium_positions(path):
    if path.name == "tissue_positions_list.csv":
        names = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
        rows = []
        for values in csv.reader(_lines(path)):
            if len(values) != len(names) or len(rows) >= MAX_CELLS:
                raise ValueError("Legacy tissue_positions_list.csv requires six columns within the observation budget.")
            rows.append(dict(zip(names, values)))
        if not rows:
            raise ValueError("Tissue positions are empty.")
    else:
        rows = _table(path)
    _required(rows, ["barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"], "Tissue positions")
    _unique([row["barcode"] for row in rows], "spot barcode")
    if any(row["in_tissue"] not in {"0", "1"} for row in rows):
        raise ValueError("in_tissue must be exactly 0 or 1.")
    return rows


def import_visium(path, destination, *, slice_id, microns_per_pixel, in_tissue_only,
                  annotations=None, matrix_dir=None, allow_unannotated=False,
                  name="Visium spot review", biological_replicates=0):
    _flag(in_tissue_only, "in_tissue_only")
    factor = _scale(microns_per_pixel)
    root = Path(path).resolve(strict=True)
    if (root / "binned_outputs").exists() or (root / "spatial" / "tissue_positions.parquet").exists():
        raise ValueError("This reader supports legacy standard Visium spots only; HD/bin data need an explicit bin-aware adapter.")
    positions = _one_file(root / "spatial", ["tissue_positions.csv", "tissue_positions_list.csv"], "legacy Visium position file")
    position_source = source_digest(positions)
    rows = _visium_positions(positions)
    selected = [row for row in rows if not in_tissue_only or row["in_tissue"] == "1"]
    if not selected:
        raise ValueError("The explicitly selected tissue subset contains no spots.")
    if matrix_dir is None:
        matrix_dir = root / ("filtered_feature_bc_matrix" if in_tissue_only else "raw_feature_bc_matrix")
        if in_tissue_only and not matrix_dir.is_dir():
            matrix_dir = root / "raw_feature_bc_matrix"
    ids, genes, counts, sources = _matrix_market(matrix_dir)
    sources.append(position_source)
    selected_ids, all_ids = {row["barcode"] for row in selected}, {row["barcode"] for row in rows}
    if set(ids) not in (selected_ids, all_ids):
        raise ValueError("Visium matrix barcodes must exactly match all positions or the explicitly requested in-tissue subset.")
    if not selected_ids <= set(ids):
        raise ValueError("Matrix lacks requested spots; choose the raw matrix for an all-spots import.")
    counts_by_id = dict(zip(ids, counts))
    labels = _annotations(annotations, selected_ids, allow_unannotated, sources)
    cells = [{"cell_id": row["barcode"],
              "x": _number(row["pxl_col_in_fullres"], "full-resolution pixel column") * factor,
              "y": _number(row["pxl_row_in_fullres"], "full-resolution pixel row") * factor,
              "label": labels[row["barcode"]], "region": "in_tissue" if row["in_tissue"] == "1" else "off_tissue",
              "counts": counts_by_id[row["barcode"]]} for row in selected]
    return _finish(destination, cells, name=name, slice_id=slice_id, coordinate_system="visium_fullres_xy_scaled_micrometer",
                   genes=genes, sources=sources, source_kind="visium", platform="Visium", observation_unit="spot",
                   biological_replicates=biological_replicates,
                   parameters={"microns_per_pixel": factor, "axis_mapping": "x=pxl_col_in_fullres; y=pxl_row_in_fullres",
                               "in_tissue_only": in_tissue_only, "source_position_count": len(rows),
                               "imported_spot_count": len(selected), "excluded_off_tissue_count": len(rows) - len(selected),
                               "matrix_directory": str(Path(matrix_dir).resolve()), "allow_unannotated": allow_unannotated})


def _wide_counts(path, *, id_columns, id_function, gene_columns=None, expected_context=None):
    rows = _csv_rows(path)
    first = next(rows, None)
    if first is None:
        raise ValueError("Counts CSV must be nonempty.")
    _required([first], id_columns, "Counts CSV")
    if gene_columns is None:
        genes = [column for column in first if column not in id_columns]
    else:
        if not isinstance(gene_columns, list):
            raise ValueError("gene_columns must be an explicit list of RNA gene columns.")
        genes = gene_columns
    genes = _unique(genes, "gene name")
    if not 0 < len(genes) <= MAX_GENES or set(genes) & set(id_columns):
        raise ValueError("Declare 1..20000 gene columns distinct from identifier columns.")
    _required([first], genes, "Counts CSV")
    counts, seen = {}, set()
    contexts = {}
    total_values, nnz = 0, 0
    from itertools import chain
    for row in chain([first], rows):
        for key in ("slice_id", "slide_ID", "slide_id", "sample_id", "Run_Tissue_name", "assay_type"):
            if key in row:
                value = _text(row[key], key)
                if key in contexts and value != contexts[key]:
                    raise ValueError(f"Mixed {key} in count CSV; import one slice/assay.")
                contexts[key] = value
                if expected_context and key in expected_context and value != expected_context[key]:
                    raise ValueError(f"Count CSV {key} does not match metadata.")
        identifier = id_function(row)
        if identifier in seen:
            raise ValueError(f"Duplicate count observation ID: {identifier}.")
        seen.add(identifier)
        total_values += len(genes)
        if total_values > MAX_CSV_COUNT_VALUES:
            raise ValueError("CSV count values exceed the streaming alpha budget; import a defined subset.")
        values = {}
        for gene in genes:
            value = _raw_count(row[gene])
            if value:
                values[gene] = value
                nnz += 1
                if nnz > MAX_NNZ:
                    raise ValueError("CSV nonzero counts exceed the alpha resource budget.")
        counts[identifier] = values
    return genes, counts, [column for column in first if column not in set(genes) | set(id_columns)]


def _cosmx_id(row):
    fov, cell = _text(row.get("fov"), "fov"), _text(row.get("cell_ID"), "cell_ID")
    return f"fov={quote(fov, safe='')};cell={quote(cell, safe='')}"


def _cosmx_file(root, supplied, patterns, standard, description):
    if supplied is not None:
        return Path(supplied).resolve(strict=True)
    candidates = {candidate.resolve() for pattern in patterns for candidate in root.glob(pattern) if candidate.is_file()}
    candidates.update((root / name).resolve() for name in standard if (root / name).is_file())
    if len(candidates) != 1:
        raise ValueError(f"Provide an explicit {description} path; matching files are missing or ambiguous.")
    return next(iter(candidates))


def import_cosmx(path, destination, *, slice_id, microns_per_pixel, gene_columns,
                 metadata_path=None, counts_path=None, annotations=None, allow_unannotated=False,
                 exclude_unassigned=False, name="CosMx cell review", biological_replicates=0):
    factor = _scale(microns_per_pixel)
    _flag(exclude_unassigned, "exclude_unassigned")
    root = Path(path).resolve(strict=True)
    metadata_path = _cosmx_file(root, metadata_path, ["*_metadata_file.csv", "*_metadata_file.csv.gz"], ["metadata.csv"], "metadata CSV")
    counts_path = _cosmx_file(root, counts_path, ["*_exprMat_file.csv", "*_exprMat_file.csv.gz"], ["exprMat.csv"], "counts CSV")
    sources = [source_digest(metadata_path), source_digest(counts_path)]
    rows = _table(metadata_path)
    _required(rows, ["fov", "cell_ID", "CenterX_global_px", "CenterY_global_px"], "CosMx global metadata (FOV-local coordinates cannot be merged)")
    _single_slice(rows)
    if "assay_type" in rows[0] and any(row["assay_type"] != "RNA" for row in rows):
        raise ValueError("CosMx reader accepts RNA assays only; mixed/protein assays need a separate adapter.")
    identifiers = _unique([_cosmx_id(row) for row in rows], "FOV/cell ID")
    genes, counts, excluded_features = _wide_counts(counts_path, id_columns=["fov", "cell_ID"],
                                                    id_function=_cosmx_id, gene_columns=gene_columns,
                                                    expected_context=rows[0])
    unassigned_metadata = {identifier for identifier, row in zip(identifiers, rows) if row["cell_ID"] == "0"}
    unassigned_counts = {identifier for identifier in counts if identifier.endswith(";cell=0")}
    if (unassigned_metadata or unassigned_counts) and not exclude_unassigned:
        raise ValueError("CosMx cell_ID=0 contains unassigned transcripts; explicitly set exclude_unassigned=true.")
    selected_rows = [row for row in rows if not exclude_unassigned or row["cell_ID"] != "0"]
    ids = [_cosmx_id(row) for row in selected_rows]
    count_ids = set(counts) - unassigned_counts if exclude_unassigned else set(counts)
    if not ids or set(ids) != count_ids:
        raise ValueError("CosMx metadata and count FOV/cell IDs must match exactly after the declared unassigned exclusion.")
    def annotation_key(row):
        return _cosmx_id(row) if "fov" in row and "cell_ID" in row else row.get("cell_id")
    labels = _annotations(annotations, ids, allow_unannotated, sources, key=annotation_key)
    cells = [{"cell_id": _cosmx_id(row), "x": _number(row["CenterX_global_px"], "global x") * factor,
              "y": _number(row["CenterY_global_px"], "global y") * factor,
              "label": labels[_cosmx_id(row)], "counts": counts[_cosmx_id(row)]} for row in selected_rows]
    return _finish(destination, cells, name=name, slice_id=slice_id, coordinate_system="cosmx_global_xy_scaled_micrometer",
                   genes=genes, sources=sources, source_kind="cosmx_csv", platform="CosMx", biological_replicates=biological_replicates,
                   parameters={"microns_per_pixel": factor, "coordinate_columns": ["CenterX_global_px", "CenterY_global_px"],
                               "identifier_rule": "fov=<percent-encoded fov>;cell=<percent-encoded cell_ID>",
                               "gene_columns": genes, "excluded_count_columns": excluded_features,
                               "exclude_unassigned": exclude_unassigned, "unassigned_metadata_rows": len(unassigned_metadata),
                               "unassigned_count_rows": len(unassigned_counts), "allow_unannotated": allow_unannotated})


def import_merscope(path, destination, *, slice_id, annotations=None, allow_unannotated=False,
                    metadata_id="EntityID", counts_id="cell", gene_columns=None,
                    name="MERSCOPE cell review", biological_replicates=0):
    root = Path(path).resolve(strict=True)
    _text(metadata_id, "metadata_id")
    _text(counts_id, "counts_id")
    metadata_path = _one_file(root, ["cell_metadata.csv", "cell_metadata.csv.gz"], "MERSCOPE cell metadata")
    counts_path = _one_file(root, ["cell_by_gene.csv", "cell_by_gene.csv.gz"], "MERSCOPE cell-by-gene counts")
    sources = [source_digest(metadata_path), source_digest(counts_path)]
    rows = _table(metadata_path)
    _required(rows, [metadata_id, "center_x", "center_y"], "MERSCOPE metadata")
    _single_slice(rows)
    ids = _unique([row[metadata_id] for row in rows], "MERSCOPE observation ID")
    genes, counts, excluded_features = _wide_counts(counts_path, id_columns=[counts_id],
        id_function=lambda row: _text(row[counts_id], "MERSCOPE count ID"), gene_columns=gene_columns)
    if set(ids) != set(counts):
        raise ValueError("MERSCOPE metadata and count IDs must match exactly.")
    labels = _annotations(annotations, ids, allow_unannotated, sources)
    cells = [{"cell_id": row[metadata_id], "x": _number(row["center_x"], "global center_x"),
              "y": _number(row["center_y"], "global center_y"), "label": labels[row[metadata_id]],
              "counts": counts[row[metadata_id]]} for row in rows]
    return _finish(destination, cells, name=name, slice_id=slice_id, coordinate_system="merscope_global_xy",
                   genes=genes, sources=sources, source_kind="merscope_csv", platform="MERSCOPE", biological_replicates=biological_replicates,
                   parameters={"metadata_id": metadata_id, "counts_id": counts_id,
                               "coordinate_columns": ["center_x", "center_y"], "input_units": "micrometer",
                               "excluded_count_columns": excluded_features, "allow_unannotated": allow_unannotated})
