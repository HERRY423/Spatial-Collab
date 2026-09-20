"""Local real-format native napari acceptance; never a biological adjudication.

Read a matching Xenium tiny output bundle, explicitly derive the cell-mask
mapping from its label_id/cell_id CSV columns, preserve original bytes, and
register a bounded DAPI plane plus a lossless NPY copy of masks/1. No annotation
proposal, revision, researcher identity or scientific approval is created.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np

from spatial_collab import assets
from spatial_collab.importers import import_xenium
from spatial_collab.napari_adapter import (NapariController, add_asset_to_viewer, capture_asset_layer,
                                           capture_points, make_dock_widget)
from spatial_collab.store import SpatialError


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SpatialError("Use a fresh output directory; preserve earlier acceptance records.")
    import napari
    import tifffile
    import zarr
    args.output.mkdir(parents=True)
    outs = args.outs.resolve(strict=True)
    experiment_path = outs / "experiment.xenium"
    experiment = json.loads(experiment_path.read_text())
    if experiment.get("run_name") != "Xenium_4_Tiny_Test_Data" or experiment.get("region_name") != "hKidney_sample":
        raise SpatialError("This bounded acceptance is specifically for the declared kidney tiny-format bundle.")
    project = import_xenium(outs, None, args.output / "kidney-project", allow_unannotated=True,
                            slice_id="Xenium_4_Tiny_Test_Data:hKidney_sample:Slide_1",
                            name="Kidney tiny real-format acceptance; no biological conclusions")
    before = project.context()
    controller = NapariController(project.root)
    binding = {key: controller.binding[key] for key in ("source_sha256", "slice_id", "coordinate_system", "units")}
    image_path = outs / "morphology_focus" / "ch0000_dapi.ome.tif"
    with tifffile.TiffFile(image_path) as image:
        pixels = ET.fromstring(image.ome_metadata).find(".//{*}Pixels")
        sx, sy = float(pixels.attrib["PhysicalSizeX"]), float(pixels.attrib["PhysicalSizeY"])
        if pixels.attrib["PhysicalSizeXUnit"] != "µm" or pixels.attrib["PhysicalSizeYUnit"] != "µm":
            raise SpatialError("OME physical units are not explicitly micrometers.")
    if not np.allclose([sx, sy], experiment["pixel_size"], rtol=1e-7, atol=0):
        raise SpatialError("OME and experiment pixel scales disagree.")
    image_record = assets.register_asset(project, image_path, kind="image", **binding, tiff_page=0,
        pixel_to_world=[[sx, 0, 0], [0, sy, 0], [0, 0, 1]], name="Kidney tiny · DAPI focus (source page 0)",
        origin_sources=[assets.file_digest(experiment_path)],
        registration_note="Matching kidney tiny bundle: experiment.xenium points to this morphology_focus file; "
                          "OME PhysicalSizeX/Y and experiment.pixel_size both 0.2125 micrometers. Native origin, "
                          "no crop or registration estimated; real-format check only.")
    mask_path = outs / "cells.zarr.zip"
    boundaries = outs / "cell_boundaries.csv.gz"
    mapping = {}
    with gzip.open(boundaries, "rt") as stream:
        for row in csv.DictReader(stream):
            label, cid = row["label_id"], row["cell_id"]
            if label in mapping and mapping[label] != cid:
                raise SpatialError("Native label_id maps to conflicting exact cell IDs.")
            mapping[label] = cid
    with zarr.storage.ZipStore(str(mask_path), mode="r") as store:
        group = zarr.open_group(store=store, mode="r")
        native = group["masks/1"]
        assets.validate_array_budget(native.shape, native.dtype)
        physical_to_image = group["masks/homogeneous_transform"][:]
        if physical_to_image.shape != (4, 4) or not np.allclose(physical_to_image[2:], [[0, 0, 1, 0], [0, 0, 0, 1]]):
            raise SpatialError("Unsupported native transform; refuse a guessed 2D reduction.")
        matrix = physical_to_image[np.ix_([0, 1, 3], [0, 1, 3])]
        image_to_physical = np.linalg.inv(assets.validate_affine(matrix))
        if not np.allclose(image_to_physical, [[sx, 0, 0], [0, sy, 0], [0, 0, 1]], rtol=1e-6, atol=1e-8):
            raise SpatialError("Mask transform and image scale/origin disagree.")
        mask = native[:]
    mask_record = assets.register_array_snapshot(project, mask, kind="segmentation", **binding,
        pixel_to_world=image_to_physical, label_to_cell_id=mapping, name="Kidney tiny · native cell mask",
        origin_sources=[assets.file_digest(p) for p in (mask_path, boundaries, experiment_path)],
        registration_note="Lossless bounded masks/1 array from this exact cells.zarr.zip. Explicit positive "
                          "label_id -> cell_id pairs read from cell_boundaries.csv.gz, never inferred from row order. "
                          "Physical-to-image homogeneous_transform inverted explicitly; no biological approval.")
    del mask
    checks = []
    viewer = napari.Viewer(show=False)
    try:
        unrelated = viewer.add_points([[0, 0]], name="unrelated acceptance layer")
        dock = make_dock_widget(viewer)
        viewer.window.add_dock_widget(dock)
        dock.opened(controller)
        for record in (image_record, mask_record):
            data, kwargs, kind = controller.asset_layer(record["asset_id"])
            layer = add_asset_to_viewer(viewer, (data, kwargs, kind))
            dock.asset_layers[record["asset_id"]] = layer
            dock.points.selected_data = set(range(min(100, len(controller.coordinates))))
            inspection = controller.inspect_asset(record["asset_id"], capture_asset_layer(layer), capture_points(dock.points))
            checks.append({"asset": record, "native_layer_kind": kind, "inspection": inspection,
                           "native_data_to_world_yx": list(layer.data_to_world((100., 200.))),
                           "source_snapshot_read_only": not data.flags.writeable,
                           "detached_render_buffer_writeable": layer.data.flags.writeable,
                           "editable": layer.editable if kind == "labels" else None})
            if kind == "labels":
                layer.translate = (10, 0)
                try:
                    controller.inspect_asset(record["asset_id"], capture_asset_layer(layer), capture_points(dock.points))
                except SpatialError as exc:
                    checks[-1]["moved_layer_rejected"] = "transformed" in str(exc)
                else:
                    raise AssertionError("Moved visible segmentation unexpectedly accepted")
        assert unrelated in viewer.layers
    finally:
        viewer.close()
    after = project.context()
    assert before == after
    all_checks = []
    ids = [cell["cell_id"] for cell in project.cells()]
    for start in range(0, len(ids), 100):
        all_checks.extend(assets.inspect_asset(project, mask_record["asset_id"], ids[start:start + 100])["observations"])
    receipt = {"acceptance": "REAL_FORMAT_NATIVE_HIDDEN_VIEWER_ONLY", "scientific_authorization": "NOT_ESTABLISHED",
               "biological_conclusion": "NOT_EVALUATED", "visible_gui_acceptance": False,
               "source_bundle": str(outs), "experiment": {key: experiment[key] for key in
                 ("run_name", "region_name", "slide_id", "experiment_uuid", "analysis_uuid", "pixel_size")},
               "project": str(project.root), "head_revision": controller.revision_id,
               "source_sha256": controller.binding["source_sha256"], "observation_count": len(ids),
               "exact_mapping_source": assets.file_digest(boundaries), "mask_mapping_count": len(mapping),
               "native_transform_physical_to_image": physical_to_image.tolist(),
               "centroid_correspondence_counts": dict(Counter(c["correspondence"] for c in all_checks)),
               "missing_or_mismatched_centroid_examples": [c for c in all_checks if c["correspondence"] != "same_declared_cell"][:20],
               "assets": checks, "shared_state_before": before, "shared_state_after": after,
               "shared_state_unchanged": before == after, "source_files_modified": False,
               "limitations": ["Tiny input tests actual formats and host coordinates, not representative biological performance.",
                               "Centroid/mask correspondence is descriptive and is not a segmentation accuracy metric.",
                               "No merged-cell, doublet or true-coexpression verdict was generated."]}
    (args.output / "native-assets-receipt.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"receipt": str(args.output / "native-assets-receipt.json"), "cells": len(ids),
                      "correspondence": receipt["centroid_correspondence_counts"], "unchanged": before == after}))


if __name__ == "__main__":
    main()
