"""Native tiled TIFF pyramids, bounded chunk decoding and explicit XY transforms."""
import base64
from io import BytesIO
from pathlib import Path

import numpy as np

from . import objects
from .assets import validate_affine
from .atlas import digest
from .store import SpatialError


def register(project, path, *, atlas_id, pixel_to_world, series=0, selectors=None, contrast=None):
    import tifffile
    atlas = objects.get(project, atlas_id, "atlas")
    affine = validate_affine(pixel_to_world).tolist()
    path = Path(path).resolve(strict=True)
    if type(series) is not int or series < 0:
        raise SpatialError("Choose a nonnegative TIFF series.")
    selectors = selectors or {}
    before = digest(path)
    levels = []
    with tifffile.TiffFile(path) as tif:
        if series >= len(tif.series):
            raise SpatialError("Unknown TIFF series.")
        for level in tif.series[series].levels:
            axes, shape = level.axes, level.shape
            if axes.count("X") != 1 or axes.count("Y") != 1:
                raise SpatialError("Image must expose explicit X and Y axes.")
            expected = {a for a in axes if a not in "YXS"}
            if set(selectors) != expected:
                raise SpatialError(f"Explicit plane indices required for axes: {sorted(expected)}")
            for a in expected:
                if type(selectors[a]) is not int or not 0 <= selectors[a] < shape[axes.index(a)]:
                    raise SpatialError("Plane selector is outside source axis.")
            if "S" in axes and (axes[-1] != "S" or shape[-1] not in (3, 4)):
                raise SpatialError("Only explicitly interleaved RGB/RGBA sample axes are supported.")
            if not level.pages[0].is_tiled:
                raise SpatialError("Streaming requires tiled TIFF. Convert strip TIFF to a tiled pyramid first.")
            chunk_pixels = level.pages[0].tilewidth * level.pages[0].tilelength
            if chunk_pixels * level.dtype.itemsize * 4 > 32*1024**2:
                raise SpatialError("Source tile exceeds bounded decoding memory.")
            levels.append({"axes": axes, "shape": list(shape), "width": shape[axes.index("X")], "height": shape[axes.index("Y")], "dtype": level.dtype.str})
        if contrast is None:
            if np.dtype(levels[0]["dtype"]) != np.dtype("uint8"):
                raise SpatialError("Non-uint8 images require an explicit common [low,high] contrast range.")
            contrast = [0, 255]
    if len(contrast) != 2 or not np.isfinite(contrast).all() or contrast[0] >= contrast[1]:
        raise SpatialError("Contrast must be a finite increasing pair.")
    if digest(path) != before:
        raise SpatialError("Image changed during registration.")
    stat = path.stat()
    return objects.put(project, "pyramid", {"schema": "spatial-collab.pyramid.v1", "atlas_id": atlas_id, "atlas_sha256": atlas["object_sha256"], "path": str(path), "sha256": before, "size": stat.st_size, "mtime_ns": stat.st_mtime_ns, "series": series, "selectors": selectors, "levels": levels, "pixel_to_world": affine, "contrast": list(contrast), "tile_size": 256, "coordinate_convention": "affine maps full-resolution pixel centers; tile corners offset -0.5", "registration_status": "researcher_declared_transform_not_validated_registration"})


def catalog(project, atlas_id):
    return {"images": [r for oid in objects.catalog(project, "pyramid") if (r := objects.get(project, oid, "pyramid"))["atlas_id"] == atlas_id]}


def tile(project, image_id, level, x, y):
    import tifffile
    import zarr
    from PIL import Image
    record = objects.get(project, image_id, "pyramid")
    if any(type(v) is not int or v < 0 for v in (level, x, y)) or level >= len(record["levels"]):
        raise SpatialError("Invalid pyramid tile indices.")
    path = Path(record["path"])
    stat = path.stat()
    if (stat.st_size, stat.st_mtime_ns) != (record["size"], record["mtime_ns"]):
        raise SpatialError("Image source changed; re-register the source.")
    info = record["levels"][level]
    size = record["tile_size"]
    x0, y0 = x*size, y*size
    x1, y1 = min(x0+size, info["width"]), min(y0+size, info["height"])
    if x0 >= x1 or y0 >= y1:
        raise SpatialError("Tile lies outside the image.")
    with tifffile.TiffFile(path) as tif:
        store = tif.aszarr(series=record["series"], level=level)
        try:
            array = zarr.open(store, mode="r")
            selection = tuple(slice(y0, y1) if a == "Y" else slice(x0, x1) if a == "X" else slice(None) if a == "S" else record["selectors"][a] for a in info["axes"])
            pixels = np.asarray(array[selection])
            remaining = [a for a in info["axes"] if a in "YXS"]
            pixels = np.transpose(pixels, [remaining.index(a) for a in ("YXS" if "S" in remaining else "YX")])
        finally:
            store.close()
    lo, hi = record["contrast"]
    pixels = np.clip((pixels.astype(float)-lo)/(hi-lo), 0, 1)
    if not np.isfinite(pixels).all():
        raise SpatialError("Image tile contains nonfinite intensities.")
    output = BytesIO()
    Image.fromarray(np.rint(pixels*255).astype("uint8")).save(output, format="PNG")
    base = record["levels"][0]
    sx, sy = base["width"]/info["width"], base["height"]/info["height"]
    # Pyramid pixel centers align at (u+.5)*scale-.5 in full-resolution space.
    corners = np.array([[x0*sx-.5, y0*sy-.5, 1], [x1*sx-.5, y0*sy-.5, 1], [x0*sx-.5, y1*sy-.5, 1]])
    world = (np.asarray(record["pixel_to_world"]) @ corners.T).T[:, :2].tolist()
    return {"image_id": image_id, "level": level, "x": x, "y": y, "width": x1-x0, "height": y1-y0, "world_corners": world, "data_url": "data:image/png;base64," + base64.b64encode(output.getvalue()).decode(), "decoded_scope": "intersecting source chunks only", "source_sha256": record["sha256"]}
