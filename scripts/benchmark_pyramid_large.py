"""A genuine 32768 x 32768 tiled pyramid generated without a whole-image array."""
import base64
from io import BytesIO
import json
from pathlib import Path
import time
import sys

import numpy as np
from PIL import Image
import psutil
import tifffile

from spatial_collab.pyramid import register, tile
from spatial_collab.store import Project


def pixels(x, y, step, height=256, width=256):
    xs = (np.arange(x, x+width, dtype=np.uint32)*step)[None, :]
    ys = (np.arange(y, y+height, dtype=np.uint32)*step)[:, None]
    return ((xs//256 + ys//256 + xs % 256//4 + ys % 256//8) % 256).astype("uint8")


def chunks(size, step):
    for y in range(0, size, 256):
        for x in range(0, size, 256):
            yield pixels(x, y, step)


def main():
    root = Path("output/scalability/pyramid-tmp")
    read_only = "--read-only" in sys.argv
    root.mkdir(parents=True, exist_ok=read_only)
    scale = json.loads(Path("output/scalability/million-optimized.json").read_text())
    project = Project(scale["project"])
    path = root / "synthetic-gigapixel.ome.tif"
    before = time.perf_counter()
    if not read_only:
        with tifffile.TiffWriter(path, ome=True, bigtiff=True) as writer:
            for level in range(8):
                size = 32768 >> level
                writer.write(chunks(size, 2**level), shape=(size, size), dtype=np.uint8, tile=(256, 256), compression="deflate", subifds=7 if level == 0 else None, subfiletype=0 if level == 0 else 1, metadata={"axes": "YX"} if level == 0 else None, maxworkers=2)
    created = time.perf_counter()-before
    r = register(project, path, atlas_id=scale["atlas_id"], pixel_to_world=[[1000/32768, 0, 0], [0, 1000/32768, 0], [0, 0, 1]])
    checks = []
    for level in range(8):
        index = (32768 >> level)//256-1
        start = time.perf_counter()
        result = tile(project, r["object_id"], level, index, index)
        seconds = time.perf_counter()-start
        actual = np.asarray(Image.open(BytesIO(base64.b64decode(result["data_url"].split(",")[1]))))
        np.testing.assert_array_equal(actual, pixels(index*256, index*256, 2**level))
        checks.append({"level": level, "tile": [index, index], "seconds": seconds, "pixels_equal": True})
    mem = psutil.Process().memory_info()
    report = {"scope": "synthetic gigapixel tiled I/O and affine control; not pathology or image registration validation", "shape": [32768, 32768], "full_resolution_pixels": 32768**2, "levels": 8, "encoded_file_bytes": path.stat().st_size, "generation_seconds": created, "peak_working_set_bytes": getattr(mem, "peak_wset", None), "rss_after_bytes": mem.rss, "whole_image_array_created": False, "image_id": r["object_id"], "atlas_id": scale["atlas_id"], "tile_checks": checks}
    report["operation"] = "fresh_process_register_and_tile_read" if read_only else "generate_register_and_read"
    if read_only:
        report["generation_seconds"] = None
    Path("output/scalability/gigapixel-read.json" if read_only else "output/scalability/gigapixel-validation.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
