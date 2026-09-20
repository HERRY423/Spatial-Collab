"""Deterministic synthetic fixture, never presented as human tissue."""
from __future__ import annotations

import random


def create_demo(destination):
    from .store import Project
    rng = random.Random(20260920)
    genes = ["CD3D", "CD3E", "LYZ", "LST1", "EPCAM", "KRT19", "COL1A1", "DCN"]
    cells = []
    # A few mixed neighborhoods intentionally contain label/marker disagreements.
    # These are planted examples, not a benchmark of annotation quality.
    for row in range(20):
        for col in range(24):
            x, y = col * 18 + rng.uniform(-4, 4), row * 18 + rng.uniform(-4, 4)
            label = ("Tumor" if col >= 15 else "Stromal" if col < 5 else
                     "T cell" if (row + col) % 3 == 0 else "Myeloid")
            true_marker_group = "Tumor" if label == "Myeloid" and 10 <= col <= 13 and 7 <= row <= 12 else label
            markers = {"T cell": ["CD3D", "CD3E"], "Myeloid": ["LYZ", "LST1"],
                       "Tumor": ["EPCAM", "KRT19"], "Stromal": ["COL1A1", "DCN"]}[true_marker_group]
            counts = {g: float(rng.randint(8, 25) if g in markers else rng.randint(0, 2)) for g in genes}
            cells.append({"cell_id": f"demo-{row:02d}-{col:02d}", "x": round(x, 3), "y": round(y, 3),
                          "label": label, "included": True, "region": "edge" if 8 <= col <= 15 else "other",
                          "counts": counts})
    metadata = {
        "name": "合成肿瘤边缘 · 协作流程演示", "slice_id": "synthetic-slice-01",
        "coordinate_system": "synthetic_xy", "units": "micrometer", "panel_genes": genes,
        "source_kind": "synthetic", "source_files": [], "biological_replicates": 0,
        "limitations": ["SYNTHETIC DEMO: not measured tissue or biological validation.",
                        "Planted marker/label disagreements illustrate review; they are not annotation ground truth.",
                        "Centroids only; no histology, segmentation masks or transcript locations.",
                        "Descriptive sensitivity only; no p-values, mechanisms or patient inference."],
        "generator": {"seed": 20260920, "version": 1},
    }
    return Project.create(destination, cells, metadata)

