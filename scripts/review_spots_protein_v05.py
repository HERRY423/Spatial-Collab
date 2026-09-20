"""Full paired SPOTS ADT/RNA check against original sparse matrices; no inferred annotations."""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
from scipy.stats import spearmanr

from spatial_collab.identity import resolve_feature
from spatial_collab.proteomics import compare_protein_regions, get_assay, inspect_protein, register_protein
from spatial_collab.replay import verify_bundle
from spatial_collab.store import Project

ROOT = Path(__file__).resolve().parents[1]
INPUT = Path(r"C:\Plugin\ExoWarrant\output\sc-spatial-multiomics-20260911-01\inputs")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sample", type=int, choices=[1, 2])
    args = parser.parse_args()
    sample = f"spots-spleen-{args.sample}"
    output = ROOT / "projects/v05-real-data"
    output.mkdir(exist_ok=True)
    old = Project(ROOT / "projects/v04-real-data" / sample)
    original = old.summary()
    destination = output / sample
    if destination.exists():
        raise RuntimeError("Destination exists; preserve this run and use its receipt, do not overwrite.")
    project = Project.create(destination, old.cells(), original["metadata"])
    assert project.summary()["source_sha256"] == original["source_sha256"]
    del old
    gc.collect()
    path = INPUT / f"spleen{args.sample}_adata_ADT.h5ad"
    source_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    assay_meta = register_protein(project, path, sample_id=sample, source_sha256=original["source_sha256"],
        name=f"SPOTS spleen {args.sample} · 21 antibody channels", measurement_type="antibody_count",
        coordinate_system=original["metadata"]["coordinate_system"], units=original["metadata"]["units"],
        registration_note="Previously downloaded paired SPOTS spleen RNA/ADT; exact original barcode and XY verified; no physical scale inferred.")
    assay = get_assay(project, assay_meta["assay_id"])
    data = ad.read_h5ad(path)
    matrix = data.X.toarray() if hasattr(data.X, "toarray") else np.asarray(data.X)
    cells = project.cells()
    raw_ids = {c["source_cell_id"]: c for c in cells}
    actual = np.asarray([assay["values"][raw_ids[i]["cell_id"]] for i in data.obs_names])
    assert np.array_equal(actual, matrix)
    assert [f["feature_id"] for f in assay["features"]] == list(data.var_names)
    xy = np.asarray(data.obsm["spatial"])
    foreground = xy[:, 0] <= np.median(xy[:, 0])
    revision = project.summary()["head_revision"]
    selection = project.set_selection(revision, cell_ids=[raw_ids[c]["cell_id"] for c in data.obs_names[foreground]],
        name="Geometry-only x <= median source x; frozen before protein comparison")
    run = compare_protein_regions(project, assay["assay_id"], revision, revision, selection["selection_id"],
                                  list(data.var_names), rna_gene="Cd79a")
    rna_id = resolve_feature(original["metadata"], "Cd79a")
    # Independent check directly from the original sparse RNA file, not plugin summaries.
    rna = ad.read_h5ad(INPUT / f"spleen{args.sample}_adata_RNA.h5ad")
    assert list(rna.obs_names) == list(data.obs_names)
    with __import__("h5py").File(INPUT / f"spleen{args.sample}_adata_RNA.h5ad", "r") as handle:
        from anndata.io import read_elem
        var = read_elem(handle["var"])
    candidates = [i for i in range(len(var)) if str(var.index[i]) == rna_id or
                  any(str(var.iloc[i][column]) == rna_id for column in var.columns)]
    assert len(candidates) == 1, (rna_id, candidates)
    rna_values = np.asarray(rna.X[:, candidates[0]].toarray()).ravel()
    libraries = np.asarray(rna.X.sum(axis=1)).ravel()
    checks = []
    for i, feature in enumerate(data.var_names):
        item = run["before"][feature]
        for key, mask in (("foreground", foreground), ("background", ~foreground)):
            assert np.isclose(item[key]["mean_raw"], matrix[mask, i].mean())
            assert item[key]["zero_count"] == int(np.sum(matrix[mask, i] == 0))
            assert item[key]["missing_count"] == 0
        rho = float(spearmanr(rna_values[foreground], matrix[foreground, i]).statistic)
        normalized = float(spearmanr(rna_values[foreground] / libraries[foreground] * 10000, matrix[foreground, i]).statistic)
        assert np.isclose(item["paired_rna"]["raw"]["rho"], rho)
        assert np.isclose(item["paired_rna"]["rna_panel_per_10000"]["rho"], normalized)
        assert run["comparison"][feature]["contrast_delta"] == 0
        checks.append({"protein": feature, "foreground_mean": item["foreground"]["mean_raw"],
                       "background_mean": item["background"]["mean_raw"], "Cd79a_raw_rho": rho,
                       "Cd79a_rna_panel_per_10000_rho": normalized})
    inspected = inspect_protein(project, assay["assay_id"], revision, "CD19", selection["selection_id"])
    assert len(inspected["points"]) == len(cells)
    assert project.summary()["source_sha256"] == original["source_sha256"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == source_hash
    bundle = project.export_bundle(run["run_id"], compact=True)
    del rna, cells, raw_ids, data, matrix, actual, project
    gc.collect()
    verification = verify_bundle(bundle["export_path"])
    assert verification["integrity_verified"] and verification["recompute_matches"]
    receipt = {"sample": sample, "observations": assay["observation_count"], "protein_features": len(assay["features"]),
               "adt_counts": sum(sum(v) for v in assay["values"].values()), "full_matrix_exact_match": True,
               "foreground_count": selection["cell_count"], "units": "array_index", "rna_source_unchanged": True,
               "assay": assay_meta, "run_id": run["run_id"], "independent_source_checks": checks,
               "export": bundle, "replay": verification,
               "limits": "All 21 channels retained, exploratory single-slice descriptive correlations. No expert annotations, protein background correction, cross-sample model, physical distances or biological replicate inference."}
    receipt_path = output / f"{sample}-protein-review.json"
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"receipt": str(receipt_path), "observations": receipt["observations"], "features": 21,
                      "recompute_matches": verification["recompute_matches"], "channels": checks}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
