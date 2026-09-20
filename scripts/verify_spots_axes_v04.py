"""Check full HDF5 axes against SQLite JSON, without invoking the importer."""
from __future__ import annotations

import json
from pathlib import Path
import sqlite3

import h5py
from anndata.io import read_elem

from spatial_collab.identity import describe_feature_query


ROOT = Path(__file__).resolve().parents[1] / "projects" / "v04-real-data"
SOURCES = Path(r"C:\Plugin\ExoWarrant\output\sc-spatial-multiomics-20260911-01\inputs")


def main():
    rows, original_sets, qualified_sets = [], [], []
    for sample in (1, 2):
        with h5py.File(SOURCES / f"spleen{sample}_adata_RNA.h5ad", "r") as handle:
            obs, var = read_elem(handle["obs"]), read_elem(handle["var"])
            originals = set(obs.index)
            expected = [{"feature_id": gene, "symbol": symbol} for gene, symbol in zip(var["gene_ids"], var.index)]
        database = ROOT / f"spots-spleen-{sample}" / "project.sqlite3"
        with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True) as connection:
            metadata = json.loads(connection.execute("SELECT json_extract(payload,'$.metadata') FROM source WHERE id=1").fetchone()[0])
            identities = connection.execute("SELECT json_extract(value,'$.cell_id'),json_extract(value,'$.source_cell_id'),"
                "json_extract(value,'$.sample_id') FROM source,json_each(source.payload,'$.cells')").fetchall()
        assert metadata["features"] == expected
        assert set(metadata["panel_genes"]) == set(var["gene_ids"])
        assert len(identities) == len(originals)
        assert {row[1] for row in identities} == originals
        assert {row[2] for row in identities} == {f"spots-spleen-{sample}"}
        qualified = {row[0] for row in identities}
        assert len(qualified) == len(identities)
        ambiguous = var.index[var.index.duplicated()].tolist()
        resolution = describe_feature_query(metadata, ambiguous[0])
        assert resolution["status"] == "ambiguous" and len(resolution["candidate_feature_ids"]) > 1
        exact = describe_feature_query(metadata, "feature_id:" + resolution["candidate_feature_ids"][0])
        assert exact["status"] == "resolved"
        rows.append({"sample": sample, "full_feature_id_symbol_axis_exact_match": True,
                     "panel_axis_exact_match_including_zero_expression_features": True,
                     "source_barcode_axis_exact_match": True, "observation_count": len(originals),
                     "feature_count": len(expected), "real_ambiguous_symbol_example": resolution,
                     "explicit_id_resolution": exact})
        original_sets.append(originals)
        qualified_sets.append(qualified)
    shared_originals = len(original_sets[0] & original_sets[1])
    shared_qualified = len(qualified_sets[0] & qualified_sets[1])
    assert shared_originals == 2158 and shared_qualified == 0
    receipt = {"inspection": "Direct read-only HDF5 axes and SQLite JSON; no importer invocation.",
               "slices": rows, "cross_slice_shared_original_barcodes": shared_originals,
               "cross_slice_shared_qualified_observation_ids": shared_qualified,
               "evidence_ceiling": "Local identity/feature-axis preservation, not independent biological validation."}
    (ROOT / "spots-axis-identity-verification.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
