"""Re-index an already verified public SPOTS ROI and compare every sparse count."""
import csv
import json
from pathlib import Path
import time

import numpy as np

from spatial_collab import atlas
from spatial_collab.demo import create_demo
from spatial_collab.store import Project
from spatial_collab.workflow_inputs import load_counts


def main():
    root = Path("output/scalability/public-tmp")
    root.mkdir(parents=True, exist_ok=False)
    original = Project("projects/v05-real-data/spots-spleen-1")
    ids = json.loads(Path("output/v07/inputs-v2.json").read_text())
    record, matrix = load_counts(original, ids["roi"])
    def write(name, fields, rows):
        with (root / name).open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(fields)
            writer.writerows(rows)
    # Source IDs are preserved verbatim; qualification in the index is explicit.
    write("cells.csv", ["id", "x", "y", "label"], ((cid, *record["coordinates"][i], record["labels"][i]) for i, cid in enumerate(record["observation_ids"])))
    write("features.csv", ["id", "symbol"], zip(record["features"], record["feature_symbols"]))
    coo = matrix.tocoo()
    write("counts.csv", ["cell_id", "feature", "value"], ((record["observation_ids"][i], record["features"][j], v) for i, j, v in zip(coo.row, coo.col, coo.data)))
    project = create_demo(root / "project")
    spec = {"metadata": {"name": "Public SPOTS spleen-1 exact 300-spot ROI", "sample_id": record["sample_id"], "coordinate_system": record["provenance"]["coordinate_system"], "units": record["provenance"]["units"], "observation_unit": "spot", "platform": "SPOTS", "species": record["species"], "source_input_id": record["object_id"], "source_input_sha256": record["object_sha256"]}}
    for kind, columns in (("cells", ["id", "x", "y", "label"]), ("features", ["id", "symbol"]), ("counts", ["cell_id", "feature", "value"])):
        spec[kind] = {"path": str(root / f"{kind}.csv"), "columns": {k: k for k in columns}}
    before = time.perf_counter()
    indexed = atlas.build(project, spec)
    imported = time.perf_counter()-before
    frozen = atlas.extract(project, indexed["object_id"], indexed["bounds"])
    replayed, actual = load_counts(project, frozen["object_id"])
    order = [record["features"].index(g) for g in replayed["features"]]
    delta = matrix[:, order]-actual
    assert not delta.nnz
    np.testing.assert_array_equal(record["coordinates"], replayed["coordinates"])
    from spatial_collab.identity import qualify_observation_id
    assert replayed["observation_ids"] == [qualify_observation_id(record["sample_id"], cid) for cid in record["observation_ids"]]
    report = {"scope": "public SPOTS exact sparse counts, coordinates and producer-ID retention; no actual molecule or histology data in this dataset", "source_input_id": record["object_id"], "source_input_sha256": record["object_sha256"], "shape": list(matrix.shape), "nonzero_counts": int(matrix.nnz), "exact_sparse_difference_nnz": int(delta.nnz), "coordinates_equal": True, "identities_equal_under_declared_qualification": True, "atlas_id": indexed["object_id"], "import_seconds": imported, "project": str(project.root), "totals": indexed["totals"]}
    Path("output/scalability/public-validation.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
