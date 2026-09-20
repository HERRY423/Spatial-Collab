from copy import deepcopy
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from spatial_collab.identity import qualify_observation_id
from spatial_collab.proteomics import (compare_protein_regions, inspect_protein, list_assays, register_protein)
from spatial_collab.replay import verify_bundle
from spatial_collab.server import ToolService
from spatial_collab.store import Project, SpatialError


@pytest.fixture
def paired(tmp_path):
    cells = [{"cell_id": qualify_observation_id("s", str(i)), "source_cell_id": str(i), "sample_id": "s",
              "x": i, "y": 0, "label": "Working", "counts": {"R": i + 1}} for i in range(4)]
    project = Project.create(tmp_path / "p", cells, {"name": "paired fixture", "sample_id": "s", "slice_id": "s",
        "source_kind": "synthetic", "units": "array_index", "coordinate_system": "array", "panel_genes": ["R"],
        "identity_scope": "sample_qualified"})
    data = ad.AnnData(np.array([[0., 10.], [2., 20.], [4., 30.], [6., 40.]]),
                     obs=pd.DataFrame(index=[str(i) for i in range(4)]), var=pd.DataFrame(index=["P", "Q"]))
    data.obsm["spatial"] = np.array([[i, 0] for i in range(4)])
    path = tmp_path / "adt.h5ad"
    data.write_h5ad(path)
    options = dict(sample_id="s", source_sha256=project.summary()["source_sha256"], name="ADT fixture",
                   measurement_type="antibody_count", coordinate_system="array", units="array_index",
                   registration_note="Known paired fixture IDs and coordinates")
    return project, path, options, data


def test_identity_shuffled_rows_share_selection_without_touching_rna(paired):
    project, path, options, data = paired
    data[[3, 1, 0, 2]].write_h5ad(path)
    before = project.summary()
    old_cursor = project.context()["cursor"]
    record = register_protein(project, path, **options)
    assert project.context(old_cursor)["changed"]
    after = project.summary()
    assert after["source_sha256"] == before["source_sha256"] and after["head_revision"] == before["head_revision"]
    selected = project.set_selection(before["head_revision"], cell_ids=[qualify_observation_id("s", "0")])
    result = inspect_protein(project, record["assay_id"], before["head_revision"], "P", selected["selection_id"])
    assert result["summary"]["mean_raw"] == 0 and result["summary"]["zero_count"] == 1
    assert {r["cell_id"]: r["value"] for r in result["points"]}[qualify_observation_id("s", "3")] == 6
    service = ToolService(project.root)
    assert service.get_overview()["protein_assays"][0]["observation_count"] == 4
    assert "values" not in service.list_assays()["assays"][0]


@pytest.mark.parametrize("change", ["sample", "source", "units", "coordinates", "fraction", "foreign_id"])
def test_wrong_pairing_and_measurement_contracts_reject(paired, change):
    project, path, options, data = paired
    if change == "sample":
        options["sample_id"] = "another"
    elif change == "source":
        options["source_sha256"] = "0" * 64
    elif change == "units":
        options["units"] = "micrometer"
    else:
        if change == "coordinates":
            data.obsm["spatial"][0, 0] = 99
        elif change == "fraction":
            data.X[0, 0] = .2
        else:
            data.obs_names = ["foreign", "1", "2", "3"]
        data.write_h5ad(path)
    with pytest.raises(SpatialError):
        register_protein(project, path, **options)
    assert list_assays(project)["assays"] == []


def test_partial_coverage_missing_intensity_not_zero_and_transforms(paired):
    project, path, options, data = paired
    data = data[:3].copy()
    data.X[0, 0] = np.nan
    data.X[1, 0] = 2.5
    data.write_h5ad(path)
    options["measurement_type"] = "intensity"
    with pytest.raises(SpatialError, match="allow_partial"):
        register_protein(project, path, **options)
    record = register_protein(project, path, allow_partial=True, **options)
    result = inspect_protein(project, record["assay_id"], project.summary()["head_revision"], "P", scale="asinh", cofactor=5.)
    assert result["summary"]["missing_count"] == 2
    assert result["summary"]["zero_count"] == 0
    assert result["summary"]["mean_transformed"] == pytest.approx(np.arcsinh(np.array([2.5, 4]) / 5).mean())


def test_frozen_roi_revision_denominators_correlations_and_offline_replay(paired):
    project, path, options, _ = paired
    record = register_protein(project, path, **options)
    base = project.summary()["head_revision"]
    ids = [qualify_observation_id("s", str(i)) for i in range(3)]
    selection = project.set_selection(base, cell_ids=ids)
    run = compare_protein_regions(project, record["assay_id"], base, base, selection["selection_id"], ["P"], rna_gene="R")
    assert run["before"]["P"]["paired_rna"]["raw"]["rho"] == pytest.approx(1)
    assert run["before"]["P"]["paired_rna"]["rna_panel_per_10000"]["rho"] is None
    one = project.set_selection(base, cell_ids=[ids[0]])
    proposal = project.propose_revision(base, one["selection_id"], {"included": False}, "Synthetic denominator test")
    target = project.apply_revision(proposal["proposal_id"], base, "TEST_ONLY", True)["revision_id"]
    changed = compare_protein_regions(project, record["assay_id"], base, target, selection["selection_id"], ["P"], rna_gene="R")
    assert changed["before"]["P"]["foreground"]["measured_count"] == 3
    assert changed["after"]["P"]["foreground"]["measured_count"] == 2
    assert changed["comparison"]["P"]["contrast_delta"] == 1
    folder = Path(project.export_bundle(changed["run_id"], compact=True)["export_path"])
    assert verify_bundle(folder)["recompute_matches"]
    assert (folder / "protein_assays.json").exists()
    source = next((project.root / "assays").glob("protein_*.json"))
    content = json.loads(source.read_text())
    content["values"][ids[0]][0] = 999
    source.write_text(json.dumps(content))
    with pytest.raises(SpatialError, match="integrity"):
        list_assays(project)
    assert verify_bundle(folder)["recompute_matches"]  # no dependence on the changed local assay


def test_wide_csv_intensity_and_duplicate_symbol_mapping(paired):
    project, path, options, data = paired
    csv_path = path.with_suffix(".csv")
    csv_path.write_text("cell_id,x,y,P\n0,0,0,1.2\n1,1,0,\n2,2,0,0\n3,3,0,4.5\n")
    record = register_protein(project, csv_path, **{**options, "format_id": "csv", "measurement_type": "intensity"})
    assert record["missing_measurements"] == 1
    other = deepcopy(options)
    data.var["stable"] = ["antibody-A", "antibody-B"]
    data.var_names = ["Repeated", "Repeated"]
    data.write_h5ad(path)
    record = register_protein(project, path, **other, feature_id_key="stable")
    with pytest.raises(SpatialError, match="Ambiguous"):
        inspect_protein(project, record["assay_id"], project.summary()["head_revision"], "Repeated")


def test_protein_home_navigation_real_service_dom(paired):
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
    import shutil
    import subprocess
    from threading import Thread
    if not shutil.which("node"):
        pytest.skip("Node runtime unavailable")
    project, path, options, _ = paired
    register_protein(project, path, **options)
    project.set_selection(project.summary()["head_revision"],
                          cell_ids=[qualify_observation_id("s", str(i)) for i in range(3)])
    service = ToolService(project.root)
    class Handler(BaseHTTPRequestHandler):
        def do_POST(self):
            args = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
            try:
                response = {"ok": True, "result": service.call(self.path.rsplit("/", 1)[-1], args)}
                self.send_response(200)
            except Exception as exc:
                response = {"ok": False, "error": str(exc)}
                self.send_response(400)
            content = json.dumps(response).encode()
            self.send_header("Content-Length", str(len(content)))
            self.end_headers()
            self.wfile.write(content)

        def log_message(self, *_args):
            pass
    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    root = Path(__file__).parents[1]
    try:
        result = subprocess.run([shutil.which("node"), str(Path(__file__).with_name("workbench_protein_smoke.cjs")),
            str(root / "src/spatial_collab/static/workbench.js"), str(root / "src/spatial_collab/static/workbench.html"),
            f"http://127.0.0.1:{server.server_port}"], capture_output=True, text=True, encoding="utf-8", timeout=60)
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)
    assert result.returncode == 0, result.stdout + result.stderr
