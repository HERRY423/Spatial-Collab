"""Production JS against actual local ToolService; DOM fixture, not browser acceptance."""
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import shutil
import subprocess
from threading import Thread

import numpy as np
import pytest

from spatial_collab import assets, hypotheses
from spatial_collab.server import ToolService
from spatial_collab.store import Project


@pytest.mark.skipif(shutil.which("node") is None, reason="Node runtime unavailable for DOM service smoke")
def test_workbench_plan_workflow_exact_versions_and_unknown_rows(tmp_path, monkeypatch):
    gene = "gene A,1"
    metadata = {"name": "Synthetic plan UI", "slice_id": "s", "coordinate_system": "s-xy",
                "source_kind": "synthetic", "units": "micrometer", "panel_genes": [gene, "G2"] + [f"other{i}" for i in range(52)]}
    metadata["features"] = [{"feature_id": key, "symbol": "G" if key in {gene, "G2"} else key} for key in metadata["panel_genes"]]
    cells = [{"cell_id": "s", "x": 0, "y": 0, "label": "T", "counts": {gene: 1}, "attributes": {"qc": 2}},
             {"cell_id": "t", "x": 1, "y": 0, "label": "B", "counts": {gene: 1}, "attributes": {"qc": 2}},
             {"cell_id": "b", "x": 10, "y": 0, "label": "B", "counts": {}, "attributes": {"qc": 0}},
             {"cell_id": "o", "x": 20, "y": 0, "label": "O", "counts": {}}]
    project = Project.create(tmp_path / "project", cells, metadata)
    revision = project.summary()["head_revision"]
    selection = project.set_selection(revision, cell_ids=["s"])
    mask = np.zeros((1, 21), dtype=np.uint32)
    mask[0, [0, 1, 10, 20]] = [91, 7, 5, 20]
    mask_path = tmp_path / "mask.npy"
    np.save(mask_path, mask)
    assets.register_asset(project, mask_path, kind="segmentation", source_sha256=project.summary()["source_sha256"],
                           slice_id="s", coordinate_system="s-xy", units="micrometer", pixel_to_world=np.eye(3),
                           registration_note="Synthetic UI fixture, not a segmentation correctness assertion",
                           label_to_cell_id={"91": "s", "7": "s", "5": "b"})
    # Make only the large radius exceed the declared computation budget so the
    # actual result contains computed, unknown and failed rows together.
    monkeypatch.setattr(hypotheses, "MAX_DIRECTED_SOURCE_EDGES", 2)
    service = ToolService(project.root)
    calls = []

    class Handler(BaseHTTPRequestHandler):
        def do_POST(self):
            body = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
            name = self.path.rsplit("/", 1)[-1]
            calls.append({"name": name, "arguments": body})
            try:
                result = {"ok": True, "result": service.call(name, body)}
                self.send_response(200)
            except Exception as exc:
                result = {"ok": False, "error": str(exc)}
                self.send_response(400)
            encoded = json.dumps(result).encode()
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(encoded)))
            self.end_headers()
            self.wfile.write(encoded)

        def log_message(self, *_args):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    root = Path(__file__).parents[1]
    try:
        result = subprocess.run([shutil.which("node"), str(Path(__file__).with_name("workbench_plan_smoke.cjs")),
                                 str(root / "src/spatial_collab/static/workbench.js"),
                                 str(root / "src/spatial_collab/static/workbench.html"),
                                 f"http://127.0.0.1:{server.server_port}"],
                                capture_output=True, text=True, encoding="utf-8", timeout=90)
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
    assert result.returncode == 0, result.stdout + result.stderr
    assert all(json.loads(result.stdout).values())
    assert project.summary()["head_revision"] == revision
    assert len(project.summary()["revisions"]) == 1
    assert not set(c["name"] for c in calls) & {"propose_revision", "apply_revision", "revert_revision"}
    run_calls = [c["arguments"] for c in calls if c["name"] == "run_hypothesis_plan"]
    assert run_calls and all(set(c) == {"plan_id", "version"} for c in run_calls)
    records = hypotheses.list_hypothesis_plans(project)["plans"]
    assert len(records) == 1
    frozen = hypotheses.get_hypothesis_plan(project, records[0]["plan_id"], 2)
    assert frozen["spec"]["selection_id"] == selection["selection_id"]
    assert frozen["resolved"]["variants"][0]["cell_ids"] == ["t"]
    assert frozen["resolved"]["variants"][1]["cell_ids"] == ["s"]
    assert frozen["resolved"]["variants"][2]["unknown_count"] == 4
