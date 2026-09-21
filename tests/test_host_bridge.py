"""Exercise the actual JS bridge/file buttons against a real local tool service."""
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import shutil
import subprocess
from threading import Thread

import pytest

from spatial_collab.demo import create_demo
from spatial_collab.server import ToolService


@pytest.mark.skipif(shutil.which("node") is None, reason="Node unavailable")
def test_embedded_bridge_upload_download_and_reconnect(tmp_path):
    project = create_demo(tmp_path / "project")
    service = ToolService(project.root)
    snapshot = tmp_path / "synthetic.json"
    snapshot.write_text(json.dumps({"cells": project.cells()[:3], "metadata": project.summary()["metadata"]}), encoding="utf-8")
    class Handler(BaseHTTPRequestHandler):
        def do_POST(self):
            args = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
            try:
                with service.call_lock:
                    response = {"ok": True, "result": service.call(self.path[1:], args)}
            except Exception as exc:
                response = {"ok": False, "error": str(exc)}
            encoded = json.dumps(response).encode()
            self.send_response(200)
            self.send_header("Content-Length", str(len(encoded)))
            self.end_headers()
            self.wfile.write(encoded)
        def log_message(self, *_):
            pass
    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    root = Path(__file__).resolve().parents[1]
    try:
        result = subprocess.run([shutil.which("node"), str(Path(__file__).with_name("workbench_host_bridge.cjs")),
            str(root / "src/spatial_collab/static/workbench.js"), str(root / "src/spatial_collab/static/workbench.html"),
            str(root / "src/spatial_collab/static/transfer-ui.js"), f"http://127.0.0.1:{server.server_port}", str(snapshot)],
            capture_output=True, text=True, encoding="utf-8", timeout=60)
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
    assert result.returncode == 0, result.stdout+result.stderr
    receipt = json.loads(result.stdout)
    assert receipt.pop("real_chatgpt") is False
    assert all(receipt.values())
