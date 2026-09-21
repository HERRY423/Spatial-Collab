"""Private delivery boundaries: file identity, project isolation and real signed JWTs."""
import asyncio
import base64
import hashlib
import io
import json
import time
from types import SimpleNamespace
import zipfile

import jwt
import pytest
from cryptography.hazmat.primitives.asymmetric import rsa
from starlette.testclient import TestClient

from spatial_collab.demo import create_demo
from spatial_collab.gateway import JWTVerifier, create_gateway
from spatial_collab.private_runtime import configure, doctor
from spatial_collab.server import ToolService, create_server
from spatial_collab.store import Project, SpatialError
from spatial_collab.transfers import CHUNK_BYTES, Transfers, check_h5ad


def upload(transfers, data, filename="source.json"):
    file = transfers.begin(filename, len(data), hashlib.sha256(data).hexdigest())
    for offset in range(0, len(data), CHUNK_BYTES):
        transfers.append(file["file_id"], offset, base64.b64encode(data[offset:offset+CHUNK_BYTES]).decode())
    return transfers.complete(file["file_id"])


@pytest.fixture
def project(tmp_path):
    return create_demo(tmp_path / "project")


def test_upload_retry_integrity_quota_and_confinement(project):
    transfers = Transfers(project)
    data = b"abc123"
    file = transfers.begin("my.json", len(data), hashlib.sha256(data).hexdigest())
    args = (file["file_id"], 0, base64.b64encode(data[:3]).decode())
    assert transfers.append(*args)["offset"] == 3
    assert transfers.append(*args)["offset"] == 3
    with pytest.raises(SpatialError, match="Conflicting"):
        transfers.append(file["file_id"], 0, base64.b64encode(b"xyz").decode())
    with pytest.raises(SpatialError, match="incomplete"):
        transfers.complete(file["file_id"])
    transfers.append(file["file_id"], 3, base64.b64encode(data[3:]).decode())
    assert transfers.complete(file["file_id"])["state"] == "ready"
    other = Transfers(create_demo(project.root.parent / "other"))
    with pytest.raises(SpatialError, match="absent"):
        other.status(file["file_id"])
    with pytest.raises(SpatialError):
        transfers.status("../../source")
    with pytest.raises(SpatialError):
        transfers.begin("../bad.json", 4, "0"*64)
    with pytest.raises(SpatialError):
        transfers.begin("big.json", 129*1024**2, "0"*64)
    bad = transfers.begin("bad.json", 3, "0"*64)
    transfers.append(bad["file_id"], 0, base64.b64encode(b"bad").decode())
    with pytest.raises(SpatialError, match="checksum"):
        transfers.complete(bad["file_id"])


def test_import_new_project_and_explicit_routing(project):
    source = json.dumps({"cells": project.cells()[:2], "metadata": project.summary()["metadata"]}).encode()
    before = project.summary()
    transfers = Transfers(project)
    file = upload(transfers, source)
    imported = transfers.import_project(file["file_id"], {})
    assert project.summary() == before
    assert imported["summary"]["cell_count"] == 2
    service = ToolService(project.root)
    current = service.call("open_project", {"project_id": imported["project_id"]})
    assert current["cell_count"] == 2
    assert service.call("open_project", {})["cell_count"] == 480
    with pytest.raises(SpatialError, match="not available"):
        service.call("open_project", {"project_id": "someone-elses-project"})
    assert len(service.call("list_projects")["projects"]) == 2
    transfers.discard(file["file_id"])
    assert service.call("open_project", {"project_id": imported["project_id"]})["cell_count"] == 2
    assert next((project.root / "imported-projects").glob("*/sources/source.json")).read_bytes() == source


def test_zip_download_roundtrip_replay(project):
    from spatial_collab.replay import verify_bundle
    transfers = Transfers(project)
    descriptor = transfers.export()
    content, offset = bytearray(), 0
    while offset < descriptor["size"]:
        part = transfers.read(descriptor["file_id"], offset, 2048)
        content.extend(base64.b64decode(part["data_base64"]))
        offset = part["next_offset"]
    assert hashlib.sha256(content).hexdigest() == descriptor["sha256"]
    with zipfile.ZipFile(io.BytesIO(content)) as archive:
        assert "source_snapshot.json" in archive.namelist()
        # Only our own generated ZIP is extracted in this test, never untrusted input.
        destination = project.root.parent / "downloaded"
        archive.extractall(destination)
    verify_bundle(destination)


def test_h5ad_external_links_rejected(tmp_path):
    h5py = pytest.importorskip("h5py")
    file = tmp_path / "evil.h5ad"
    with h5py.File(file, "w") as h:
        h["outside"] = h5py.ExternalLink("private.h5", "/")
    with pytest.raises(SpatialError, match="links"):
        check_h5ad(file)


def test_h5ad_upload_preserves_original_and_ids(project, tmp_path):
    ad = pytest.importorskip("anndata")
    import numpy as np
    from scipy import sparse
    data = ad.AnnData(X=sparse.csr_matrix([[2, 0], [0, 4]]))
    data.obs_names = ["a", "b"]
    data.var_names = ["G1", "G2"]
    data.obs["cell_type"] = ["unknown", "unknown"]
    data.obsm["spatial"] = np.array([[0, 1], [2, 3]], dtype=float)
    source = tmp_path / "source.h5ad"
    data.write_h5ad(source)
    transfers = Transfers(project)
    file = upload(transfers, source.read_bytes(), "source.h5ad")
    imported = transfers.import_project(file["file_id"], {"counts_layer": "X", "slice_id": "synthetic-slice",
        "coordinate_system": "synthetic", "units": "unknown", "observation_unit": "spot"})
    assert imported["summary"]["cell_count"] == 2
    child = ToolService(project.root).call("open_project", {"project_id": imported["project_id"]})
    assert child["metadata"]["units"] == "unknown"
    transfers.discard(file["file_id"])
    assert next((project.root / "imported-projects").glob("*/sources/source.h5ad")).read_bytes() == source.read_bytes()


def test_private_runtime_empty_is_not_demo(tmp_path, monkeypatch):
    monkeypatch.setenv("SPATIAL_COLLAB_HOME", str(tmp_path / "private"))
    monkeypatch.delenv("SPATIAL_COLLAB_PROJECT", raising=False)
    settings = configure(create_empty=True)
    assert Project(settings["project"]).summary()["cell_count"] == 0
    assert doctor()["local_ready"]
    assert doctor()["chatgpt_host_verified"] is False


@pytest.fixture
def credentials():
    private = rsa.generate_private_key(public_exponent=65537, key_size=2048)
    keys = SimpleNamespace(get_signing_key_from_jwt=lambda token: SimpleNamespace(key=private.public_key()))
    def issue(**overrides):
        payload = {"sub": "researcher", "iss": "https://id.example", "aud": "https://spatial.example/projects/one/mcp",
                   "iat": int(time.time()), "exp": int(time.time())+600, "scope": "spatial:read spatial:write", **overrides}
        return jwt.encode(payload, private, algorithm="RS256")
    return keys, issue


def test_jwt_validation_rejects_foreign_tokens(credentials):
    keys, issue = credentials
    verifier = JWTVerifier("https://id.example", "https://spatial.example/projects/one/mcp", "https://id.example/jwks", {"researcher": "write"}, key_provider=keys)
    assert asyncio.run(verifier.verify_token(issue())).subject == "researcher"
    for change in ({"aud": "https://other.example"}, {"iss": "https://evil.example"}, {"sub": "intruder"},
                   {"exp": int(time.time())-1}, {"scope": "spatial:write"}, {"exp": int(time.time())+7200}):
        assert asyncio.run(verifier.verify_token(issue(**change))) is None
    assert asyncio.run(verifier.verify_token("not-a-token")) is None


def test_authenticated_http_project_and_scope_isolation(project, credentials):
    keys, issue = credentials
    other = create_demo(project.root.parent / "other")
    app = create_gateway({"public_origin": "https://spatial.example", "issuer": "https://id.example",
        "jwks_url": "https://id.example/jwks", "projects": {
            "one": {"path": str(project.root), "subjects": {"researcher": "write", "reader": "read"}},
            "two": {"path": str(other.root), "subjects": {"different-person": "write"}}}}, key_provider=keys)
    headers = {"Accept": "application/json, text/event-stream"}
    def rpc(name, arguments):
        return {"jsonrpc": "2.0", "id": 1, "method": "tools/call", "params": {"name": name, "arguments": arguments}}
    with TestClient(app, base_url="https://spatial.example") as client:
        assert client.post("/projects/one/mcp", headers=headers, json=rpc("open_project", {})).status_code == 401
        meta = client.get("/.well-known/oauth-protected-resource/projects/one/mcp").json()
        assert meta["resource"] == "https://spatial.example/projects/one/mcp"
        good = {**headers, "Authorization": "Bearer " + issue()}
        assert client.post("/projects/one/mcp", headers=good, json=rpc("open_project", {})).json()["result"]["structuredContent"]["cell_count"] == 480
        assert client.post("/projects/two/mcp", headers=good, json=rpc("open_project", {})).status_code == 401
        denied = client.post("/projects/one/mcp", headers=good, json=rpc("open_project", {"project_id": other.summary()["project_id"]}))
        assert denied.json()["result"]["isError"]
        read = {**headers, "Authorization": "Bearer " + issue(sub="reader", scope="spatial:read")}
        assert not client.post("/projects/one/mcp", headers=read, json=rpc("open_project", {})).json()["result"].get("isError")
        write = client.post("/projects/one/mcp", headers=read, json=rpc("begin_upload", {"filename": "x.json", "size": 3, "sha256": "0"*64}))
        assert write.json()["result"]["isError"]
        assert client.post("/projects/one/mcp", headers=good, content=b"x"*2_000_001).status_code == 413
        assert client.get("/healthz", headers={"Host": "evil.example"}).status_code == 400
    audit = (project.root / "access-audit.jsonl").read_text()
    assert "denied_write" in audit and "Bearer" not in audit and issue() not in audit


def test_tools_expose_explicit_project_identity(project):
    tools = asyncio.run(create_server(project.root).list_tools())
    assert all("project_id" in t.inputSchema["properties"] for t in tools)
    assert next(t for t in tools if t.name == "read_download").meta["ui"]["visibility"] == ["app"]


def test_transfer_reservations_expiry_and_cleanup(project, monkeypatch):
    import spatial_collab.transfers as module
    monkeypatch.setattr(module, "STORAGE_BYTES", 10)
    transfers = Transfers(project)
    first = transfers.begin("one.json", 6, "0" * 64)
    with pytest.raises(SpatialError, match="quota"):
        transfers.begin("two.json", 6, "0" * 64)
    with transfers.db() as db:
        db.execute("UPDATE files SET expires=0 WHERE id=?", (first["file_id"],))
    with pytest.raises(SpatialError, match="expired"):
        transfers.status(first["file_id"])
    transfers.begin("two.json", 6, "0" * 64)
    assert not transfers.path(first["file_id"]).exists()


def test_request_body_deadline():
    from spatial_collab.gateway import RequestLimits
    messages = []
    async def downstream(scope, receive, send):
        raise AssertionError("Timed-out request must not reach MCP")
    async def slow_receive():
        await asyncio.sleep(1)
        return {"type": "http.request", "body": b"", "more_body": False}
    async def send(message):
        messages.append(message)
    asyncio.run(RequestLimits(downstream, body_timeout=.01)({"type": "http"}, slow_receive, send))
    assert messages[0]["status"] == 408
