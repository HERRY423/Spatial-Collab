import asyncio
import base64
import hashlib
import re

import pytest
from starlette.testclient import TestClient

from spatial_collab.server import MAX_BODY_BYTES, MAX_VIEW_CELLS, UI_URI, ToolService, create_app, create_server
from spatial_collab.store import Project, SpatialError


@pytest.fixture
def project_root(tmp_path):
    root = tmp_path / "project"
    Project.create(root, [
        {"cell_id": "t1", "x": 0, "y": 0, "label": "T cell", "counts": {"CD3D": 4}},
        {"cell_id": "m1", "x": 1, "y": 0, "label": "Myeloid", "counts": {"LST1": 5}},
        {"cell_id": "e1", "x": 100, "y": 100, "label": "Epithelial", "counts": {}},
    ], {"name": "Synthetic review", "slice_id": "synthetic-slice", "coordinate_system": "slice_um",
        "units": "micrometer", "panel_genes": ["CD3D", "LST1"], "source_kind": "synthetic",
        "source_files": [], "biological_replicates": 0, "limitations": ["Synthetic fixture"]})
    return root


def test_tools_and_shared_ui_resource(project_root):
    server = create_server(project_root)
    tools = asyncio.run(server.list_tools())
    assert {"get_context", "query_observations", "get_proposal", "get_capabilities", "run_region_comparison"} <= {tool.name for tool in tools}
    by_name = {t.name: t for t in tools}
    assert by_name["open_project"].meta["ui"]["resourceUri"] == UI_URI
    assert by_name["open_project"].annotations.readOnlyHint is True
    assert by_name["apply_revision"].annotations.readOnlyHint is False
    assert "confirmation" in by_name["apply_revision"].inputSchema["required"]
    assert all("path" not in t.inputSchema.get("properties", {}) for t in tools)
    resources = asyncio.run(server.read_resource(UI_URI))
    resource = list(resources)[0]
    assert resource.mime_type == "text/html;profile=mcp-app"
    assert resource.meta["ui"]["csp"]["connectDomains"] == []
    assert "ui/initialize" in resource.content
    assert "ui/update-model-context" in resource.content
    assert "script src=" not in resource.content
    assert "/*__SCRIPT__*/" not in resource.content
    assert 'content="__CSRF_TOKEN__"' not in resource.content


def test_mcp_strict_errors_do_not_coerce_confirmation(project_root):
    server = create_server(project_root)
    summary = asyncio.run(server.call_tool("open_project", {}))
    assert not summary.isError
    result = asyncio.run(server.call_tool("apply_revision", {
        "proposal_id": "none", "expected_revision": "none", "reviewer": "Researcher", "confirmation": "true",
    }))
    assert result.isError
    assert "boolean" in result.structuredContent["error"]
    result = asyncio.run(server.call_tool("open_project", {"path": "C:/other"}))
    assert result.isError
    assert "Unexpected keyword" in result.structuredContent["error"]
    result = asyncio.run(server.call_tool("invented", {}))
    assert result.isError


def test_shared_api_workflow_and_stale_revision(project_root):
    first, second = ToolService(project_root), ToolService(project_root)
    baseline = first.call("open_project")["head_revision"]
    selection = first.call("set_selection", {"expected_revision": baseline, "cell_ids": ["m1"]})
    assert second.call("get_selection")["selection"]["cell_ids"] == ["m1"]
    evidence = second.call("inspect_selection", {"genes": ["CD3D", "ABSENT"]})
    assert evidence["expression"]["CD3D"]["measured"]
    assert not evidence["expression"]["ABSENT"]["measured"]
    proposal = second.call("propose_revision", {"expected_revision": baseline, "selection_id": selection["selection_id"],
        "changes": {"label": "T cell"}, "rationale": "Review marker evidence"})
    assert first.call("open_project")["head_revision"] == baseline
    with pytest.raises(SpatialError):
        first.call("apply_revision", {"proposal_id": proposal["proposal_id"], "expected_revision": baseline,
            "reviewer": "Researcher", "confirmation": False})
    applied = first.call("apply_revision", {"proposal_id": proposal["proposal_id"], "expected_revision": baseline,
        "reviewer": "Researcher", "confirmation": True})
    assert second.call("get_selection")["selection"]["stale"]
    with pytest.raises(SpatialError):
        second.call("set_selection", {"expected_revision": baseline, "cell_ids": ["t1"]})
    result = first.call("run_comparison", {"base_revision": baseline, "target_revision": applied["revision_id"],
        "selection_id": selection["selection_id"]})
    assert result["comparison"]["status"] == "indeterminate"
    assert first.call("export_review_bundle", {"run_id": result["run_id"]})["export_path"]


def test_view_does_not_silently_sample(project_root, monkeypatch):
    service = ToolService(project_root)
    monkeypatch.setattr(service.project, "cells", lambda revision_id: [
        {"cell_id": str(i), "x": float(i), "y": float(i), "label": "T", "included": True, "region": ""}
        for i in range(MAX_VIEW_CELLS + 1)])
    view = service.call("open_project")["view"]
    assert view["cells"] == [] and view["complete"] is False
    assert view["matching_cells"] == MAX_VIEW_CELLS + 1
    smaller = service.call("open_project", {"bounds": [0.0, 0.0, 5.0, 5.0]})["view"]
    assert smaller["complete"] and len(smaller["cells"]) == 6
    with pytest.raises(SpatialError):
        service.call("open_project", {"bounds": [0.0, 0.0, 0.0, 0.0]})


@pytest.fixture
def client(project_root):
    with TestClient(create_app(project_root), base_url="http://127.0.0.1:8765") as value:
        yield value


def browser_headers(client):
    html = client.get("/").text
    token = re.search(r'name="spatial-csrf" content="([^"]+)"', html).group(1)
    return {"Origin": "http://127.0.0.1:8765", "X-Spatial-CSRF": token}


def test_local_page_csp_and_tool_api(client):
    response = client.get("/")
    assert response.status_code == 200
    script = re.search(r"<script>([\s\S]*)</script>", response.text).group(1)
    digest = base64.b64encode(hashlib.sha256(script.encode()).digest()).decode()
    assert f"'sha256-{digest}'" in response.headers["content-security-policy"]
    assert "frame-ancestors 'none'" in response.headers["content-security-policy"]
    assert response.headers["cache-control"] == "no-store"
    result = client.post("/api/tool/open_project", json={}, headers=browser_headers(client))
    assert result.status_code == 200 and result.json()["ok"]
    assert result.json()["result"]["metadata"]["source_kind"] == "synthetic"
    invalid = client.post("/api/tool/open_project", json={"unknown": 1}, headers=browser_headers(client))
    assert invalid.status_code == 400 and not invalid.json()["ok"]


@pytest.mark.parametrize("headers", [
    {}, {"Origin": "https://evil.invalid"}, {"Origin": "null"},
    {"Origin": "http://localhost:8765"}, {"Host": "evil.invalid"},
])
def test_rebinding_and_csrf_are_rejected(client, headers):
    response = client.post("/api/tool/open_project", json={}, headers=headers)
    assert response.status_code == 403


def test_browser_api_rejects_forms_and_oversized_bodies(client):
    headers = browser_headers(client)
    assert client.post("/api/tool/open_project", data="{}", headers=headers).status_code == 415
    assert client.post("/api/tool/open_project", content="x" * (MAX_BODY_BYTES + 1),
                       headers={**headers, "Content-Type": "application/json"}).status_code == 413


def test_streamable_http_negotiation_and_resource(client):
    headers = {"Accept": "application/json, text/event-stream"}
    initialized = client.post("/mcp", json={"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {
        "protocolVersion": "2025-11-25", "capabilities": {}, "clientInfo": {"name": "integration-test", "version": "1"}
    }}, headers=headers)
    assert initialized.status_code == 200
    assert initialized.json()["result"]["serverInfo"]["name"] == "Spatial Collab"
    response = client.post("/mcp", json={"jsonrpc": "2.0", "id": 2, "method": "tools/call",
        "params": {"name": "open_project", "arguments": {}}}, headers=headers)
    assert response.status_code == 200
    assert response.json()["result"]["structuredContent"]["view"]["complete"]
    resource = client.post("/mcp", json={"jsonrpc": "2.0", "id": 4, "method": "resources/read",
        "params": {"uri": UI_URI}}, headers=headers).json()["result"]["contents"][0]
    assert resource["mimeType"] == "text/html;profile=mcp-app"
    assert resource["_meta"]["ui"]["csp"]["connectDomains"] == []
    strict = client.post("/mcp", json={"jsonrpc": "2.0", "id": 5, "method": "tools/call", "params": {
        "name": "apply_revision", "arguments": {"proposal_id": "none", "expected_revision": "none",
            "reviewer": "Researcher", "confirmation": "true"}}}, headers=headers)
    assert strict.json()["result"]["isError"]
    assert "boolean" in strict.json()["result"]["structuredContent"]["error"]
    forbidden = client.post("/mcp", json={"jsonrpc": "2.0", "id": 3, "method": "tools/list"},
        headers={**headers, "Origin": "https://evil.invalid"})
    assert forbidden.status_code == 403
