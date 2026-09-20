"""Exercise actual MCP client/server transports, not merely direct Python calls."""
import asyncio
import os
from pathlib import Path
import socket
import subprocess
import sys

import httpx
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from mcp.client.streamable_http import streamable_http_client

from spatial_collab.demo import create_demo

ROOT = Path(__file__).resolve().parents[1]


def process_env(project):
    return {**os.environ, "PYTHONPATH": str(ROOT / "src"), "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONIOENCODING": "utf-8", "SPATIAL_COLLAB_PROJECT": str(project)}


async def exercise(session):
    await session.initialize()
    catalog = await session.list_tools()
    assert {"get_context", "query_observations", "get_proposal", "get_capabilities", "run_region_comparison"} <= {tool.name for tool in catalog.tools}
    assert any(tool.name == "open_project" for tool in catalog.tools)
    opened = await session.call_tool("open_project", {})
    assert not opened.isError
    summary = opened.structuredContent
    assert summary["cell_count"] == 480
    base = summary["head_revision"]
    resource = await session.read_resource("ui://spatial-collab/workbench.html")
    assert "mcp-app" in resource.contents[0].mimeType
    assert "Spatial Collab" in resource.contents[0].text
    picked = await session.call_tool("set_selection", {"expected_revision": base,
        "cell_ids": ["demo-08-10", "demo-08-11"], "name": "protocol-test"})
    assert not picked.isError
    inspected = await session.call_tool("inspect_selection", {"genes": ["CD3D", "UNMEASURED"]})
    assert inspected.structuredContent["expression"]["UNMEASURED"]["measured"] is False
    proposal = await session.call_tool("propose_revision", {"expected_revision": base,
        "selection_id": picked.structuredContent["selection_id"], "changes": {"label": "Uncertain"},
        "rationale": "Synthetic protocol test, not biological adjudication."})
    proposal_id = proposal.structuredContent["proposal_id"]
    context = await session.call_tool("get_context", {})
    assert context.structuredContent["latest_proposal_id"] == proposal_id
    handoff = await session.call_tool("get_proposal", {"proposal_id": proposal_id})
    assert handoff.structuredContent["deltas"] == proposal.structuredContent["deltas"]
    query = await session.call_tool("query_observations", {"revision_id": base, "limit": 5})
    assert len(query.structuredContent["observations"]) == 5 and query.structuredContent["next_cursor"]
    refused = await session.call_tool("apply_revision", {"proposal_id": proposal_id,
        "expected_revision": base, "reviewer": "synthetic-test", "confirmation": False})
    assert refused.isError
    applied = await session.call_tool("apply_revision", {"proposal_id": proposal_id,
        "expected_revision": base, "reviewer": "synthetic-test", "confirmation": True})
    assert not applied.isError
    target = applied.structuredContent["revision_id"]
    run = await session.call_tool("run_comparison", {"base_revision": base, "target_revision": target})
    assert not run.isError
    assert run.structuredContent["scientific_authorization"] == "NOT_ESTABLISHED"
    exported = await session.call_tool("export_review_bundle", {"run_id": run.structuredContent["run_id"]})
    assert not exported.isError
    assert Path(exported.structuredContent["manifest_path"]).is_file()
    regional = await session.call_tool("run_region_comparison", {"base_revision": base, "target_revision": target,
        "selection_id": picked.structuredContent["selection_id"], "genes": ["CD3D", "UNMEASURED"]})
    assert not regional.isError
    assert regional.structuredContent["before"]["foreground"]["markers"]["UNMEASURED"]["status"] == "unmeasured"
    stale_apply = await session.call_tool("apply_revision", {"proposal_id": proposal_id,
        "expected_revision": base, "reviewer": "synthetic-test", "confirmation": True})
    assert stale_apply.isError


def test_real_stdio_workflow(tmp_path):
    project = tmp_path / "stdio"
    create_demo(project)

    async def run():
        params = StdioServerParameters(command=sys.executable, args=[str(ROOT / "scripts/run_server.py")],
                                       env=process_env(project))
        async with stdio_client(params) as (read, write):
            async with ClientSession(read, write) as session:
                await exercise(session)
    asyncio.run(run())


def test_real_streamable_http_workflow(tmp_path):
    project = tmp_path / "http"
    create_demo(project)
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    with (tmp_path / "server.log").open("w", encoding="utf-8") as log:
        child = subprocess.Popen([sys.executable, "-m", "spatial_collab", "serve", "--project", str(project),
                                  "--port", str(port)], env=process_env(project), cwd=ROOT,
                                 stdout=log, stderr=log)
        async def run():
            async with httpx.AsyncClient() as client:
                for _ in range(100):
                    if child.poll() is not None:
                        raise AssertionError("HTTP test server stopped before becoming ready")
                    try:
                        response = await client.get(f"http://127.0.0.1:{port}/", timeout=.5)
                        if response.status_code == 200:
                            break
                    except httpx.TransportError:
                        pass
                    await asyncio.sleep(.1)
                else:
                    raise AssertionError("HTTP test server did not become ready in 10 seconds")
            async with streamable_http_client(f"http://127.0.0.1:{port}/mcp") as (read, write, _):
                async with ClientSession(read, write) as session:
                    await exercise(session)
        try:
            asyncio.run(run())
        finally:
            child.terminate()
            child.wait(timeout=10)
