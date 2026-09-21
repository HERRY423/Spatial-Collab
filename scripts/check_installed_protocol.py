"""Bounded installed-wheel stdio smoke with per-request diagnostics."""
import asyncio
from datetime import timedelta
import os
from pathlib import Path
import sys

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from spatial_collab.demo import create_demo


async def main():
    project = Path(sys.argv[1]).resolve()
    create_demo(project)
    params = StdioServerParameters(command=sys.executable,
        args=["-m", "spatial_collab", "mcp", "--project", str(project)],
        env={**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUNBUFFERED": "1"})
    if os.environ.get("SPATIAL_PROTOCOL_DEBUG") == "1":
        params.args = ["-c", "import faulthandler; faulthandler.dump_traceback_later(10); from spatial_collab.__main__ import main; import sys; main(['mcp','--project',sys.argv[1]])", str(project)]
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write, read_timeout_seconds=timedelta(seconds=20)) as session:
            await session.initialize()
            catalog = await session.list_tools()
            assert {"submit_integration", "compare_integrations", "get_integration_sensitivity", "submit_analysis", "list_analysis_jobs", "inspect_analysis_result"} <= {t.name for t in catalog.tools}
            for name in ("open_project", "get_context", "list_integrations", "get_multimodal_context", "get_analysis_catalog", "list_analysis_inputs", "list_analysis_jobs", "list_analysis_results"):
                print("CALL", name, flush=True)
                r = await session.call_tool(name, {})
                assert not r.isError, r
                print("PASS", name, flush=True)
            context = (await session.call_tool("get_context", {})).structuredContent
            frozen = await session.call_tool("prepare_analysis_input", {"revision_id": context["head_revision"], "species": "human"})
            assert not frozen.isError, frozen
            job = await session.call_tool("submit_analysis", {"spec": {"method": "moran_svg", "input_ids": [frozen.structuredContent["input_id"]], "parameters": {"permutations": 19}, "seed": 19}})
            assert not job.isError, job
            job_id = job.structuredContent["id"]
            for _ in range(60):
                status = (await session.call_tool("get_analysis_job", {"job_id": job_id})).structuredContent
                if status["status"] in {"succeeded", "failed", "cancelled"}:
                    break
                await asyncio.sleep(0.5)
            assert status["status"] == "succeeded", status
            inspected = await session.call_tool("inspect_analysis_result", {"result_id": status["result_id"], "limit": 3})
            assert not inspected.isError and inspected.structuredContent["rows"], inspected
            print("PASS analysis_freeze_submit_poll_inspect", flush=True)
    print("STDIO_CLOSED", flush=True)


if __name__ == "__main__":
    asyncio.run(main())
