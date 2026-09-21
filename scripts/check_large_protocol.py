"""Actual MCP transport acceptance of disk atlas, TIFF tile and study model."""
import asyncio
from datetime import timedelta
import json
import os
from pathlib import Path
import sys

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client


async def main():
    scale = json.loads(Path(sys.argv[1]).read_text())
    env = {**os.environ, "PYTHONIOENCODING": "utf-8", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    params = StdioServerParameters(command=sys.executable, args=["-m", "spatial_collab", "mcp", "--project", scale["project"]], env=env)
    checks = []
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write, read_timeout_seconds=timedelta(seconds=45)) as session:
            await session.initialize()
            async def call(name, args):
                r = await session.call_tool(name, args)
                assert not r.isError, r
                checks.append(name)
                return r.structuredContent
            await call("list_atlases", {})
            census = await call("get_atlas_view", {"atlas_id": scale["atlas_id"]})
            assert census["aggregate_total"] == 1000000
            shared = await call("share_atlas_view", {"atlas_id": scale["atlas_id"], "bounds": [-.5, -.5, 1.5, 1.5], "layer": "transcripts", "feature": "G0", "image_id": scale["image_id"]})
            context = await call("get_context", {})
            assert context["atlas_view"] == shared
            molecules = await call("get_atlas_view", {"atlas_id": scale["atlas_id"], "layer": "transcripts", "feature": "G0"})
            assert molecules["total"] > 0
            image = await call("get_pyramid_tile", {"image_id": scale["image_id"], "level": 0, "x": 0, "y": 0})
            assert image["data_url"].startswith("data:image/png;base64,")
            frozen = await call("freeze_atlas_roi", {"atlas_id": scale["atlas_id"], "bounds": [-.5, -.5, 1.5, 1.5]})
            assert frozen["shape"] == [4, 30001]
            folder = Path("output/scalability/study-tmp")
            design = await call("register_study_design", {"spec": json.loads((folder / "mixedlm-design.json").read_text())})
            result = await call("run_study_inference", {"study_id": design["object_id"], "records": json.loads((folder / "mixedlm-outcomes.json").read_text()), "plan": json.loads((folder / "mixedlm-plan.json").read_text())})
            assert result["results"][0]["status"] == "tested"
            await call("list_study_analyses", {})
    report = {"transport": "actual MCP stdio; not host adoption or independent scientific validation", "checks": checks, "all_passed": True, "census": census["aggregate_total"], "frozen_shape": frozen["shape"], "mixedlm_status": result["results"][0]["status"]}
    Path("output/scalability/protocol-validation.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    asyncio.run(main())
