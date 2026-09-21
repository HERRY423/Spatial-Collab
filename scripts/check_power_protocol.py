"""Installed-runtime MCP round trip for prospective declarations and responsibility."""
import asyncio
from datetime import timedelta
import json
import os
from pathlib import Path
import sys

import numpy as np
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from spatial_collab.store import Project
from spatial_collab.workflow_inputs import register_counts


async def main():
    folder = Path(sys.argv[1]).resolve()
    p = Project.create_workspace(folder, "Power declaration numerical control")
    r = register_counts(p, np.random.default_rng(9).poisson(5, (64, 3))+1,
                        [str(i) for i in range(64)], ["A", "B", "C"], sample_id="synthetic",
                        species="mouse", kind="spatial", coordinates=[[i % 8, i // 8] for i in range(64)],
                        provenance={"synthetic": True})
    params = StdioServerParameters(command=sys.executable, args=["-m", "spatial_collab", "mcp", "--project", str(folder)],
                                  env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    checks = []
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write, read_timeout_seconds=timedelta(seconds=65)) as session:
            await session.initialize()
            async def call(name, args):
                answer = await session.call_tool(name, args)
                assert not answer.isError, answer
                checks.append(name)
                return answer.structuredContent
            catalog = await call("get_analysis_catalog", {})
            spec = {"method": "moran_svg", "input_ids": [r["object_id"]], "parameters": {"permutations": 99}}
            plan = await call("create_analysis_power_plan", {"spec": spec, "assumptions": {
                "model": "gaussian_sar", "simulations": 80, "effect_grid": [0, .4, .8, .95],
                "rationale": "Prospective synthetic transport control, not biological calibration."}})
            assert plan["timing"] == "before_target_run_in_this_workspace"
            reread = await call("get_analysis_power_plan", {"plan_id": plan["object_id"]})
            assert reread["object_sha256"] == plan["object_sha256"]
            job = await call("submit_analysis", {"spec": {**spec, "power_plan_id": plan["object_id"]}})
            for _ in range(100):
                status = await call("get_analysis_job", {"job_id": job["id"]})
                if status["status"] in {"succeeded", "failed", "cancelled"}:
                    break
                await asyncio.sleep(.2)
            assert status["status"] == "succeeded", status
            result = await call("inspect_analysis_result", {"result_id": status["result_id"]})
            assert result["metadata"]["power_evidence"]["plan"]["plan_id"] == plan["object_id"]
            assert result["method_contract"]["tier"] == "internal_numerical"
    import spatial_collab
    report = {"runtime": spatial_collab.__file__, "project": str(folder), "checks": checks,
              "all_passed": True, "tiers": {r["method"]: r["method_contract"]["tier"] for r in catalog["methods"]},
              "power_plan": plan, "result_id": status["result_id"], "environment_sha256": result["execution_receipt"]["environment_sha256"],
              "scope": "Actual MCP stdio and configured child process; synthetic prospective ordering, not external preregistration or host adoption."}
    Path("output/power/protocol-validation.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({"all_passed": True, "runtime": report["runtime"], "checks": sorted(set(checks))}))


if __name__ == "__main__":
    asyncio.run(main())
