"""Verify exact private release bytes and real installed MCP transports.

Run with the dedicated runtime's Python -I. No ChatGPT acceptance is inferred.
"""
import argparse
import asyncio
import base64
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import sys
import tempfile
import zipfile

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

import spatial_collab
from spatial_collab.demo import create_demo
from spatial_collab.store import Project

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


async def protocol(params, expected_cells):
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            tools = (await session.list_tools()).tools
            async def call(name, **args):
                result = await session.call_tool(name, args)
                assert not result.isError, (name, result)
                return result.structuredContent
            opened = await call("open_project")
            assert opened["cell_count"] == expected_cells
            assert {"begin_upload", "import_uploaded_project", "read_download"} <= {t.name for t in tools}
            resource = (await session.read_resource("ui://spatial-collab/workbench.html")).contents[0]
            assert "uploadAndImport" in resource.text and "ui/initialize" in resource.text
            download = await call("prepare_download", project_id=opened["project_id"])
            data = bytearray()
            while len(data) < download["size"]:
                chunk = await call("read_download", file_id=download["file_id"], offset=len(data), project_id=opened["project_id"])
                data.extend(base64.b64decode(chunk["data_base64"]))
            assert hashlib.sha256(data).hexdigest() == download["sha256"]
            await call("discard_transfer", file_id=download["file_id"])
            return {"tool_count": len(tools), "cell_count": expected_cells, "ui_resource": True,
                    "download_sha256_verified": True, "real_chatgpt": False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--home", required=True)
    parser.add_argument("--plugin-root", default=str(ROOT))
    parser.add_argument("--test-log", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    installed = Path(spatial_collab.__file__).parent
    assert not installed.is_relative_to(ROOT / "src"), "Run with the dedicated Python -I, not checkout imports."
    version = importlib.metadata.version("spatial-collab")
    manifest = json.loads((ROOT / "docs/source-manifest.json").read_text())
    for name, value in manifest["files"].items():
        assert digest(ROOT / name) == value, name
    wheel = ROOT / "dist" / f"spatial_collab-{version}-py3-none-any.whl"
    archive = ROOT / "dist" / f"spatial-collab-{manifest['version']}.zip"
    count = 0
    with zipfile.ZipFile(wheel) as bundle:
        for name in bundle.namelist():
            if name.startswith("spatial_collab/") and not name.endswith("/"):
                relative = name.removeprefix("spatial_collab/")
                assert bundle.read(name) == (ROOT / "src" / name).read_bytes() == (installed / relative).read_bytes(), name
                count += 1
    zip_count = 0
    with zipfile.ZipFile(archive) as bundle:
        for name in bundle.namelist():
            relative = name.removeprefix("spatial-collab/")
            assert name != relative and bundle.read(name) == (ROOT / relative).read_bytes(), name
            assert bundle.read(name) == (Path(args.plugin_root) / relative).read_bytes(), name
            zip_count += 1
    log = Path(args.test_log)
    match = re.search(r"(\d+) passed(?:, (\d+) skipped)?.* in ([\d.]+)s", log.read_text())
    assert match and not re.search(r"\d+ failed", log.read_text()), "A completed, passing full suite is required."
    settings = json.loads((Path(args.home) / "runtime.json").read_text())
    assert Project(settings["project"]).summary()["cell_count"] == 0, "The configured private workspace should remain empty."
    # Both modes use separate synthetic projects; the user's workspace is read only.
    env = dict(os.environ, SPATIAL_COLLAB_HOME=args.home, PYTHONIOENCODING="utf-8")
    results = {}
    with tempfile.TemporaryDirectory(prefix="private-protocol-", dir=Path(args.home)) as temporary:
        for mode in ("installed_wheel", "plugin_launcher"):
            project = Path(temporary) / mode
            create_demo(project)
            child_env = dict(env, SPATIAL_COLLAB_PROJECT=str(project))
            if mode == "installed_wheel":
                params = StdioServerParameters(command=sys.executable, args=["-I", "-m", "spatial_collab.private_runtime", "mcp"], env=child_env)
            else:
                params = StdioServerParameters(command="cmd", args=["/d", "/c", str(Path(args.plugin_root) / "scripts/launch_private.cmd")], env=child_env)
            results[mode] = asyncio.run(protocol(params, 480))
    report = {"version": version, "source_manifest_files_verified": len(manifest["files"]),
              "source_wheel_installed_files_matched": count, "source_zip_plugin_files_matched": zip_count,
              "runtime": sys.executable, "plugin_root": str(Path(args.plugin_root).resolve()),
              "private_workspace_empty": True, "protocol": results,
              "tests": {"passed": int(match[1]), "skipped": int(match[2] or 0), "seconds": float(match[3]), "log_sha256": digest(log)},
              "artifacts": {p.name: digest(p) for p in (ROOT / "dist").iterdir() if p.suffix in {".whl", ".zip", ".gz"} and (version in p.name or manifest["version"] in p.name)},
              "chatgpt": {"status": "BLOCKED_ACCOUNT_PERMISSIONS", "live_iframe_accepted": False,
                          "reason": "User reports no developer mode or tunnel permissions."},
              "public_service": "RESOURCE_SERVER_TESTED_LOCALLY_NOT_DEPLOYED",
              "evidence_scope": "Local engineering, exact package bytes and real MCP transports; no biological validation or live ChatGPT acceptance."}
    destination = Path(args.output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
