"""Private installation state and connection diagnostics; no credentials in receipts."""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import sys


def home():
    return Path(os.environ.get("SPATIAL_COLLAB_HOME", str(Path(os.environ.get("LOCALAPPDATA", Path.home())) / "SpatialCollab"))).expanduser().resolve()


def config():
    path = home() / "runtime.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}


def configure(project=None, *, create_empty=False):
    from .store import Project, SpatialError
    root = home()
    root.mkdir(parents=True, exist_ok=True)
    if project:
        target = Path(project).expanduser().resolve(strict=True)
        Project(target)
    elif create_empty:
        target = root / "projects" / "private"
        if not (target / "project.sqlite3").exists():
            Project.create_workspace(target, "My private spatial workspace")
        Project(target)
    else:
        raise SpatialError("Choose an existing project or explicitly create an empty private workspace.")
    data = {"schema": "spatial-collab.private-runtime.v1", "python": sys.executable, "project": str(target),
            "data_policy": "Tool results reach the host; raw sources remain local unless explicitly transferred.",
            "chatgpt": {"status": "not_connected", "reason": "Requires developer-mode/tunnel permissions and live host acceptance."}}
    destination = root / "runtime.json"
    temporary = root / "runtime.json.tmp"
    temporary.write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8")
    temporary.replace(destination)
    return data


def project_path():
    value = os.environ.get("SPATIAL_COLLAB_PROJECT") or config().get("project")
    if not value:
        raise ValueError("Run scripts/setup_private.py first to create or choose a private project.")
    return Path(value).expanduser().resolve(strict=True)


def doctor():
    from . import __version__
    versions = {}
    for name in ("numpy", "scipy", "mcp", "uvicorn", "starlette", "anndata", "h5py", "PyJWT", "cryptography"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    try:
        from .store import Project
        summary = Project(project_path()).summary()
        project = {"ready": True, "project_id": summary["project_id"], "name": summary["metadata"]["name"]}
    except (OSError, ValueError) as exc:
        project = {"ready": False, "reason": str(exc)}
    return {"version": __version__, "python": sys.executable, "packages": versions, "project": project,
            "local_ready": project["ready"] and all(versions[n] for n in ("numpy", "scipy", "mcp", "uvicorn", "starlette")),
            "tunnel_client_available": shutil.which("tunnel-client") is not None,
            "tunnel_key_present": bool(os.environ.get("CONTROL_PLANE_API_KEY")),
            "chatgpt_host_verified": False,
            "notice": "This diagnostic does not authenticate to ChatGPT or establish host acceptance."}


def main(argv=None):
    parser = argparse.ArgumentParser(description="Configure or inspect a private Spatial Collab installation")
    parser.add_argument("action", choices=["configure", "doctor", "mcp"])
    parser.add_argument("--project")
    parser.add_argument("--create-empty", action="store_true")
    args = parser.parse_args(argv)
    if args.action == "configure":
        result = configure(args.project, create_empty=args.create_empty)
    elif args.action == "doctor":
        result = doctor()
    else:
        from .server import create_server
        create_server(project_path()).run(transport="stdio")
        return
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
