"""Install a dedicated runtime without changing shell execution policies or system Python."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import venv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--home", default=str(Path(os.environ.get("LOCALAPPDATA", Path.home())) / "SpatialCollab"))
    parser.add_argument("--project")
    parser.add_argument("--create-empty", action="store_true")
    parser.add_argument("--with-import", action="store_true")
    parser.add_argument("--with-gateway", action="store_true")
    parser.add_argument("--wheel", help="Install a verified wheel instead of building this source tree")
    args = parser.parse_args()
    if sys.version_info < (3, 11):
        parser.error("Python 3.11+ is required for initial setup.")
    if bool(args.project) == args.create_empty:
        parser.error("Choose exactly one of --project and --create-empty.")
    state = Path(args.home).expanduser().resolve()
    runtime = state / "runtime"
    executable = runtime / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    if not executable.is_file():
        venv.EnvBuilder(with_pip=True).create(runtime)
    package = str(Path(args.wheel).resolve(strict=True)) if args.wheel else str(Path(__file__).resolve().parents[1])
    extras = (["import"] if args.with_import else []) + (["gateway"] if args.with_gateway else [])
    if extras:
        package += "[" + ",".join(extras) + "]"
    subprocess.run([str(executable), "-I", "-m", "pip", "install", package], check=True)
    env = dict(os.environ, SPATIAL_COLLAB_HOME=str(state))
    configure = [str(executable), "-I", "-m", "spatial_collab.private_runtime", "configure"]
    configure += ["--project", args.project] if args.project else ["--create-empty"]
    subprocess.run(configure, env=env, check=True)
    check = subprocess.run([str(executable), "-I", "-m", "spatial_collab.private_runtime", "doctor"],
                           env=env, check=True, capture_output=True, text=True)
    report = json.loads(check.stdout)
    if not report["local_ready"]:
        raise SystemExit("Private runtime diagnostic failed; inspect runtime.json and dependencies.")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
