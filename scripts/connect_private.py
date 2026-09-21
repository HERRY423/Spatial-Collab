"""Prepare an official Secure MCP Tunnel after account permissions are granted."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tunnel-id", required=True)
    parser.add_argument("--profile", default="spatial-collab-private")
    parser.add_argument("--home", default=str(Path(os.environ.get("LOCALAPPDATA", Path.home())) / "SpatialCollab"))
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    if not re.fullmatch(r"tunnel_[A-Za-z0-9]+", args.tunnel_id):
        parser.error("Provide a real tunnel_id from OpenAI Platform.")
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.profile):
        parser.error("Invalid profile name.")
    client = shutil.which("tunnel-client")
    if not client:
        parser.error("Install the official tunnel-client after obtaining account permissions. No tunnel was created.")
    if not os.environ.get("CONTROL_PLANE_API_KEY"):
        parser.error("Set CONTROL_PLANE_API_KEY securely in this process; never place it in chat or source files.")
    env = dict(os.environ, SPATIAL_COLLAB_HOME=str(Path(args.home).resolve()))
    python = Path(args.home).resolve() / "runtime" / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    diagnostic = subprocess.run([str(python), "-I", "-m", "spatial_collab.private_runtime", "doctor"],
                                env=env, check=True, capture_output=True, text=True)
    if not json.loads(diagnostic.stdout)["local_ready"]:
        parser.error("Run setup_private.py to initialize a valid local runtime and project first.")
    # Use the installed wheel, with isolated imports. The client's child inherits
    # SPATIAL_COLLAB_HOME; no key is written into its command or this script.
    command = subprocess.list2cmdline([str(python), "-I", "-m", "spatial_collab.private_runtime", "mcp"])
    subprocess.run([client, "init", "--sample", "sample_mcp_stdio_local", "--profile", args.profile,
                    "--tunnel-id", args.tunnel_id, "--mcp-command", command], env=env, check=True)
    subprocess.run([client, "doctor", "--profile", args.profile, "--explain"], env=env, check=True)
    if args.run:
        subprocess.run([client, "run", "--profile", args.profile], env=env, check=True)
    else:
        print(f"Profile prepared. Start with: tunnel-client run --profile {args.profile}")


if __name__ == "__main__":
    main()
