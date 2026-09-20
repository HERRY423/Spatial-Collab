"""Portable plugin launcher; project path is explicitly inherited from host environment."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from spatial_collab.__main__ import main

if __name__ == "__main__":
    raise SystemExit(main(["mcp", *sys.argv[1:]]))

