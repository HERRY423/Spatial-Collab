"""Create a deterministic plugin ZIP containing source and host adapters, no project data."""
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    version = json.loads((ROOT / ".codex-plugin/plugin.json").read_text(encoding="utf-8"))["version"]
    paths = []
    for directory in ("src/spatial_collab", "scripts", "skills", "docs", "tests", "requirements", ".codex-plugin", ".claude-plugin"):
        paths.extend(p for p in (ROOT / directory).rglob("*") if p.is_file() and "__pycache__" not in p.parts)
    paths.extend(ROOT / name for name in ("README.md", "pyproject.toml", "MANIFEST.in", ".mcp.json", ".gitignore",
                                         "IMPLEMENTATION_CONTRACT.md", "V02_CONTRACT.md"))
    paths = sorted(set(paths))
    # Freeze the files actually delivered; never reuse a previous release's
    # source manifest or accidentally include synthetic projects/test installs.
    source_paths = [p for p in paths if p.relative_to(ROOT).parts[0] != "docs"]
    source_manifest = {"version": version, "scope": "source and contracts; excludes runtime data, docs, generated wheels",
                       "algorithm": "sha256", "files": {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                                                          for p in source_paths}}
    (ROOT / "docs" / "source-manifest.json").write_text(json.dumps(source_manifest, indent=2) + "\n", encoding="utf-8")
    destination = ROOT / "dist" / f"spatial-collab-{version}.zip"
    destination.parent.mkdir(exist_ok=True)
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in paths:
            info = zipfile.ZipInfo("spatial-collab/" + path.relative_to(ROOT).as_posix(), (2026, 9, 20, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, path.read_bytes())
    sha = hashlib.sha256(destination.read_bytes()).hexdigest()
    artifacts = sorted(p for p in destination.parent.iterdir() if p.is_file() and p.suffix in {".zip", ".whl", ".gz"})
    (destination.parent / "SHA256SUMS.txt").write_text("".join(hashlib.sha256(p.read_bytes()).hexdigest() + "  " + p.name + "\n"
                                                               for p in artifacts), encoding="utf-8")
    print(json.dumps({"path": str(destination), "sha256": sha, "files": len(paths), "contains_project_data": False}, indent=2))


if __name__ == "__main__":
    main()
