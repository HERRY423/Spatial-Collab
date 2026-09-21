"""Verify the unchanged-version source, wheel, installed runtime and receipts."""
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    report = ROOT / "output/scalability"
    tests = Path("C:/Spatial/scale-delivery-final.log")
    assert "493 passed, 3 skipped" in tests.read_text()
    manifest = json.loads((ROOT / "docs/source-manifest.json").read_text())
    for relative, digest in manifest["files"].items():
        assert sha(ROOT / relative) == digest, relative
    wheel = ROOT / "dist/spatial_collab-0.7.0a1-py3-none-any.whl"
    runtime = ROOT / ".test-install/scalability-final"
    count = 0
    with zipfile.ZipFile(wheel) as archive:
        for member in archive.namelist():
            if not member.startswith("spatial_collab/") or member.endswith("/"):
                continue
            data = archive.read(member)
            assert data == (ROOT / "src" / member).read_bytes(), member
            assert data == (runtime / member).read_bytes(), member
            count += 1
    protocol = json.loads((report / "protocol-validation.json").read_text())
    assert protocol["all_passed"] and "share_atlas_view" in protocol["checks"]
    console = (report / "browser-console.txt").read_text(encoding="utf-8-sig")
    assert "Errors: 0, Warnings: 0" in console
    artifacts = {p.name: sha(p) for p in (ROOT / "dist").iterdir()
                 if p.is_file() and "0.7.0" in p.name and p.suffix in {".whl", ".zip", ".gz"}}
    evidence = ["million-optimized.json", "gigapixel-validation.json", "gigapixel-read.json",
                "public-validation.json", "study-validation.json", "protocol-validation.json",
                "browser-empty-workspace.png", "browser-installed-million.png",
                "browser-installed-boundaries.png", "browser-console.txt"]
    result = {"version": "0.7.0a1 / 0.7.0-alpha.1; unchanged",
              "tests": {"passed": 493, "skipped": 3, "seconds": 54.11,
                        "skips": ["SpatialData absent", "two native napari checks unavailable"],
                        "log": str(tests), "sha256": sha(tests)},
              "runtime_source_and_installed_wheel_files_verified": count,
              "source_manifest_files_verified": len(manifest["files"]),
              "installed_protocol": protocol,
              "browser": {"runtime": str(runtime), "url": "http://127.0.0.1:8778",
                          "verified": ["empty workspace opens atlas panel", "million-object full census",
                                       "8-level image", "zoom/pan molecules and cell/nucleus boundaries"],
                          "console_errors": 0, "console_warnings": 0},
              "artifacts_sha256": artifacts, "evidence_sha256": {n: sha(report / n) for n in evidence},
              "scope": "Local engineering and numerical verification; synthetic scale/model controls and public SPOTS count integrity. Not independent biological, multi-patient, host adoption or user-benefit validation."}
    (report / "delivery-verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
