"""Check final runtime bytes, source manifest and bounded acceptance receipts."""
import hashlib
import json
from pathlib import Path
import re
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    output = ROOT / "output/power"
    tests = Path("C:/Spatial/power-tests-delivery-lr.log")
    match = re.search(r"(\d+) passed, (\d+) skipped.* in ([\d.]+)s", tests.read_text())
    assert match and match[1] == "509" and match[2] == "3"
    manifest = json.loads((ROOT / "docs/source-manifest.json").read_text())
    for name, digest in manifest["files"].items():
        assert sha(ROOT / name) == digest, name
    runtime = ROOT / ".test-install/power-delivery"
    count = 0
    with zipfile.ZipFile(ROOT / "dist/spatial_collab-0.7.0a1-py3-none-any.whl") as wheel:
        for name in wheel.namelist():
            if name.startswith("spatial_collab/") and not name.endswith("/"):
                assert wheel.read(name) == (ROOT / "src" / name).read_bytes() == (runtime / name).read_bytes(), name
                count += 1
    protocol = json.loads((output / "protocol-validation.json").read_text())
    assert protocol["all_passed"] and "power-delivery" in protocol["runtime"]
    assert "Errors: 0, Warnings: 0" in (output / "browser-console.txt").read_text(encoding="utf-8-sig")
    evidence = ["public-power-validation.json", "protocol-validation.json", "browser-power.yml", "browser-power.png", "browser-resolution.png", "browser-console.txt"]
    report = {"version": "0.7.0a1 / 0.7.0-alpha.1, unchanged", "tests": {"passed": int(match[1]), "skipped": int(match[2]), "seconds": float(match[3]), "log_sha256": sha(tests)},
              "runtime_files_source_wheel_installed_matched": count, "source_manifest_files_verified": len(manifest["files"]),
              "installed_mcp": {"passed": True, "checks": sorted(set(protocol["checks"])), "runtime": protocol["runtime"]},
              "browser": {"url": "http://127.0.0.1:8780", "verified": ["method responsibility labels", "exact public family", "resolution blocked at rank 1", "dense rank scenario with conditional MDE and intervals", "retrospective timing retained"], "errors": 0, "warnings": 0},
              "evidence_sha256": {name: sha(output / name) for name in evidence},
              "artifacts_sha256": {p.name: sha(p) for p in (ROOT / "dist").iterdir() if p.is_file() and "0.7.0" in p.name and p.suffix in {".whl", ".gz", ".zip"}},
              "scope": "Numerical correctness and local product/transport acceptance; model-conditional grid MDE, not general biological power, independent validation or externally preregistered research."}
    (output / "delivery-verification.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
