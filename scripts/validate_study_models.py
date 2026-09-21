"""Developer-generated repeated-section controls, with full inputs for replay."""
import json
from pathlib import Path
import runpy

from spatial_collab.demo import create_demo
from spatial_collab.study_inference import infer
from spatial_collab.replay import verify_bundle


def main():
    root = Path("output/scalability/study-tmp")
    root.mkdir(parents=True, exist_ok=False)
    project = create_demo(root / "project")
    controls = runpy.run_path("tests/test_study_inference.py")
    receipts = []
    for method, paired, batches in (("welch", False, False), ("paired_t", True, False), ("mixedlm", True, True)):
        design, rows = controls["design_data"](project, paired=paired, batches=batches, n=24)
        plan = controls["plan"](method)
        result = infer(project, design["object_id"], rows, plan)
        assert result["results"][0]["status"] == "tested", result
        for name, value in (("design", {"study_id": method+"_synthetic_control", "sections": design["sections"]}), ("outcomes", rows), ("plan", plan), ("result", result)):
            (root / f"{method}-{name}.json").write_text(json.dumps(value, indent=2), encoding="utf-8")
        receipts.append({"method": method, "result_id": result["object_id"], "results": result["results"], "scope": "synthetic known generating model, not a biological cohort"})
    bundle = project.export_bundle(compact=True)
    verified = verify_bundle(bundle["export_path"])
    assert all(r["recompute_matches"] for r in verified["study_verification"])
    Path("output/scalability/study-validation.json").write_text(json.dumps({"controls": receipts, "replay": verified, "project": str(project.root)}, indent=2))
    print(json.dumps({"methods": len(receipts), "all_tested": True, "offline_replays_matched": len(verified["study_verification"])}, indent=2))


if __name__ == "__main__":
    main()
