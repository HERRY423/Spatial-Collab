"""Acceptance outputs must reflect actual product results and retained failures."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from spatial_collab import benchmark
from spatial_collab.store import SpatialError


def test_end_to_end_tasks_have_observed_artifacts_and_explicit_evidence_ceiling(tmp_path):
    report = benchmark.run_benchmark(tmp_path)
    assert report["task_count"] == 8
    assert report["all_passed"], report["results"]
    assert report["passed_count"] == 8
    assert report["failed_count"] == 0
    assert report["external_validation"] == "NOT_ESTABLISHED"
    assert report["model_evaluation"] == "NOT_PERFORMED"
    assert report["human_efficiency_evaluation"] == "NOT_PERFORMED"
    assert report["spatialbench_evaluation"] == "NOT_PERFORMED"
    assert report["execution_mode"] == "scripted_reference_driver"
    folder = Path(report["output_dir"])
    manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
    for name, receipt in manifest["files"].items():
        content = (folder / name).read_bytes()
        assert receipt["bytes"] == len(content)
        assert receipt["sha256"] == hashlib.sha256(content).hexdigest()
    assert all(result["trace"] for result in report["results"])
    assert all(Path(result["result_path"]).is_file() for result in report["results"])
    selected = next(result for result in report["results"] if result["task_id"] == "annotation_sensitivity")
    assert selected["answer"]["before_effect"] == 0.5
    assert selected["answer"]["after_effect"] == -0.25
    assert selected["answer"]["region_only_status"] == "stable"


def test_reruns_never_overwrite_existing_projects_or_results(tmp_path):
    first = benchmark.run_benchmark(tmp_path)
    first_report = Path(first["output_dir"]) / "report.json"
    original = first_report.read_bytes()
    second = benchmark.run_benchmark(tmp_path)
    assert first["output_dir"] != second["output_dir"]
    assert first_report.read_bytes() == original
    assert first["all_passed"] and second["all_passed"]


def test_numeric_regression_cannot_be_hidden_by_correct_status_labels(tmp_path, monkeypatch):
    original = benchmark.analysis.compare

    def corrupted(*args, **kwargs):
        result = original(*args, **kwargs)
        if result["before"]["excess_over_abundance"] is not None:
            result["before"]["excess_over_abundance"] += 0.25
        return result

    monkeypatch.setattr(benchmark.analysis, "compare", corrupted)
    report = benchmark.run_benchmark(tmp_path)
    assert not report["all_passed"]
    results = {result["task_id"]: result for result in report["results"]}
    assert not results["annotation_sensitivity"]["passed"]
    assert not results["neighbor_context"]["passed"]
    assert results["spatial_object_identity"]["passed"]
    assert report["failed_count"] >= 2


def test_failed_task_is_retained_and_does_not_hide_remaining_tasks(tmp_path, monkeypatch):
    original = benchmark._execute

    def failed(task_id, folder, trace):
        if task_id == "panel_evidence":
            trace.append({"operation": "synthetic_successful_step_before_failure"})
            raise RuntimeError("synthetic injected failure")
        return original(task_id, folder, trace)

    monkeypatch.setattr(benchmark, "_execute", failed)
    report = benchmark.run_benchmark(tmp_path)
    assert report["task_count"] == 8
    assert report["passed_count"] == 7
    assert report["failed_count"] == 1
    assert report["results"][0]["error"]["type"] == "RuntimeError"
    assert report["results"][0]["trace"][0]["operation"] == "synthetic_successful_step_before_failure"
    assert Path(report["results"][0]["result_path"]).is_file()


@pytest.mark.parametrize("answer", [{}, {"observation_unit": "cell"},
                                     {"observation_unit": "spot", "observation_count": True},
                                     {"observation_unit": "spot", "observation_count": float("nan")}])
def test_grader_does_not_reward_missing_fields_wrong_units_boolean_counts_or_nan(answer):
    assert not benchmark.grade_answer("spot_observation_semantics", answer)["passed"]


def test_grader_preserves_missing_versus_null():
    answer = {"status": "indeterminate", "after_classification": "inconclusive",
              "reason": "target_label_absent_from_neighbor_universe"}
    assert not benchmark.grade_answer("missing_population", answer)["passed"]
    answer["effect_delta"] = None
    assert benchmark.grade_answer("missing_population", answer)["passed"]


def test_grader_rejects_duplicate_identity_and_unknown_tasks():
    answer = {"selected_ids": ["near-other", "near-target", "source", "source"],
              "units": "micrometer", "coordinate_system": "synthetic_xy",
              "revision_matches": True, "source_matches": True}
    assert not benchmark.grade_answer("spatial_object_identity", answer)["passed"]
    with pytest.raises(SpatialError, match="Unknown synthetic task"):
        benchmark.grade_answer("external_scientific_truth", {})
    with pytest.raises(SpatialError, match="JSON object"):
        benchmark.grade_answer("missing_population", [])
