"""Create unfilled, counterbalanced independent-user trial materials."""
import argparse
import csv
import itertools
import json
from pathlib import Path

ARMS = ("notebook_viewer", "generic_agent_tools", "agent_spatial_collab")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--participants", type=int, default=6)
    args = parser.parse_args()
    if not 6 <= args.participants <= 60:
        parser.error("Use 6..60 planned participants. Enrollment is not simulated.")
    args.output.mkdir(parents=True, exist_ok=True)
    sequences = list(itertools.permutations(ARMS))
    with (args.output / "assignment.csv").open("x", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["participant_code", "period", "workflow", "task_variant", "enrolled"])
        for i in range(args.participants):
            for period, arm in enumerate(sequences[i % 6]):
                w.writerow([f"P{i+1:02}", period+1, arm, f"T{(period+i//6)%3+1}", ""])
    fields = ["participant_code", "period", "workflow", "task_variant", "dataset_sha256", "model_version", "backend_version",
              "hardware", "active_seconds", "waiting_seconds", "total_seconds", "correctly_completed", "rework_count",
              "omitted_steps", "replay_success", "peak_memory_mb", "compute_seconds", "critical_identity_error",
              "result_bundle", "reviewer_code", "reviewer_independence", "notes"]
    with (args.output / "observations.csv").open("x", newline="", encoding="utf-8") as f:
        csv.writer(f).writerow(fields)
    protocol = {
        "status": "prepared_not_enrolled", "planned_participants": args.participants, "observed_participants": 0,
        "arms": ARMS, "assignment": "All six workflow orders balanced in blocks of six; three disjoint task variants, fixed across arms.",
        "tasks": [
            {"id": "T1", "question": "Compare RNA, protein and joint domains; inspect a high-disagreement region, cite raw measured features and decide retain/uncertain/exclude with a reason."},
            {"id": "T2", "question": "Use a distinct predefined region. Apply a justified inclusion hypothesis in a working copy; compare fixed-model filtering with refitting and explain denominators."},
            {"id": "T3", "question": "Resume an interrupted review on another predefined region, recover exact versions, diagnose a supplied stale-result or foreign-ID case and replay the exported computation."}],
        "before_enrollment": ["Independent scientist fixes three task ROIs and evidence rubrics before seeing trial outputs.",
             "Separate training task and washout; match backend, model and machine settings across arms.",
             "Obtain participant consent; use anonymous participant codes. Do not record patient information.",
             "Recruit 6-10 researchers independently; developer-operated runs cannot fill human observations."],
        "rubric": {"identity_and_version": 1, "measured_vs_unmeasured": 1, "raw_evidence_cited": 1,
                   "fixed_vs_refit_distinguished": 1, "denominators_explained": 1, "replay_success": 1,
                   "no_unsupported_biological_claim": 1},
        "critical_failure": "Silent object mismatch, version mixing, invented measurements or unapproved annotation mutation.",
        "primary_endpoints": ["Within-participant active time ratio at non-decreased independently graded correctness.",
             "Critical failures and correctly completed tasks."],
        "secondary_endpoints": ["Waiting and total time separately", "Rework", "Omissions", "Replay success", "Resource consumption"],
        "target_not_claim": "Approximately 30% less active time without lower correctness; exploratory pilot, not proven benefit.",
        "scientific_authorization": "NOT_ESTABLISHED"}
    (args.output / "protocol.json").write_text(json.dumps(protocol, ensure_ascii=False, indent=2), encoding="utf-8")
    guide = Path(__file__).resolve().parents[1] / "docs" / "V06_PILOT_GUIDE.md"
    (args.output / "README.md").write_text(guide.read_text(encoding="utf-8"), encoding="utf-8")
    (args.output / "frozen_tasks.template.json").write_text(json.dumps({
        "status": "awaiting_independent_scientist", "reviewer_code": None,
        "tasks": [{"id": f"T{i}", "source_sha256": None, "revision_id": None,
                   "observation_ids": [], "grading_evidence": [], "frozen_at": None} for i in range(1, 4)]
    }, indent=2), encoding="utf-8")
    print(json.dumps({"pilot_dir": str(args.output.resolve()), "observed_participants": 0}))


if __name__ == "__main__":
    main()
