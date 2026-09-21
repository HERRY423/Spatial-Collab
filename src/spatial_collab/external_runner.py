"""Explicit export/import contract for externally executed paired integration.

Exports frozen input contracts and ingests validated external integration
results with provenance and object identity verification. This bridge does
not execute external code; workflow_jobs provides the actual method runner.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
from typing import Any

from .integration import input_identity, register_result, SCHEMA
from .proteomics import get_assay
from .store import Project, SpatialError, _json


def export_input_contract(
    project: Project,
    destination: str | Path,
    *,
    revision_id: str | None = None,
    assay_id: str | None = None,
    observation_ids: list[str] | None = None,
    rna_features: list[str] | None = None,
    protein_features: list[str] | None = None,
) -> dict[str, Any]:
    """Export standardized inputs for external integration pipelines."""
    dest = Path(destination).resolve()
    if not assay_id:
        raise SpatialError(
            "The legacy integration bridge requires a paired protein assay. Use analysis snapshot/import for RNA-only scientific workflows."
        )
    if dest.exists() and any(dest.iterdir()):
        raise SpatialError("Use a new empty export directory; frozen input contracts are not overwritten.")

    summary = project.summary()
    rev_id = revision_id or summary["head_revision"]
    cells = [c for c in project.cells(rev_id) if c["included"]]
    if observation_ids is not None:
        if (
            not observation_ids
            or len(observation_ids) != len(set(observation_ids))
            or not set(observation_ids) <= {c["cell_id"] for c in cells}
        ):
            raise SpatialError("Explicit export IDs must be unique included observations at this revision.")
        chosen = set(observation_ids)
        cells = [c for c in cells if c["cell_id"] in chosen]
    if not cells:
        raise SpatialError("No observations match the declared export scope.")

    panel = summary["metadata"]["panel_genes"]
    selected_rna = rna_features or panel[: min(128, len(panel))]

    selected_assay = None
    selected_protein = []
    if assay_id:
        selected_assay = get_assay(project, assay_id)
        selected_protein = protein_features or [f["feature_id"] for f in selected_assay["features"]]

    ids = [c["cell_id"] for c in cells]
    identity = input_identity(project, rev_id, assay_id or "rna", selected_rna, selected_protein, ids)
    dest.mkdir(parents=True, exist_ok=True)

    # 1. metadata.json
    meta_path = dest / "input_contract.json"
    meta_path.write_text(
        _json(
            {
                "schema": SCHEMA,
                "input": identity,
                "observation_count": len(cells),
                "slice_id": summary["metadata"]["slice_id"],
                "units": summary["metadata"]["units"],
            }
        ),
        encoding="utf-8",
    )

    # 2. coordinates & labels
    obs_path = dest / "observations.csv"
    with obs_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["cell_id", "x", "y", "label", "included"])
        for c in cells:
            writer.writerow([c["cell_id"], c["x"], c["y"], c["label"], c["included"]])

    # 3. RNA counts
    rna_path = dest / "rna_counts.csv"
    with rna_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["cell_id", *selected_rna])
        for c in cells:
            writer.writerow([c["cell_id"], *[c["counts"].get(g, 0.0) for g in selected_rna]])

    # 4. Protein counts (if paired)
    if selected_assay and selected_protein:
        p_path = dest / "protein_counts.csv"
        paxis = [f["feature_id"] for f in selected_assay["features"]]
        pindex = [paxis.index(feat) for feat in selected_protein]
        with p_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["cell_id", *selected_protein])
            for cid in ids:
                vals = selected_assay["values"].get(cid, [None] * len(paxis))
                writer.writerow([cid, *[vals[idx] for idx in pindex]])

    return {
        "status": "exported",
        "destination": str(dest),
        "observation_count": len(cells),
        "observation_ids": ids,
        "rna_feature_count": len(selected_rna),
        "protein_feature_count": len(selected_protein),
        "input_identity": identity,
        "input": identity,
    }


def ingest_external_result(
    project: Project, result_file_or_dict: str | Path | dict, *, execution_evidence: dict | None = None
) -> dict[str, Any]:
    """Validate and register external integration results into the project."""
    if isinstance(result_file_or_dict, (str, Path)):
        path = Path(result_file_or_dict).resolve(strict=True)
        spec = json.loads(path.read_text(encoding="utf-8"))
    else:
        spec = dict(result_file_or_dict)
    spec.setdefault("output_origin", "model_prediction")
    method = spec.setdefault("method", {})
    if not isinstance(method.get("environment"), dict) or not method["environment"]:
        raise SpatialError(
            "Declare the actual external method environment; a bridge name is not execution evidence."
        )
    if type(method.get("uses_annotations")) is not bool:
        raise SpatialError("Explicitly declare whether the external fit used annotations.")
    if execution_evidence is not None:
        prov = spec.setdefault("provenance", {})
        prov["execution_evidence"] = execution_evidence
    registered = register_result(project, spec)
    return {
        "status": "ingested",
        "result_id": registered["object_id"],
        "object_sha256": registered["object_sha256"],
        "registered": registered,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="External pipeline runner and bridge.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    export_cmd = subparsers.add_parser("export", help="Export frozen input contract for external scripts.")
    export_cmd.add_argument("--project", required=True, help="Path to project directory.")
    export_cmd.add_argument("--output", required=True, help="Destination directory for inputs.")
    export_cmd.add_argument("--assay", help="Protein assay ID if paired.")
    export_cmd.add_argument("--revision", help="Revision ID to freeze.")

    ingest_cmd = subparsers.add_parser("ingest", help="Validate and register external integration output.")
    ingest_cmd.add_argument("--project", required=True, help="Path to project directory.")
    ingest_cmd.add_argument("--result", required=True, help="Path to result JSON.")

    args = parser.parse_args(argv)
    proj = Project(Path(args.project).resolve())

    if args.command == "export":
        res = export_input_contract(proj, args.output, revision_id=args.revision, assay_id=args.assay)
        print(json.dumps(res, indent=2))
    elif args.command == "ingest":
        res = ingest_external_result(proj, args.result)
        print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
