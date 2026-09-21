"""Explicit local file imports for the disk atlas and native image pyramid."""
import argparse
import json
from pathlib import Path

from .store import Project


def main(argv=None):
    parser = argparse.ArgumentParser(prog="spatial-collab large-data")
    parser.add_argument("operation", choices=["init", "import", "image", "verify"])
    parser.add_argument("--project", required=True)
    parser.add_argument("--spec", help="Explicit JSON file/column/coordinate contract")
    parser.add_argument("--atlas-id")
    parser.add_argument("--name", default="Spatial atlas workspace")
    args = parser.parse_args(argv)
    if args.operation == "init":
        project = Project.create_workspace(args.project, args.name)
        print(json.dumps(project.summary(), ensure_ascii=False, indent=2))
        return 0
    if args.operation in {"import", "image"} and not args.spec:
        parser.error("--spec is required for import/image")
    if args.operation == "verify" and not args.atlas_id:
        parser.error("--atlas-id is required for verify")
    project = Project(args.project)
    if args.operation == "verify":
        from .atlas import connect
        with connect(project, args.atlas_id, verify=True) as (_, record):
            result = {"atlas_id": args.atlas_id, "integrity": "verified", "totals": record["totals"]}
    else:
        spec = json.loads(Path(args.spec).read_text(encoding="utf-8-sig"))
        if args.operation == "import":
            from .atlas import build
            result = build(project, spec)
        else:
            from .pyramid import register
            result = register(project, **spec)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0
