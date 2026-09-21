"""Local CLI for explicit scientific input registration and executable workflows."""

import argparse
import json
from pathlib import Path

from .store import Project


def main(argv=None):
    parser = argparse.ArgumentParser(prog="spatial-collab analysis")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("catalog")
    for command in ("configure", "import", "snapshot", "run", "status", "result"):
        p = sub.add_parser(command)
        p.add_argument("--project", required=True)
        if command == "configure":
            p.add_argument("--python", required=True)
        elif command == "import":
            p.add_argument("path")
            p.add_argument("--sample-id", required=True)
            p.add_argument("--species", required=True, choices=["human", "mouse"])
            p.add_argument("--kind", required=True, choices=["spatial", "reference"])
            p.add_argument("--counts-layer", default="counts")
            p.add_argument("--label-key")
            p.add_argument("--batch-key")
            p.add_argument("--condition-key")
            p.add_argument("--feature-id-key", help="var column with stable gene IDs; default var index")
            p.add_argument("--feature-symbol-key", help="var column with gene symbols; _index for var index")
            p.add_argument("--spatial-key", default="spatial")
            p.add_argument(
                "--units", default="unknown", choices=["unknown", "array_index", "pixel", "micrometer"]
            )
        elif command == "snapshot":
            p.add_argument("--species", required=True, choices=["human", "mouse"])
            p.add_argument("--revision-id")
            p.add_argument("--observation-ids", help="JSON file containing explicit included IDs")
        elif command == "run":
            p.add_argument(
                "--spec", required=True, help="JSON recipe with method, input_ids, parameters, seed"
            )
            p.add_argument("--background", action="store_true")
        elif command == "status":
            p.add_argument("job_id")
        else:
            p.add_argument("result_id")
            p.add_argument("--output", help="Save complete result to a new JSON file instead of paging")
            p.add_argument("--offset", type=int, default=0)
    args = parser.parse_args(argv)
    from . import workflow_inputs as inputs, workflow_jobs as jobs, workflow_methods as methods

    if args.command == "catalog":
        result = methods.catalog()
    else:
        project = Project(args.project)
        if args.command == "configure":
            result = jobs.configure(project, args.python)
        elif args.command == "import":
            kwargs = {k: v for k, v in vars(args).items() if k not in {"command", "project", "path"}}
            record = inputs.import_h5ad(project, args.path, **kwargs)
            result = {k: record[k] for k in ("object_id", "shape", "kind", "species")}
        elif args.command == "snapshot":
            ids = (
                json.loads(Path(args.observation_ids).read_text(encoding="utf-8-sig"))
                if args.observation_ids
                else None
            )
            record = inputs.snapshot(
                project, args.revision_id or project.context()["head_revision"], args.species, ids
            )
            result = {k: record[k] for k in ("object_id", "shape", "kind", "species")}
        elif args.command == "run":
            spec = json.loads(Path(args.spec).read_text(encoding="utf-8-sig"))
            result = jobs.submit(project, spec, launch=args.background)
            if not args.background and result["status"] == "queued":
                result = jobs.execute(project, result["id"])
        elif args.command == "status":
            result = jobs.get(project, args.job_id)
        elif args.output:
            from . import objects

            record = objects.get(project, args.result_id, "analysisresult")
            with Path(args.output).open("x", encoding="utf-8") as stream:
                json.dump(record, stream, ensure_ascii=False, allow_nan=False)
            result = {"result_id": args.result_id, "output": str(Path(args.output).resolve())}
        else:
            result = jobs.inspect(project, args.result_id, offset=args.offset)
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
    return 1 if result.get("status") in {"failed", "cancelled"} else 0


if __name__ == "__main__":
    raise SystemExit(main())
