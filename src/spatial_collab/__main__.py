"""Command line for project import, local review and host transports."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(prog="spatial-collab", description="Spatial Collab research alpha")
    commands = parser.add_subparsers(dest="command", required=True)
    workflow = commands.add_parser("analysis", help="Executable algorithms, sparse input registration and durable analysis jobs")
    workflow.add_argument("arguments", nargs=argparse.REMAINDER)
    large = commands.add_parser("large-data", help="Stream large datasets and native TIFF pyramids into disk-backed atlases")
    large.add_argument("arguments", nargs=argparse.REMAINDER)
    commands.add_parser("formats", help="List lazy reader contracts, required inputs and limits")
    probe = commands.add_parser("probe", help="Inspect source filenames/headers without loading or importing")
    probe.add_argument("path")
    reader = commands.add_parser("import", help="Use an explicitly selected reader and a JSON options file")
    reader.add_argument("format_id")
    reader.add_argument("path")
    reader.add_argument("project")
    reader.add_argument("--options", required=True, help="JSON file with explicit reader arguments; see formats")
    benchmark = commands.add_parser("benchmark", help="Run developer-authored synthetic collaboration scenarios")
    benchmark.add_argument("output_dir")
    demo = commands.add_parser("demo", help="Create an explicitly synthetic 480-cell project")
    demo.add_argument("project")
    verify = commands.add_parser("verify-bundle", help="Verify and recompute an exported review bundle offline")
    verify.add_argument("bundle")
    verify.add_argument("--integrity-only", action="store_true")
    verify.add_argument("--include-workflows", action="store_true", help="Also refit task-specific algorithms; may require optional environment and long training")
    asset = commands.add_parser("register-asset", help="Register a local image/mask using an explicit source/frame/mapping JSON declaration")
    asset.add_argument("--project", required=True)
    asset.add_argument("--path", required=True)
    asset.add_argument("--options", required=True)
    for name in ("serve", "mcp", "inspect", "export"):
        sub = commands.add_parser(name)
        sub.add_argument("--project", default=os.environ.get("SPATIAL_COLLAB_PROJECT"))
        if name == "serve":
            sub.add_argument("--port", type=int, default=8765)
        elif name == "mcp":
            sub.add_argument("--transport", choices=["stdio", "streamable-http"], default="stdio")
            sub.add_argument("--port", type=int, default=8766)
        elif name == "export":
            sub.add_argument("--run-id")
            sub.add_argument("--compact", action="store_true", help="Export source plus overlays without redundant count copies")
    h5 = commands.add_parser("import-h5ad")
    h5.add_argument("path")
    h5.add_argument("project")
    h5.add_argument("--label-key", default="cell_type")
    h5.add_argument("--spatial-key", default="spatial")
    h5.add_argument("--counts-layer", default="counts")
    h5.add_argument("--slice-id", required=True)
    h5.add_argument("--slice-key")
    h5.add_argument("--coordinate-system", required=True)
    h5.add_argument("--units", required=True, choices=["micrometer", "pixel", "array_index", "unknown"])
    h5.add_argument("--name", default="Spatial transcriptomics")
    h5.add_argument("--biological-replicates", type=int, default=0, choices=[0, 1])
    xenium = commands.add_parser("import-xenium")
    xenium.add_argument("outs")
    xenium.add_argument("annotations")
    xenium.add_argument("project")
    xenium.add_argument("--slice-id", required=True)
    xenium.add_argument("--name", default="Xenium review")
    xenium.add_argument("--biological-replicates", type=int, default=0, choices=[0, 1])
    protein = commands.add_parser("register-protein", help="Attach a paired protein H5AD or wide CSV through explicit ID, sample and frame declarations")
    protein.add_argument("path")
    protein.add_argument("--project", required=True)
    protein.add_argument("--options", required=True)
    integration_import = commands.add_parser("import-integration", help="Register an external IntegrationResult JSON with exact source/axis checks")
    integration_import.add_argument("path")
    integration_import.add_argument("--project", required=True)
    integrate = commands.add_parser("integrate", help="Run bounded paired baselines, retaining single-modality controls")
    integrate.add_argument("--project", required=True)
    integrate.add_argument("--assay-id", required=True)
    integrate.add_argument("--revision-id")
    integrate.add_argument("--backend", choices=["balanced_pca", "smopca"], default="balanced_pca")
    integrate.add_argument("--components", type=int, default=10)
    integrate.add_argument("--clusters", type=int, default=6)
    integrate.add_argument("--seed", type=int, default=0)
    integrate.add_argument("--options", help="Optional exact feature/ROI/processing arguments in JSON")
    args = parser.parse_args(argv)
    try:
        if args.command == "large-data":
            from .large_data import main as large_main
            return large_main(args.arguments)
        if args.command == "analysis":
            from .workflows import main as workflow_main
            return workflow_main(args.arguments)
        if args.command == "import-integration":
            from .integration import register_result
            from .store import Project
            path = Path(args.path)
            if path.stat().st_size > 128 * 1024**2:
                raise ValueError("Integration result exceeds 128 MiB.")
            r = register_result(Project(args.project), json.loads(path.read_text(encoding="utf-8-sig")))
            result = {"result_id": r["object_id"], "object_sha256": r["object_sha256"]}
        elif args.command == "integrate":
            from .integration import run_baseline
            from .store import Project
            project = Project(args.project)
            options = json.loads(Path(args.options).read_text(encoding="utf-8-sig")) if args.options else {}
            r = run_baseline(project, args.revision_id or project.context()["head_revision"], args.assay_id,
                             backend=args.backend, components=args.components, clusters=args.clusters, seed=args.seed, **options)
            result = {"result_id": r["object_id"], "object_sha256": r["object_sha256"]}
        elif args.command == "register-protein":
            from .proteomics import register_protein
            from .store import Project
            option_path = Path(args.options)
            if option_path.stat().st_size > 65536:
                raise ValueError("Protein options exceed 64 KiB.")
            result = register_protein(Project(args.project), args.path, **json.loads(option_path.read_text(encoding="utf-8-sig")))
        elif args.command == "register-asset":
            from .assets import register_asset
            from .store import Project
            options_path = Path(args.options)
            if options_path.stat().st_size > 32 * 1024**2:
                raise ValueError("Asset declaration exceeds 32 MiB.")
            options = json.loads(options_path.read_text(encoding="utf-8-sig"))
            result = register_asset(Project(args.project), args.path, **options)
        elif args.command == "formats":
            from .import_registry import list_formats
            result = list_formats()
        elif args.command == "probe":
            from .import_registry import probe_source
            result = probe_source(args.path)
        elif args.command == "import":
            from .import_registry import import_source
            options_path = Path(args.options)
            if options_path.stat().st_size > 65536:
                raise ValueError("Import options exceed 64 KiB.")
            options = json.loads(options_path.read_text(encoding="utf-8-sig"))
            if not isinstance(options, dict):
                raise ValueError("Import options must be a JSON object.")
            result = import_source(args.format_id, args.path, args.project, **options).summary()
        elif args.command == "benchmark":
            from .benchmark import run_benchmark
            result = run_benchmark(args.output_dir)
            print(json.dumps(result, ensure_ascii=False, indent=2))
            return 0 if result["all_passed"] else 1
        elif args.command == "demo":
            from .demo import create_demo
            result = create_demo(args.project).summary()
        elif args.command == "verify-bundle":
            from .replay import verify_bundle
            result = verify_bundle(args.bundle, recompute=not args.integrity_only, recompute_workflows=args.include_workflows)
            print(json.dumps(result, ensure_ascii=False, indent=2))
            return 1 if result.get("recompute_matches") is False else 0
        elif args.command == "import-h5ad":
            from .importers import import_h5ad
            values = vars(args).copy()
            values.pop("command")
            values["destination"] = values.pop("project")
            result = import_h5ad(**values).summary()
        elif args.command == "import-xenium":
            from .importers import import_xenium
            values = vars(args).copy()
            values.pop("command")
            values["destination"] = values.pop("project")
            result = import_xenium(**values).summary()
        else:
            if not args.project:
                parser.error("Set --project or SPATIAL_COLLAB_PROJECT to one existing project directory.")
            from .store import Project
            project = Project(Path(args.project).resolve(strict=True))
            if args.command == "inspect":
                result = project.summary()
            elif args.command == "export":
                result = project.export_bundle(args.run_id, compact=args.compact)
            elif args.command == "serve":
                import uvicorn
                from .server import create_app
                uvicorn.run(create_app(project.root), host="127.0.0.1", port=args.port)
                return 0
            else:
                from .server import create_server
                server = create_server(project.root)
                server.settings.host = "127.0.0.1"
                server.settings.port = args.port
                server.run(transport=args.transport)
                return 0
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    except (ValueError, FileNotFoundError, FileExistsError, KeyError, ImportError, TypeError) as error:
        print(f"Spatial Collab: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
