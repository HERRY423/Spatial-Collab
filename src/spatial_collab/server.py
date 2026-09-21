"""One pinned project, exposed through MCP and the local review workbench."""

from __future__ import annotations

import hashlib
import inspect
import json
import logging
import math
import secrets
import threading
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit

from mcp.server.fastmcp import FastMCP
from mcp.server.transport_security import TransportSecuritySettings
from mcp.types import CallToolResult, TextContent, ToolAnnotations
from pydantic import ConfigDict, ValidationError, validate_call
from starlette.concurrency import run_in_threadpool
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request
from starlette.responses import HTMLResponse, JSONResponse
from starlette.routing import Route

from . import analysis, proteomics  # Load SciPy extensions on the server's main thread.
from .store import Project, SpatialError

UI_URI = "ui://spatial-collab/workbench.html"
MAX_VIEW_CELLS = 10_000
MAX_BODY_BYTES = 2_000_000
STATIC = Path(__file__).with_name("static")
LOG = logging.getLogger(__name__)
READ_TOOLS = {"open_project", "get_selection", "inspect_selection", "get_run", "get_revision",
              "get_context", "query_observations", "get_capabilities", "get_proposal", "inspect_quality",
              "get_feature_catalog", "list_assets", "inspect_asset", "get_asset_raster", "get_hypothesis_plan", "list_hypothesis_plans"}
READ_TOOLS.update({"get_overview", "list_assays", "inspect_protein"})
READ_TOOLS.update({"list_integrations", "get_integration", "inspect_integration", "compare_integrations", "get_integration_job", "get_multimodal_context"})
READ_TOOLS.add("compare_study")
READ_TOOLS.add("get_analysis_power_plan")
READ_TOOLS.add("get_integration_sensitivity")
READ_TOOLS.update({"get_analysis_catalog", "list_analysis_inputs", "get_analysis_job", "list_analysis_jobs", "list_analysis_results", "inspect_analysis_result"})
READ_TOOLS.update({"list_atlases", "get_atlas_view", "list_pyramids", "get_pyramid_tile", "list_study_analyses"})
READ_TOOLS.update({"list_projects", "get_transfer", "read_download"})
TOOL_NAMES = (
    "list_projects", "begin_upload", "append_upload", "complete_upload", "get_transfer",
    "discard_transfer", "import_uploaded_project", "prepare_download", "read_download",
    "create_analysis_power_plan", "get_analysis_power_plan",
    "list_atlases", "get_atlas_view", "share_atlas_view", "freeze_atlas_roi", "list_pyramids", "get_pyramid_tile", "register_study_design", "run_study_inference", "list_study_analyses",
    "open_project", "get_selection", "set_selection", "inspect_selection",
    "propose_revision", "apply_revision", "revert_revision", "run_comparison",
    "get_run", "get_revision", "export_review_bundle", "get_context", "query_observations",
    "get_capabilities", "run_region_comparison", "get_proposal", "inspect_quality", "run_sensitivity",
    "get_feature_catalog", "list_assets", "inspect_asset", "get_asset_raster", "create_hypothesis_plan", "revise_hypothesis_plan",
    "freeze_hypothesis_plan", "get_hypothesis_plan", "list_hypothesis_plans", "run_hypothesis_plan",
    "get_overview", "list_assays", "inspect_protein", "run_protein_comparison",
    "register_integration", "list_integrations", "get_integration", "inspect_integration", "compare_integrations",
    "submit_integration", "get_integration_job", "cancel_integration_job", "retry_integration_job", "filter_integration",
    "register_multimodal_object", "get_multimodal_context",
    "compare_study",
    "run_integration_sensitivity", "get_integration_sensitivity",
    "get_analysis_catalog", "prepare_analysis_input", "list_analysis_inputs", "submit_analysis",
    "get_analysis_job", "list_analysis_jobs", "cancel_analysis_job", "retry_analysis_job", "list_analysis_results", "inspect_analysis_result",
)


class ToolService:
    """Shared dispatch. No tool accepts an input/output filesystem path."""

    def __init__(self, project_root: str | Path):
        self.project = Project(Path(project_root).resolve())
        self.call_lock = threading.RLock()
        self._calls = {
            name: validate_call(config=ConfigDict(strict=True))(getattr(self, name))
            for name in TOOL_NAMES
        }

    def call(self, name: str, arguments: dict | None = None) -> dict:
        if name not in self._calls:
            raise SpatialError("Unknown tool")
        if arguments is not None and not isinstance(arguments, dict):
            raise SpatialError("Tool arguments must be an object")
        arguments = dict(arguments or {})
        project_id = arguments.pop("project_id", None)
        if project_id is not None:
            if not isinstance(project_id, str):
                raise SpatialError("project_id must be a string.")
            from .transfers import projects
            matches = [p for p in projects(self.project) if p.summary()["project_id"] == project_id]
            if not matches:
                raise SpatialError("Project is not available to this connection.")
            if matches[0].root != self.project.root:
                return ToolService(matches[0].root).call(name, arguments)
        try:
            result = self._calls[name](**arguments)
            result.setdefault("project_id", self.project.summary()["project_id"])
            return result
        except ValidationError as exc:
            issues = [f"{'.'.join(str(x) for x in e['loc'])}: {e['msg']}" for e in exc.errors()]
            raise SpatialError("Invalid arguments: " + "; ".join(issues)) from None

    def list_projects(self) -> dict:
        """List only this connection's primary project and uploaded projects. Pin project_id on subsequent calls; never accept filesystem paths."""
        from .transfers import projects
        return {"projects": [{"project_id": p.summary()["project_id"], "name": p.summary()["metadata"]["name"],
                              "cell_count": p.summary()["cell_count"]} for p in projects(self.project)]}

    def begin_upload(self, filename: str, size: int, sha256: str) -> dict:
        """Reserve a bounded upload after the user selects a JSON snapshot or H5AD. Do not place file bytes in the conversation; use the workbench file picker."""
        from .transfers import Transfers
        return Transfers(self.project).begin(filename, size, sha256)

    def append_upload(self, file_id: str, offset: int, data_base64: str) -> dict:
        """Append a checked 256 KiB upload chunk; same-byte retries are safe. Intended for the UI, not model-generated data."""
        from .transfers import Transfers
        return Transfers(self.project).append(file_id, offset, data_base64)

    def complete_upload(self, file_id: str) -> dict:
        """Verify complete length and SHA256 before an explicit import. Upload does not change existing project data."""
        from .transfers import Transfers
        return Transfers(self.project).complete(file_id)

    def get_transfer(self, file_id: str) -> dict:
        """Read transfer size, checksum, state, resume offset and expiry without exposing local paths."""
        from .transfers import Transfers
        return Transfers(self.project, create=False).status(file_id)

    def discard_transfer(self, file_id: str) -> dict:
        """Delete only a temporary uploaded/download copy; imported projects and original review exports are retained."""
        from .transfers import Transfers
        return Transfers(self.project).discard(file_id)

    def import_uploaded_project(self, file_id: str, options: dict) -> dict:
        """Validate an uploaded JSON snapshot or H5AD into a NEW project. H5AD needs explicit counts_layer, slice_id, coordinate_system and units. Never overwrites the primary project. Use returned project_id on every subsequent call."""
        from .transfers import Transfers
        return Transfers(self.project).import_project(file_id, options)

    def prepare_download(self, run_id: str | None = None, compact: bool = True) -> dict:
        """Create a checksummed ZIP and a project-scoped download handle. Use the workbench Download button to retrieve bytes; a server path is not a ChatGPT attachment."""
        from .transfers import Transfers
        return Transfers(self.project).export(run_id, compact)

    def read_download(self, file_id: str, offset: int = 0, length: int = 262144) -> dict:
        """Read one bounded, project-owned ZIP chunk for the workbench download bridge."""
        from .transfers import Transfers
        return Transfers(self.project, create=False).read(file_id, offset, length)

    def open_project(self, revision_id: str | None = None, bounds: list[float] | None = None) -> dict:
        """Open shared workbench with exact revision and centroid view. At most 10000 cells; narrow bounds [xmin,ymin,xmax,ymax] if larger. No silent sampling."""
        summary = self.project.summary()
        revision_id = revision_id or summary["head_revision"]
        cells = self.project.cells(revision_id)
        full_bounds = None
        if cells:
            full_bounds = [min(c["x"] for c in cells), min(c["y"] for c in cells),
                           max(c["x"] for c in cells), max(c["y"] for c in cells)]
        if bounds is not None:
            if len(bounds) != 4 or not all(math.isfinite(v) for v in bounds):
                raise SpatialError("bounds must contain four finite numbers")
            if bounds[0] >= bounds[2] or bounds[1] >= bounds[3]:
                raise SpatialError("bounds must have positive width and height")
            visible = [c for c in cells if bounds[0] <= c["x"] <= bounds[2]
                       and bounds[1] <= c["y"] <= bounds[3]]
        else:
            visible = cells
        count = len(visible)
        limited = count > MAX_VIEW_CELLS
        points = [] if limited else [{k: c[k] for k in ("cell_id", "x", "y", "label", "included", "region")}
                                    for c in visible]
        from .exploration import semantics
        from .coordinates import coordinate_capabilities
        # Keep the complete feature axis on disk and in scientific operations;
        # agents and the viewer request only the page they need.
        host_metadata = dict(summary["metadata"])
        features = host_metadata.pop("features", None)
        if features is not None:
            from collections import Counter
            symbols = Counter(row["symbol"] for row in features if row["symbol"] is not None)
            host_metadata["feature_catalog"] = {"feature_count": len(features),
                "duplicate_symbol_count": sum(n > 1 for n in symbols.values()),
                "tool": "get_feature_catalog", "mapping_in_response": False}
        panel = host_metadata.get("panel_genes", [])
        host_metadata["panel_gene_count"] = len(panel)
        if len(panel) > 500:
            host_metadata.pop("panel_genes")
            host_metadata["panel_genes_omitted"] = True
        return {
            **summary,
            "metadata": host_metadata,
            "semantics": semantics(summary["metadata"]),
            "coordinate_capabilities": coordinate_capabilities(summary["metadata"]),
            "view": {"revision_id": revision_id, "cells": points, "total_cells": len(cells),
                     "matching_cells": count, "returned_cells": len(points), "complete": not limited,
                     "bounds": bounds or full_bounds, "full_bounds": full_bounds, "limit": MAX_VIEW_CELLS,
                     "notice": "View exceeds 10000 cells. Zoom in; no cells were sampled." if limited else "Centroids only; polygon selection uses centroid containment."},
            "scientific_authorization": "NOT_ESTABLISHED",
        }

    def create_analysis_power_plan(self, spec: dict, assumptions: dict) -> dict:
        """Freeze model/graph/family/rank assumptions before Moran/LR analysis; report exact permutation resolution and conditional simulated MDE, never universal or observed power."""
        from .power import create_configured
        return create_configured(self.project, spec, assumptions)

    def get_analysis_power_plan(self, plan_id: str) -> dict:
        """Read a frozen power declaration and uncertainty, retaining prospective/retrospective local timing."""
        from .objects import get
        plan = get(self.project, plan_id, "powerplan")
        return {**plan, "declaration": get(self.project, plan["design_id"], "powerdesign")}

    def get_analysis_catalog(self) -> dict:
        """List executable algorithms by scientific task; availability is package discovery, not validated biology."""
        from .workflow_methods import catalog
        from .workflow_jobs import runtime
        result, rt = catalog(), runtime(self.project, refresh=True)
        if "packages" in rt:
            packages = {k.lower().replace("_", "-"): v for k, v in rt["packages"].items()}
            for method in result["methods"]:
                version = packages.get(method["package"].lower().replace("_", "-"))
                method.update(available=version is not None, version=version)
        return {**result, "configured_runtime": rt}

    def list_atlases(self) -> dict:
        """List immutable disk-backed spatial datasets, preserved molecule and boundary layers."""
        from .atlas import catalog
        return catalog(self.project)

    def get_atlas_view(self, atlas_id: str, bounds: list[float] | None = None, layer: str = "cells", feature: str | None = None, limit: int = 10000, grid: int = 64) -> dict:
        """Query indexed viewport. Dense views return complete count aggregates, not sampled cells; zoom for exact IDs."""
        from .atlas import view
        return view(self.project, atlas_id, bounds, layer, feature, limit, grid)

    def freeze_atlas_roi(self, atlas_id: str, bounds: list[float], max_cells: int = 10000) -> dict:
        """Create an exact, hash-bound sparse analysis input from an atlas ROI; no implicit sampling."""
        from .atlas import extract
        result = extract(self.project, atlas_id, bounds, max_cells)
        return {"input_id": result["object_id"], "shape": result["shape"], "provenance": result["provenance"]}

    def share_atlas_view(self, atlas_id: str, bounds: list[float], layer: str = "cells", feature: str | None = None, image_id: str | None = None) -> dict:
        """Share the last active immutable atlas/viewport/layer with get_context; does not select cells or change annotations."""
        from .atlas import share_view
        return share_view(self.project, atlas_id, bounds, layer, feature, image_id)

    def list_pyramids(self, atlas_id: str) -> dict:
        """List native tiled image pyramids explicitly registered to this atlas frame."""
        from .pyramid import catalog
        return catalog(self.project, atlas_id)

    def get_pyramid_tile(self, image_id: str, level: int, x: int, y: int) -> dict:
        """Decode only intersecting TIFF source tiles for one 256-pixel image tile, with explicit world corners."""
        from .pyramid import tile
        return tile(self.project, image_id, level, x, y)

    def register_study_design(self, spec: dict) -> dict:
        """Freeze researcher-declared subject/sample/section identities, conditions, batches and ROI origins."""
        from .study import register_study
        return register_study(self.project, spec)

    def run_study_inference(self, study_id: str, records: list[dict], plan: dict) -> dict:
        """Fit subject-level Welch/paired tests or batch-adjusted subject-random-intercept LMM; preserve missing/failed tests and BH family."""
        from .study_inference import infer
        return infer(self.project, study_id, records, plan)

    def list_study_analyses(self) -> dict:
        """Read frozen study designs and inference results including failures and missing subjects."""
        from . import objects
        return {kind: [objects.get(self.project, oid, kind) for oid in objects.catalog(self.project, kind)] for kind in ("study", "studyresult")}

    def prepare_analysis_input(self, revision_id: str, species: str, observation_ids: list[str] | None = None) -> dict:
        """Freeze exact included raw RNA counts, IDs, species and coordinates for analysis. Explicit ROI never silently intersects."""
        from .workflow_inputs import snapshot
        result = snapshot(self.project, revision_id, species, observation_ids)
        return {"input_id": result["object_id"], "shape": result["shape"], "species": result["species"]}

    def list_analysis_inputs(self) -> dict:
        """Read registered sparse raw-count input metadata, including spatial samples and annotated references."""
        from .workflow_inputs import list_inputs
        return {"inputs": list_inputs(self.project)}

    def submit_analysis(self, spec: dict) -> dict:
        """Run MOFA/MEFISTO, deconvolution, Harmony, PASTE, spatial domains, SVG, spatial LR or PROGENy in a separate worker. Inspect catalog and input contracts first."""
        from .workflow_jobs import submit
        return submit(self.project, spec)

    def get_analysis_job(self, job_id: str) -> dict:
        """Read saved analysis job status, recipe and concrete failures."""
        from .workflow_jobs import get
        return get(self.project, job_id)

    def list_analysis_jobs(self) -> dict:
        """Resume the latest 50 durable analysis tasks after reconnecting; retain failures and cancellations."""
        from .workflow_jobs import list_jobs
        return list_jobs(self.project)

    def cancel_analysis_job(self, job_id: str) -> dict:
        """Request cooperative cancellation at method boundaries. A training phase can finish before stopping."""
        from .workflow_jobs import cancel
        return cancel(self.project, job_id)

    def retry_analysis_job(self, job_id: str) -> dict:
        """Retry a failed or cancelled analysis while retaining the earlier attempt."""
        from .workflow_jobs import retry
        return retry(self.project, job_id)

    def list_analysis_results(self) -> dict:
        """List computed task-specific outputs; they are not automatically annotation revisions."""
        from .workflow_jobs import results
        return results(self.project)

    def inspect_analysis_result(self, result_id: str, offset: int = 0, limit: int = 100,
                                field: str | None = None, input_index: int = 0) -> dict:
        """Page exact result rows, spatial positions and semantics. Only same-current-source rows may directly become a shared selection."""
        from .workflow_jobs import inspect
        return inspect(self.project, result_id, offset, limit, field, input_index)

    def get_selection(self, selection_id: str | None = None) -> dict:
        """Read exact server-persisted cell IDs, geometry, revision and stale status; optional historical selection ID."""
        summary = self.project.summary()
        selection = self.project.get_selection(selection_id) if selection_id else summary["selection"]
        return {"selection": selection, "head_revision": summary["head_revision"]}

    def get_overview(self) -> dict:
        """Overview of RNA and paired protein layers, shared selection, revision and recent saved analyses."""
        summary = self.project.summary()
        return {"project_id": summary["project_id"], "name": summary["metadata"]["name"],
                "sample_id": summary["metadata"].get("sample_id"), "source_sha256": summary["source_sha256"],
                "head_revision": summary["head_revision"], "observation_count": summary["cell_count"],
                "observation_unit": summary["metadata"].get("observation_unit", "cell"),
                "units": summary["metadata"]["units"], "rna_feature_count": len(summary["metadata"]["panel_genes"]),
                "protein_assays": proteomics.list_assays(self.project)["assays"], "selection_count": (summary["selection"] or {}).get("cell_count", 0),
                "revision_count": summary["revision_count"], "recent_runs": summary["runs"][-12:],
                "scope": "Empty atlas container; datasets have independent identities." if summary["metadata"]["source_kind"] == "atlas_workspace" else "One editable primary spatial sample plus independently identified atlases; no cross-sample integration inferred.",
                "scientific_authorization": "NOT_ESTABLISHED"}

    def list_assays(self) -> dict:
        """List registered protein assays and measured features without exporting their full matrices."""
        return proteomics.list_assays(self.project)

    def inspect_protein(self, assay_id: str, revision_id: str, feature: str, selection_id: str | None = None,
                        scale: str = "raw", cofactor: float = 5.0, bounds: list[float] | None = None,
                        limit: int = 10000, offset: int = 0) -> dict:
        """Read an exact protein channel at a shared revision/ROI. Missing and zero differ; raw/log1p/asinh are explicit display transforms."""
        return proteomics.inspect_protein(self.project, assay_id, revision_id, feature, selection_id, scale, cofactor, bounds, limit, offset)

    def register_integration(self, spec: dict) -> dict:
        """Register external integration outputs with exact project, source, revision, feature and observation identities. No code execution."""
        from .integration import register_result
        r = register_result(self.project, spec)
        return {"result_id": r["object_id"], "object_sha256": r["object_sha256"]}

    def list_integrations(self) -> dict:
        """List retained integration results, fit/filter mode and historical status."""
        from .integration import list_results
        return list_results(self.project)

    def get_integration(self, result_id: str) -> dict:
        """Read method, inputs and preprocessing without sending full embeddings to the model."""
        from .integration import get_result
        r = get_result(self.project, result_id)
        return {k: v for k, v in r.items() if k not in {"representations", "domains", "fitted_models", "observation_ids"}}

    def inspect_integration(self, result_id: str, partition: str = "joint", bounds: list[float] | None = None,
                            limit: int = 2000, offset: int = 0) -> dict:
        """Page exact method domain points in the shared spatial coordinates; no display sampling."""
        from .integration import inspect_result
        return inspect_result(self.project, result_id, partition, bounds, limit, offset)

    def compare_integrations(self, left_id: str, right_id: str | None = None, left_partition: str = "rna_only",
                             right_partition: str = "joint", limit: int = 100, offset: int = 0) -> dict:
        """Compare partitions using permutation-invariant co-membership and rank boundary disagreement for review, never automatic correction."""
        from .integration import compare_results
        return compare_results(self.project, left_id, right_id, left_partition, right_partition, limit, offset)

    def submit_integration(self, spec: dict) -> dict:
        """Queue paired RNA-only/protein-only/joint baseline or optional SMOPCA on an explicit revision. Worker is separate; missing measurements fail."""
        from .integration_jobs import submit
        return submit(self.project, spec)

    def get_integration_job(self, job_id: str) -> dict:
        """Read durable job status, recipe, cache key, result or failure."""
        from .integration_jobs import get_job
        return get_job(self.project, job_id)

    def cancel_integration_job(self, job_id: str) -> dict:
        """Cancel pending work cooperatively; active numerical phase may finish before stopping."""
        from .integration_jobs import cancel
        return cancel(self.project, job_id)

    def retry_integration_job(self, job_id: str) -> dict:
        """Create a new attempt for a failed/cancelled job, preserving the original record."""
        from .integration_jobs import retry
        return retry(self.project, job_id)

    def filter_integration(self, result_id: str, revision_id: str) -> dict:
        """Retain original fitted embeddings on included observations at a new revision. Explicit fixed-model filter, not refit."""
        from .integration import fixed_filter
        r = fixed_filter(self.project, result_id, revision_id)
        return {"result_id": r["object_id"], "mode": r["mode"]}

    def register_multimodal_object(self, kind: str, spec: dict) -> dict:
        """Register explicit assay descriptors, correspondences, derived layers, molecular relations or study designs; no inferred identity."""
        from . import multimodal
        from .study import register_study
        handlers = {"assay": multimodal.register_descriptor, "correspondence": multimodal.register_correspondence,
                    "layer": multimodal.register_layer, "molecular_relation": multimodal.register_molecular_relation,
                    "study": register_study}
        if kind not in handlers:
            raise SpatialError("Unknown multimodal object kind.")
        result = handlers[kind](self.project, spec)
        return {"object_id": result["object_id"], "object_sha256": result["object_sha256"]}

    def get_multimodal_context(self) -> dict:
        """List explicit assay relationships, layer semantics and study designs without their numeric matrices."""
        from . import objects
        result = {}
        for kind in ("assaydescriptor", "correspondence", "measurementlayer", "molecularrelation", "study"):
            result[kind] = [{k: v for k, v in objects.get(self.project, oid, kind).items() if k not in {"values", "edges", "observation_ids"}}
                            for oid in objects.catalog(self.project, kind)]
        return result

    def compare_study(self, study_id: str, records: list[dict]) -> dict:
        """Compare declared same-metric ROI summaries at subject level, retaining batch confounding and missing sections; no group hypothesis test."""
        from .study import compare_study
        return compare_study(self.project, study_id, records)

    def run_integration_sensitivity(self, plan_id: str, version: int, assay_id: str, alternatives: list[dict]) -> dict:
        """Extend an existing frozen hypothesis plan with all declared preprocessing/backend alternatives, retaining unknowns/failures. No annotation revision or human approval created."""
        from .integration_sensitivity import run_plan
        return run_plan(self.project, plan_id, version, assay_id, alternatives)

    def get_integration_sensitivity(self, run_id: str) -> dict:
        """Read all live job outcomes from an integration sensitivity experiment; no winning-only filtering."""
        from .integration_sensitivity import inspect_plan_run
        return inspect_plan_run(self.project, run_id)

    def run_protein_comparison(self, assay_id: str, base_revision: str, target_revision: str, selection_id: str,
                               features: list[str], background_selection_id: str | None = None,
                               scale: str = "raw", cofactor: float = 5.0, rna_gene: str | None = None) -> dict:
        """Compare frozen foreground/background protein signals before/after shared annotation/inclusion revisions; optional exact paired RNA Spearman concordance. Descriptive only."""
        return proteomics.compare_protein_regions(self.project, assay_id, base_revision, target_revision, selection_id,
                                       features, background_selection_id, scale, cofactor, rna_gene)

    def set_selection(self, expected_revision: str, cell_ids: list[str] | None = None,
                      polygon: list[list[float]] | None = None, name: str = "Selection") -> dict:
        """Persist exact IDs OR a polygon in the project's recorded coordinate frame and units, including boundary. Unknown physical scale permits selection, not micrometer distances. Requires current revision."""
        return self.project.set_selection(expected_revision, cell_ids, polygon, name)

    def inspect_selection(self, genes: list[str] | None = None) -> dict:
        """Summarize selected labels and measured raw marker counts. Unmeasured panel genes differ from measured zeros; no automated cell identity decision."""
        result = self.project.inspect_selection(genes)
        measured_ids = {(item.get("feature_resolution") or {}).get("feature_id", query)
                        for query, item in result["expression"].items() if item["measured"]}
        result["cells"] = [{**cell, "counts": {key: value for key, value in cell["counts"].items()
                                              if key in measured_ids}} for cell in result["cells"]]
        result["cell_counts_scope"] = "Requested measured features only; full counts remain in immutable source."
        result["cell_counts_feature_ids"] = sorted(measured_ids)
        return result

    def propose_revision(self, expected_revision: str, selection_id: str, changes: dict,
                         rationale: str, actor: str = "agent") -> dict:
        """Preview exact label/included/region overlay and affected results without committing. Supply scientific rationale and current selection ID."""
        return self.project.propose_revision(expected_revision, selection_id, changes, rationale, actor)

    def apply_revision(self, proposal_id: str, expected_revision: str, reviewer: str, confirmation: bool) -> dict:
        """Commit an already reviewed proposal ONLY after explicit researcher confirmation. Record their supplied name; this is attribution, not authenticated scientific authority."""
        return self.project.apply_revision(proposal_id, expected_revision, reviewer, confirmation)

    def revert_revision(self, target_revision: str, expected_revision: str, reviewer: str, confirmation: bool) -> dict:
        """After explicit researcher confirmation, restore target overlays as a NEW revision. Keeps all history and marks affected results stale."""
        return self.project.revert_revision(target_revision, expected_revision, reviewer, confirmation)

    def run_comparison(self, base_revision: str, target_revision: str, selection_id: str | None = None,
                       radius_um: float = 35.0, graph_scope: str = "roi_induced", source_label: str = "T cell",
                       target_label: str = "Myeloid", min_effect: float = 0.1) -> dict:
        """Rerun identical descriptive composition and directed neighbor metrics before/after. Frozen ROI; declared radius, graph scope and excess-over-abundance threshold. No biological significance or causality."""
        return analysis.compare(self.project, base_revision, target_revision, selection_id, radius_um,
                                graph_scope, source_label, target_label, min_effect)

    def get_run(self, run_id: str) -> dict:
        """Read a saved comparison, its exact parameters, denominators, limitations and current stale status."""
        return self.project.get_run(run_id)

    def get_revision(self, revision_id: str) -> dict:
        """Read an immutable revision with rationale, attributed reviewer and exact overlay delta."""
        return self.project.get_revision(revision_id)

    def export_review_bundle(self, run_id: str | None = None, compact: bool = False) -> dict:
        """Write a uniquely named review bundle under this project's exports, including hashes and provenance. Hash integrity is not scientific validation."""
        return self.project.export_bundle(run_id, compact=compact)

    def get_context(self, after_cursor: str | None = None) -> dict:
        """Poll lightweight shared state without reloading expression data. changed=false permits waiting; a change requires rereading affected objects. CAS still protects all edits."""
        return self.project.context(after_cursor)

    def get_proposal(self, proposal_id: str) -> dict:
        """Read the exact saved proposal, changed IDs and rationale for researcher review in another host. Reports stale/applied state; never grants approval."""
        return self.project.get_proposal(proposal_id)

    def get_capabilities(self) -> dict:
        """Discover format/host contracts and research operations before planning a task. Reader availability is not proof that optional dependencies are installed."""
        from .import_registry import list_formats
        from .coordinates import coordinate_capabilities
        return {"formats": list_formats(), "tools": list(TOOL_NAMES),
                "coordinate_capabilities": coordinate_capabilities(self.project.summary()["metadata"]),
                "hosts": ["MCP stdio", "MCP Streamable HTTP", "MCP Apps compatible viewer", "local web", "optional napari npe2"],
                "workflow": ["inspect explicit units and observation semantics", "query exact observations with pinned pagination",
                             "save shared selection", "inspect measured markers and source quality attributes", "evaluate explicit non-committing hypotheses",
                             "propose overlay", "researcher confirms concrete delta",
                             "compare identical frozen objects", "inspect paired protein channels and RNA concordance",
                             "export and independently recompute locally"],
                "protein": {"registration": "Local CLI register-protein; H5AD or wide CSV", "pairing": "exact sample/original ID/XY",
                            "measurement_types": ["antibody_count", "intensity"], "scales": ["raw", "log1p", "asinh"],
                            "max_features_per_assay": 512, "max_comparison_features": 32, "cross_sample_integration": False},
                "scientific_analysis": {"catalog_tool": "get_analysis_catalog", "runner": "submit_analysis",
                    "inputs": "immutable raw-count snapshots; explicitly registered spatial samples and annotated references",
                    "multi_input_tasks": ["reference deconvolution", "expression batch correction", "full-overlap slice correspondence"],
                    "result_review": "inspect_analysis_result", "resume": "list_analysis_jobs", "model_outputs_commit_annotations": False},
                "large_data": {"catalog": "list_atlases", "viewport": "get_atlas_view", "shared_view": "share_atlas_view", "roi_analysis": "freeze_atlas_roi", "image_tiles": "get_pyramid_tile", "storage": "SQLite R-tree and sparse counts", "dense_views": "all-object density aggregates, not sampling"},
                "study_inference": {"design": "register_study_design", "runner": "run_study_inference", "methods": ["welch", "paired_t", "mixedlm"], "independent_unit": "subject", "mixedlm_inference": "asymptotic_Wald"},
                "limits": {"single_slice": True, "single_slice_meaning": "One editable primary slice; scientific analyses and indexed atlases may use explicit additional inputs.", "centroids_only": False, "scientific_authorization": "NOT_ESTABLISHED",
                           "remote_chatgpt_requires_deployment": True, "host_conversation_acceptance": "not_established"}}

    def query_observations(self, revision_id: str, labels: list[str] | None = None, included: bool | None = None,
                           bounds: list[float] | None = None, gene: str | None = None, min_count: float = 0.0,
                           limit: int = 100, cursor: str | None = None, attribute_filter: dict | None = None) -> dict:
        """Discover exact IDs by labels, inclusion, bounds, measured raw marker count or attribute_filter={field,operator,value}. Attribute operators eq/ne/lt/le/gt/ge exclude unknown/null values, never converting them to zero. Pages bind revision and all filters; no selection change."""
        from .exploration import query_observations
        return query_observations(self.project, revision_id, labels=labels, included=included, bounds=bounds,
                                  gene=gene, min_count=min_count, limit=limit, cursor=cursor, attribute_filter=attribute_filter)

    def run_region_comparison(self, base_revision: str, target_revision: str, selection_id: str,
                              background_selection_id: str | None = None, genes: list[str] | None = None) -> dict:
        """Compare a frozen ROI against an explicit disjoint background or its whole-slice complement, before/after revision. Reports label fractions, marker detection/raw/panel-normalized means and denominators. Descriptive, no DE significance or cell identity verdict."""
        from .exploration import compare_regions
        return compare_regions(self.project, base_revision, target_revision, selection_id,
                               background_selection_id=background_selection_id, genes=genes)

    def inspect_quality(self, revision_id: str, selection_id: str | None = None) -> dict:
        """Inspect imported measurement-quality attributes and expression-library summaries at a pinned revision/selection. Missing quality fields stay unknown. No automatic doublet verdict, QC pass/fail or exclusion."""
        from .quality import inspect_quality
        return inspect_quality(self.project, revision_id, selection_id)

    def run_sensitivity(self, revision_id: str, selection_id: str, variants: list[dict], radii_um: list[float],
                        source_label: str, target_label: str, graph_scope: str = "whole_slice", min_effect: float = 0.1) -> dict:
        """Compute independent what-if label/inclusion overlays at declared radii WITHOUT committing annotation revisions or claiming human approval. Each variant has name, cell_ids, changes and rationale. Saves a replayable hypothesis result. Whole_slice means imported universe, possibly cropped."""
        from .sensitivity import run_sensitivity
        return run_sensitivity(self.project, revision_id, selection_id, variants, radii_um,
                               source_label, target_label, graph_scope, min_effect)

    def get_feature_catalog(self, query: str = "", limit: int = 100, offset: int = 0) -> dict:
        """Review stable feature IDs and original symbols. Duplicate symbols stay distinct; use feature_id:<ID> to disambiguate. This paged catalog never merges counts."""
        from .identity import validate_features
        if not 1 <= limit <= 500 or offset < 0 or len(query) > 500:
            raise SpatialError("Use limit 1..500, nonnegative offset and a bounded query.")
        summary = self.project.summary()
        features = validate_features(summary["metadata"])
        from collections import Counter
        counts = Counter(item["symbol"] for item in features if item["symbol"] is not None)
        matching = [dict(item, symbol_ambiguous=counts[item["symbol"]] > 1) for item in features
                    if query.casefold() in item["feature_id"].casefold() or query.casefold() in (item["symbol"] or "").casefold()]
        page = sorted(matching, key=lambda item: item["feature_id"])[offset:offset + limit]
        return {"source_sha256": summary["source_sha256"], "features": page, "matching_count": len(matching),
                "offset": offset, "next_offset": offset + len(page) if offset + len(page) < len(matching) else None,
                "counts_merged": False, "identity_scope": summary["metadata"].get("identity_scope", "project_local")}

    def list_assets(self) -> dict:
        """List locally registered, source-bound images and masks. Availability is not biological registration correctness; files are checked on inspection/loading."""
        from .assets import list_assets
        return list_assets(self.project)

    def inspect_asset(self, asset_id: str, cell_ids: list[str] | None = None) -> dict:
        """Inspect up to 100 exact objects against an explicitly registered image or mask. Coordinates, label IDs and identity correspondence are returned; no merge/coexpression verdict."""
        from .assets import inspect_asset
        return inspect_asset(self.project, asset_id, cell_ids)

    def get_asset_raster(self, asset_id: str, max_dimension: int = 1024) -> dict:
        """Render a verified registered image or segmentation mask plane to a bounded base64 PNG data URL with affine bounds."""
        from .assets import get_asset_raster
        return get_asset_raster(self.project, asset_id, max_dimension=max_dimension)

    def create_hypothesis_plan(self, spec: dict) -> dict:
        """Save an editable draft plan with explicit selectors, changes, radii, graph scopes and background. No annotation change or scientific approval."""
        from .hypotheses import create_hypothesis_plan
        return create_hypothesis_plan(self.project, spec)

    def revise_hypothesis_plan(self, plan_id: str, expected_version: int, spec: dict) -> dict:
        """Append a new draft version; rejects stale concurrent edits. Previous frozen plans remain immutable."""
        from .hypotheses import revise_hypothesis_plan
        return revise_hypothesis_plan(self.project, plan_id, expected_version, spec)

    def freeze_hypothesis_plan(self, plan_id: str, expected_version: int) -> dict:
        """Freeze exact IDs and all declared settings before execution. Unknown selectors remain visible. Freezing is computational prerecording, not independent preregistration or researcher approval."""
        from .hypotheses import freeze_hypothesis_plan
        return freeze_hypothesis_plan(self.project, plan_id, expected_version)

    def get_hypothesis_plan(self, plan_id: str, version: int | None = None) -> dict:
        """Read an exact draft/frozen plan version for cross-host review."""
        from .hypotheses import get_hypothesis_plan
        return get_hypothesis_plan(self.project, plan_id, version)

    def list_hypothesis_plans(self) -> dict:
        """List saved hypothesis plans; older frozen versions remain retrievable by ID/version."""
        from .hypotheses import list_hypothesis_plans
        return list_hypothesis_plans(self.project)

    def run_hypothesis_plan(self, plan_id: str, version: int | None = None) -> dict:
        """Run every frozen hypothesis/radius/scope with both edge-weighted and source-equal-weighted methods. Keep unknowns/failures and rank all differences. No annotation commits."""
        from .hypotheses import run_hypothesis_plan
        return run_hypothesis_plan(self.project, plan_id, version)


def workbench_script() -> str:
    return "\n".join((STATIC / name).read_text(encoding="utf-8") for name in ("workbench.js", "transfer-ui.js"))


def workbench_html(csrf_token: str = "") -> str:
    html = (STATIC / "workbench.html").read_text(encoding="utf-8")
    return html.replace("/*__STYLE__*/", (STATIC / "workbench.css").read_text(encoding="utf-8")).replace(
        "/*__SCRIPT__*/", workbench_script()).replace("__CSRF_TOKEN__", csrf_token)


def _safe_result(service: ToolService, name: str, arguments: dict) -> CallToolResult:
    try:
        with service.call_lock:
            result = service.call(name, arguments)
        # All exact IDs remain in structuredContent; text is bounded to avoid echoing thousands of points.
        summary = {k: v for k, v in result.items() if k not in {"view", "cells"}}
        text = json.dumps(summary, ensure_ascii=False, allow_nan=False)
        if len(text) > 24_000:
            text = json.dumps({"tool": name, "message": "Result available in structuredContent", "head_revision": result.get("head_revision"), "selection_id": result.get("selection_id")})
        return CallToolResult(content=[TextContent(type="text", text=text)], structuredContent=result)
    except (SpatialError, ValueError) as exc:
        return CallToolResult(content=[TextContent(type="text", text=str(exc))], isError=True,
                              structuredContent={"error": str(exc)})
    except Exception:
        LOG.exception("Tool failed: %s", name)
        return CallToolResult(content=[TextContent(type="text", text="Internal error; inspect local server logs.")], isError=True,
                              structuredContent={"error": "Internal error; inspect local server logs."})


class PinnedFastMCP(FastMCP):
    async def call_tool(self, name: str, arguments: dict[str, Any]) -> CallToolResult:
        # The SDK otherwise coerces strings such as "true" to booleans and drops
        # unknown arguments. Confirmation and object identity require strict input.
        return await run_in_threadpool(_safe_result, self.spatial_service, name, arguments)


def create_server(project_root: str | Path, *, service=None, token_verifier=None, auth=None,
                  transport_security=None) -> FastMCP:
    service = service or ToolService(project_root)
    server = PinnedFastMCP(
        "Spatial Collab", instructions="Explore the shared pinned spatial transcriptomics project. Always read revision and selection before proposing edits. A proposal is not approval; only apply after explicit researcher confirmation. Changes are overlays. Results are descriptive and may be stale. Do not invent reviewer identity or biological authority.",
        host="127.0.0.1", stateless_http=True, json_response=True,
        token_verifier=token_verifier, auth=auth,
        transport_security=transport_security or TransportSecuritySettings(enable_dns_rebinding_protection=True,
            allowed_hosts=["localhost:*", "127.0.0.1:*", "[::1]:*"],
            allowed_origins=["http://localhost:*", "http://127.0.0.1:*", "http://[::1]:*"]),
    )
    for name in TOOL_NAMES:
        method = getattr(service, name)
        # Route by an opaque project identity, never a model-provided path. The
        # strict dispatcher above still validates the original method arguments.
        def routed(*args, _method=method, **kwargs):
            return _method(*args, **kwargs)
        parameters = list(inspect.signature(method).parameters.values())
        parameters.append(inspect.Parameter("project_id", inspect.Parameter.KEYWORD_ONLY,
                                           default=None, annotation=str | None))
        routed.__signature__ = inspect.signature(method).replace(parameters=parameters)
        routed.__annotations__ = dict(method.__annotations__, project_id=str | None)
        meta = {"ui": {"visibility": ["model", "app"]}}
        if name in {"append_upload", "read_download"}:
            meta["ui"]["visibility"] = ["app"]
        if name == "open_project":
            meta["ui"]["resourceUri"] = UI_URI
        server.add_tool(routed, name=name, description=method.__doc__,
                        annotations=ToolAnnotations(readOnlyHint=name in READ_TOOLS,
                            destructiveHint=name == "discard_transfer", idempotentHint=name in READ_TOOLS, openWorldHint=False),
                        meta=meta, structured_output=False)

    @server.resource(UI_URI, name="spatial_collab_workbench", title="空间转录组协作工作台",
                     mime_type="text/html;profile=mcp-app",
                     meta={"ui": {"prefersBorder": True, "csp": {"connectDomains": [], "resourceDomains": []}}})
    def workbench() -> str:
        return workbench_html()

    server.spatial_service = service
    return server


class LocalOnlyMiddleware(BaseHTTPMiddleware):
    """Reject DNS rebinding and browser cross-origin access, even for reads."""

    async def dispatch(self, request: Request, call_next):
        try:
            host = request.url.hostname
            if host not in {"127.0.0.1", "localhost", "::1"}:
                return JSONResponse({"ok": False, "error": "Loopback Host required"}, status_code=403)
            origin = request.headers.get("origin")
            if origin is not None:
                parsed = urlsplit(origin)
                if parsed.scheme != request.url.scheme or parsed.netloc != request.url.netloc or parsed.path or parsed.query or parsed.fragment:
                    return JSONResponse({"ok": False, "error": "Same-origin request required"}, status_code=403)
            if request.url.path.startswith("/api/"):
                if not origin or not secrets.compare_digest(request.headers.get("x-spatial-csrf", ""), request.app.state.csrf_token):
                    return JSONResponse({"ok": False, "error": "Same-origin workbench token required"}, status_code=403)
        except ValueError:
            return JSONResponse({"ok": False, "error": "Invalid Host or Origin"}, status_code=403)
        response = await call_next(request)
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["Cache-Control"] = "no-store"
        response.headers["Referrer-Policy"] = "no-referrer"
        return response


def create_app(project_root: str | Path):
    """Return a local-only Starlette app with browser API and /mcp transport."""
    server = create_server(project_root)
    app = server.streamable_http_app()
    app.state.mcp = server
    app.state.service = server.spatial_service
    app.state.csrf_token = secrets.token_urlsafe(32)

    async def index(request):
        html = workbench_html(app.state.csrf_token)
        script = workbench_script()
        import base64
        digest = base64.b64encode(hashlib.sha256(script.encode()).digest()).decode()
        return HTMLResponse(html, headers={"Content-Security-Policy":
            f"default-src 'none'; script-src 'sha256-{digest}'; style-src 'unsafe-inline'; connect-src 'self'; img-src data:; base-uri 'none'; form-action 'none'; frame-ancestors 'none'"})

    async def call_tool(request):
        if request.headers.get("content-type", "").split(";")[0].strip() != "application/json":
            return JSONResponse({"ok": False, "error": "JSON body required"}, status_code=415)
        data = bytearray()
        async for chunk in request.stream():
            data.extend(chunk)
            if len(data) > MAX_BODY_BYTES:
                return JSONResponse({"ok": False, "error": "Request too large"}, status_code=413)
        try:
            arguments = json.loads(data)
        except (ValueError, UnicodeDecodeError):
            return JSONResponse({"ok": False, "error": "Invalid JSON"}, status_code=400)
        result = await run_in_threadpool(_safe_result, app.state.service, request.path_params["name"], arguments)
        return JSONResponse({"ok": not result.isError, "result": result.structuredContent}, status_code=400 if result.isError else 200)

    app.routes.extend([Route("/", index), Route("/api/tool/{name}", call_tool, methods=["POST"])])
    app.add_middleware(LocalOnlyMiddleware)
    return app
