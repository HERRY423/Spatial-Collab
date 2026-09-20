"""One pinned project, exposed through MCP and the local review workbench."""

from __future__ import annotations

import hashlib
import json
import logging
import math
import secrets
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

from . import analysis
from .store import Project, SpatialError

UI_URI = "ui://spatial-collab/workbench.html"
MAX_VIEW_CELLS = 10_000
MAX_BODY_BYTES = 2_000_000
STATIC = Path(__file__).with_name("static")
LOG = logging.getLogger(__name__)
READ_TOOLS = {"open_project", "get_selection", "inspect_selection", "get_run", "get_revision",
              "get_context", "query_observations", "get_capabilities", "get_proposal", "inspect_quality",
              "get_feature_catalog", "list_assets", "inspect_asset", "get_hypothesis_plan", "list_hypothesis_plans"}
READ_TOOLS.update({"get_overview", "list_assays", "inspect_protein"})
TOOL_NAMES = (
    "open_project", "get_selection", "set_selection", "inspect_selection",
    "propose_revision", "apply_revision", "revert_revision", "run_comparison",
    "get_run", "get_revision", "export_review_bundle", "get_context", "query_observations",
    "get_capabilities", "run_region_comparison", "get_proposal", "inspect_quality", "run_sensitivity",
    "get_feature_catalog", "list_assets", "inspect_asset", "create_hypothesis_plan", "revise_hypothesis_plan",
    "freeze_hypothesis_plan", "get_hypothesis_plan", "list_hypothesis_plans", "run_hypothesis_plan",
    "get_overview", "list_assays", "inspect_protein", "run_protein_comparison",
)


class ToolService:
    """Shared dispatch. No tool accepts an input/output filesystem path."""

    def __init__(self, project_root: str | Path):
        self.project = Project(Path(project_root).resolve())
        self._calls = {
            name: validate_call(config=ConfigDict(strict=True))(getattr(self, name))
            for name in TOOL_NAMES
        }

    def call(self, name: str, arguments: dict | None = None) -> dict:
        if name not in self._calls:
            raise SpatialError("Unknown tool")
        if arguments is not None and not isinstance(arguments, dict):
            raise SpatialError("Tool arguments must be an object")
        try:
            return self._calls[name](**(arguments or {}))
        except ValidationError as exc:
            issues = [f"{'.'.join(str(x) for x in e['loc'])}: {e['msg']}" for e in exc.errors()]
            raise SpatialError("Invalid arguments: " + "; ".join(issues)) from None

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

    def get_selection(self, selection_id: str | None = None) -> dict:
        """Read exact server-persisted cell IDs, geometry, revision and stale status; optional historical selection ID."""
        summary = self.project.summary()
        selection = self.project.get_selection(selection_id) if selection_id else summary["selection"]
        return {"selection": selection, "head_revision": summary["head_revision"]}

    def get_overview(self) -> dict:
        """Overview of RNA and paired protein layers, shared selection, revision and recent saved analyses."""
        from .proteomics import list_assays
        summary = self.project.summary()
        return {"project_id": summary["project_id"], "name": summary["metadata"]["name"],
                "sample_id": summary["metadata"].get("sample_id"), "source_sha256": summary["source_sha256"],
                "head_revision": summary["head_revision"], "observation_count": summary["cell_count"],
                "observation_unit": summary["metadata"].get("observation_unit", "cell"),
                "units": summary["metadata"]["units"], "rna_feature_count": len(summary["metadata"]["panel_genes"]),
                "protein_assays": list_assays(self.project)["assays"], "selection_count": (summary["selection"] or {}).get("cell_count", 0),
                "revision_count": summary["revision_count"], "recent_runs": summary["runs"][-12:],
                "scope": "One paired spatial sample; exact shared objects, no cross-sample integration inferred.",
                "scientific_authorization": "NOT_ESTABLISHED"}

    def list_assays(self) -> dict:
        """List registered protein assays and measured features without exporting their full matrices."""
        from .proteomics import list_assays
        return list_assays(self.project)

    def inspect_protein(self, assay_id: str, revision_id: str, feature: str, selection_id: str | None = None,
                        scale: str = "raw", cofactor: float = 5.0) -> dict:
        """Read an exact protein channel at a shared revision/ROI. Missing and zero differ; raw/log1p/asinh are explicit display transforms."""
        from .proteomics import inspect_protein
        return inspect_protein(self.project, assay_id, revision_id, feature, selection_id, scale, cofactor)

    def run_protein_comparison(self, assay_id: str, base_revision: str, target_revision: str, selection_id: str,
                               features: list[str], background_selection_id: str | None = None,
                               scale: str = "raw", cofactor: float = 5.0, rna_gene: str | None = None) -> dict:
        """Compare frozen foreground/background protein signals before/after shared annotation/inclusion revisions; optional exact paired RNA Spearman concordance. Descriptive only."""
        from .proteomics import compare_protein_regions
        return compare_protein_regions(self.project, assay_id, base_revision, target_revision, selection_id,
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
                "limits": {"single_slice": True, "centroids_only": True, "scientific_authorization": "NOT_ESTABLISHED",
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


def workbench_html(csrf_token: str = "") -> str:
    html = (STATIC / "workbench.html").read_text(encoding="utf-8")
    return html.replace("/*__STYLE__*/", (STATIC / "workbench.css").read_text(encoding="utf-8")).replace(
        "/*__SCRIPT__*/", (STATIC / "workbench.js").read_text(encoding="utf-8")).replace("__CSRF_TOKEN__", csrf_token)


def _safe_result(service: ToolService, name: str, arguments: dict) -> CallToolResult:
    try:
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


def create_server(project_root: str | Path) -> FastMCP:
    service = ToolService(project_root)
    server = PinnedFastMCP(
        "Spatial Collab", instructions="Explore the shared pinned spatial transcriptomics project. Always read revision and selection before proposing edits. A proposal is not approval; only apply after explicit researcher confirmation. Changes are overlays. Results are descriptive and may be stale. Do not invent reviewer identity or biological authority.",
        host="127.0.0.1", stateless_http=True, json_response=True,
        transport_security=TransportSecuritySettings(enable_dns_rebinding_protection=True,
            allowed_hosts=["localhost:*", "127.0.0.1:*", "[::1]:*"],
            allowed_origins=["http://localhost:*", "http://127.0.0.1:*", "http://[::1]:*"]),
    )
    for name in TOOL_NAMES:
        method = getattr(service, name)

        meta = {"ui": {"visibility": ["model", "app"]}}
        if name == "open_project":
            meta["ui"]["resourceUri"] = UI_URI
        server.add_tool(method, name=name, description=method.__doc__,
                        annotations=ToolAnnotations(readOnlyHint=name in READ_TOOLS,
                            destructiveHint=False, idempotentHint=name in READ_TOOLS, openWorldHint=False),
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
        script = (STATIC / "workbench.js").read_text(encoding="utf-8")
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
