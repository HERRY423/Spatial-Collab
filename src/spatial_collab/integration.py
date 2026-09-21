"""Exact integration results and bounded paired RNA/protein reference baselines.

Domain labels are method-specific partitions, never biological truth. External
outputs are validated as data; no user-provided code or shell commands execute.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import os
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.cluster.vq import kmeans2
from scipy.linalg import svd
from scipy.spatial import cKDTree
from threadpoolctl import threadpool_limits

from . import objects
from .proteomics import get_assay
from .store import SpatialError, _text

SCHEMA = "spatial-collab.integration.v1"
BACKEND_VERSION = "balanced-pca-kmeans/1"


def numerical_environment():
    """Package versions alone do not identify MKL/OpenBLAS builds or thread policy."""
    return {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
            "numpy_build": getattr(np.__config__, "CONFIG", {}),
            "scipy_build": getattr(scipy.__config__, "CONFIG", {}),
            "blas_runtime_limit": 1,
            "thread_environment": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}}


def _ids(value, field):
    if not isinstance(value, list) or not value or len(value) > 100_000:
        raise SpatialError(f"{field} must be a nonempty bounded ID list.")
    for item in value:
        _text(item, field, limit=512)
    if len(value) != len(set(value)):
        raise SpatialError(f"Duplicate {field}.")
    return value


def input_identity(project, revision_id, assay_id, rna_features, protein_features, observation_ids):
    context = project.context()
    rev = project.get_revision(revision_id)
    assay = get_assay(project, assay_id)
    return {"project_id": context["project_id"], "source_sha256": context["source_sha256"],
            "fit_revision": revision_id, "fit_revision_sha256": rev["revision_sha256"],
            "protein_assay_id": assay_id, "protein_assay_sha256": assay["assay_sha256"],
            "rna_features": rna_features, "protein_features": protein_features,
            "fit_observation_ids": observation_ids}


def register_result(project, spec):
    """Import full output with exact scope, feature axes, source and revision hashes."""
    if not isinstance(spec, dict) or spec.get("schema") != SCHEMA:
        raise SpatialError("Expected IntegrationResult schema spatial-collab.integration.v1.")
    spec = {k: v for k, v in spec.items() if k not in {"object_id", "object_sha256", "object_kind"}}
    identity = spec.get("input", {})
    required = {"project_id", "source_sha256", "fit_revision", "fit_revision_sha256", "protein_assay_id",
                "protein_assay_sha256", "rna_features", "protein_features", "fit_observation_ids"}
    if set(identity) != required:
        raise SpatialError("Integration input identity fields are missing or unknown.")
    fit_ids = _ids(identity["fit_observation_ids"], "fit_observation_ids")
    output_ids = _ids(spec.get("observation_ids"), "observation_ids")
    for key in ("rna_features", "protein_features"):
        _ids(identity[key], key)
    expected = input_identity(project, identity["fit_revision"], identity["protein_assay_id"],
                              identity["rna_features"], identity["protein_features"], fit_ids)
    if expected != identity:
        raise SpatialError("Integration source, assay or revision version does not match.")
    cells = project.cells(identity["fit_revision"])
    active = {c["cell_id"] for c in cells if c["included"]}
    if not set(fit_ids) <= active:
        raise SpatialError("Fit IDs contain excluded or foreign sample observations.")
    meta = project.summary()["metadata"]
    assay = get_assay(project, identity["protein_assay_id"])
    if not set(identity["rna_features"]) <= set(meta["panel_genes"]) or not set(identity["protein_features"]) <= {f["feature_id"] for f in assay["features"]}:
        raise SpatialError("Integration features are not measured on the declared axes.")
    scope = spec.get("fit_scope")
    if scope not in {"whole_slice", "explicit_roi"}:
        raise SpatialError("Declare whole_slice or explicit_roi fit_scope.")
    if scope == "whole_slice" and set(fit_ids) != active:
        raise SpatialError("Whole-slice fit is missing included observations.")
    mode = spec.get("mode", "fit")
    if mode not in {"fit", "fixed_model_filter"}:
        raise SpatialError("Unknown integration mode.")
    if mode == "fit" and set(output_ids) != set(fit_ids):
        raise SpatialError("Output contains missing or foreign rows relative to the declared fit scope.")
    if mode == "fixed_model_filter":
        _text(spec.get("view_revision"), "view_revision")
        parent = objects.get(project, spec.get("parent_result_id"), "integration")
        if parent["method"]["uses_annotations"] and spec["view_revision"] != parent["input"]["fit_revision"]:
            raise SpatialError("Annotation-dependent models require refitting at a different revision.")
        if parent["input"] != identity or parent["method"] != spec.get("method") or parent["preprocessing"] != spec.get("preprocessing"):
            raise SpatialError("Fixed-model filtering must preserve parent fit identity and method.")
        retained = {c["cell_id"] for c in project.cells(spec.get("view_revision")) if c["included"]}
        if set(output_ids) != set(parent["observation_ids"]) & retained:
            raise SpatialError("Fixed-model output must contain exactly the retained parent rows.")
        positions = {cid: i for i, cid in enumerate(parent["observation_ids"])}
        for key in ("representations", "domains", "uncertainty"):
            if spec.get(key, {}) != {name: [values[positions[cid]] for cid in output_ids] for name, values in parent.get(key, {}).items()}:
                raise SpatialError("Fixed-model outputs cannot change embeddings or labels.")
        if spec.get("uncertainty_meaning", {}) != parent.get("uncertainty_meaning", {}):
            raise SpatialError("Fixed-model filtering cannot reinterpret uncertainty.")
    method = spec.get("method", {})
    for key in ("name", "version", "code_reference"):
        _text(method.get(key), "method." + key)
    if not isinstance(method.get("environment"), dict) or not method["environment"] or not isinstance(method.get("parameters"), dict):
        raise SpatialError("Method environment and parameters must be explicit objects.")
    if type(method.get("seed")) is not int or type(method.get("uses_annotations")) is not bool:
        raise SpatialError("Declare integer seed and whether annotations were used for fitting.")
    if not isinstance(spec.get("preprocessing"), dict) or not spec["preprocessing"]:
        raise SpatialError("Explicit preprocessing required.")
    if spec.get("output_origin") not in {"transformed_measurements", "model_prediction"}:
        raise SpatialError("Embeddings are transformed measurements or model predictions, never raw measurements.")
    reps, domains = spec.get("representations"), spec.get("domains")
    if not isinstance(reps, dict) or not reps or not isinstance(domains, dict) or not domains:
        raise SpatialError("Representations and domain outputs required.")
    if len(reps) > 8 or len(domains) > 8:
        raise SpatialError("At most eight representations/partitions per result.")
    for name, values in reps.items():
        _text(name, "representation name", limit=128)
        array = np.asarray(values, dtype=float)
        if array.ndim != 2 or array.shape[0] != len(output_ids) or not 1 <= array.shape[1] <= 64 or not np.isfinite(array).all():
            raise SpatialError("Representation must be finite N x 1..64, exactly matching output IDs.")
    for name, values in domains.items():
        _text(name, "partition name", limit=128)
        if not isinstance(values, list) or len(values) != len(output_ids):
            raise SpatialError("Every domain partition must match all output IDs.")
        for label in values:
            _text(label, "domain label", limit=128)
    for name, values in spec.get("uncertainty", {}).items():
        a = np.asarray(values, dtype=float)
        if a.shape != (len(output_ids),) or not np.isfinite(a).all():
            raise SpatialError("Uncertainty must match all output rows.")
        _text(spec.get("uncertainty_meaning", {}).get(name), "uncertainty meaning")
    spec["mode"] = mode
    spec["scientific_authorization"] = "NOT_ESTABLISHED"
    return objects.put(project, "integration", spec)


def get_result(project, result_id):
    record = objects.get(project, result_id, "integration")
    head = project.context()["head_revision"]
    view = record.get("view_revision", record["input"]["fit_revision"])
    record["current_revision"] = head
    record["historical_view"] = head != view
    record["annotation_changed_since_fit"] = head != record["input"]["fit_revision"]
    active = {c["cell_id"] for c in project.cells(head) if c["included"]}
    fit_ids = set(record["input"]["fit_observation_ids"])
    record["fit_population_changed"] = fit_ids != active if record["fit_scope"] == "whole_slice" else bool(fit_ids - active)
    record["refit_required_for_current_population"] = record["fit_population_changed"] or record["method"]["uses_annotations"] and record["annotation_changed_since_fit"]
    record["interpretation"] = "Descriptive method partitions; disagreement is a review priority, not an error probability."
    return record


def list_results(project):
    rows = []
    head = project.context()["head_revision"]
    for rid in objects.catalog(project, "integration"):
        item = objects.get(project, rid, "integration")
        rows.append({"result_id": rid, "method": item["method"]["name"], "mode": item["mode"],
                     "observation_count": len(item["observation_ids"]), "fit_revision": item["input"]["fit_revision"],
                     "historical_view": head != item.get("view_revision", item["input"]["fit_revision"]), "partitions": list(item["domains"])})
    return {"results": rows}


def _standardize(matrix):
    means = matrix.mean(axis=0)
    scales = matrix.std(axis=0)
    scales[scales == 0] = 1
    scaled = (matrix - means) / scales
    return scaled / np.sqrt(matrix.shape[1]), {"mean": means.tolist(), "std": scales.tolist(), "block_divisor": float(np.sqrt(matrix.shape[1]))}


@threadpool_limits.wrap(limits=1, user_api="blas")
def run_baseline(project, revision_id, assay_id, *, rna_features=None, protein_features=None,
                 observation_ids=None, components=10, clusters=6, seed=0, protein_transform="log1p",
                 checkpoint=None, backend="balanced_pca", feature_selection="variance"):
    start = time.perf_counter()
    check = checkpoint or (lambda: None)
    if type(seed) is not int or not 0 <= seed < 2**32 or type(components) is not int or not 1 <= components <= 32 or type(clusters) is not int or not 2 <= clusters <= 30:
        raise SpatialError("Use components 1..32, clusters 2..30, and uint32 seed.")
    if protein_transform not in {"log1p", "asinh"} or backend not in {"balanced_pca", "smopca", "spatial_graph"}:
        raise SpatialError("Unsupported preprocessing or backend.")
    if feature_selection not in {"variance", "morans_i"}:
        raise SpatialError("feature_selection must be variance or morans_i.")
    check()
    cells = [c for c in project.cells(revision_id) if c["included"]]
    if observation_ids is not None:
        chosen = set(_ids(observation_ids, "observation_ids"))
        if not chosen <= {c["cell_id"] for c in cells}:
            raise SpatialError("Fit ROI contains unknown/excluded observations.")
        cells = [c for c in cells if c["cell_id"] in chosen]
    if len(cells) < max(clusters + 1, 4):
        raise SpatialError("Insufficient observations for the requested cluster count.")
    assay = get_assay(project, assay_id)
    ids = [c["cell_id"] for c in cells]
    panel = project.summary()["metadata"]["panel_genes"]
    xy_coords = np.array([[c["x"], c["y"]] for c in cells])
    if rna_features is None:
        if feature_selection == "morans_i":
            from .integration_runners import select_spatial_features
            rna_features = select_spatial_features(cells, panel, xy_coords, max_features=128)
        else:
            # Rank measured features by raw count variance; the exact selected axis is recorded.
            sums, squares = Counter(), Counter()
            for c in cells:
                for f, value in c["counts"].items():
                    sums[f] += value
                    squares[f] += value * value
            rna_features = sorted(sums, key=lambda f: (-(squares[f] / len(cells) - (sums[f] / len(cells))**2), f))[:128]
    _ids(rna_features, "rna_features")
    protein_features = protein_features or [f["feature_id"] for f in assay["features"]]
    _ids(protein_features, "protein_features")
    paxis = [f["feature_id"] for f in assay["features"]]
    if not set(rna_features) <= set(panel) or not set(protein_features) <= set(paxis):
        raise SpatialError("Requested feature is unmeasured.")
    if len(rna_features) > 512 or len(cells) * (len(rna_features) + len(protein_features)) > 10_000_000:
        raise SpatialError("Baseline budget: at most 512 RNA features and 10M selected matrix entries.")
    pindex = [paxis.index(f) for f in protein_features]
    protein = np.array([[assay["values"].get(cid, [None] * len(paxis))[j] for j in pindex] for cid in ids], dtype=float)
    if not np.isfinite(protein).all():
        raise SpatialError("Paired baseline requires complete selected protein measurements; no implicit imputation/intersection.")
    rna = np.array([[c["counts"].get(f, 0) for f in rna_features] for c in cells], dtype=float)
    library = np.array([sum(c["counts"].values()) for c in cells])
    if np.any(library == 0):
        raise SpatialError("Zero RNA library: explicitly exclude/review before fitting.")
    rna, rna_scale = _standardize(np.log1p(rna / library[:, None] * 10000))
    if backend == "smopca":
        protein_library = protein.sum(axis=1)
        if np.any(protein_library == 0):
            raise SpatialError("SMOPCA library preprocessing requires nonzero protein library sizes.")
        protein = protein / protein_library[:, None] * 10000
    protein, protein_scale = _standardize(np.log1p(protein) if protein_transform == "log1p" else np.arcsinh(protein / 5))
    representations, domains, models = {}, {}, {}
    if backend == "spatial_graph":
        from .integration_runners import fit_spatial_graph
        representations, domains, models = fit_spatial_graph(
            rna, protein, xy_coords, components, clusters, seed, checkpoint=check, alpha=0.65
        )
    else:
        for name, matrix in (("rna_only", rna), ("protein_only", protein), ("joint", np.concatenate([rna, protein], axis=1))):
            check()
            u, s, vt = svd(matrix, full_matrices=False, check_finite=False)
            n = min(components, len(s))
            embedding = u[:, :n] * s[:n]
            if backend == "smopca" and name == "joint":
                from .smopca_adapter import fit_smopca
                embedding = fit_smopca(rna, protein, xy_coords, n, seed, check)
            check()
            if len(np.unique(embedding, axis=0)) < clusters:
                raise SpatialError(f"{name}: fewer distinct representations than requested clusters.")
            centers, labels = kmeans2(embedding, clusters, iter=50, minit="++", seed=seed, missing="raise")
            representations[name] = embedding.tolist()
            domains[name] = [str(int(x)) for x in labels]
            models[name] = {"loadings": vt[:n].tolist(), "centers": centers.tolist(), "singular_values": s[:n].tolist()} if not (backend == "smopca" and name == "joint") else {"centers": centers.tolist()}
    check()
    method_name = "spatial_graph_kmeans" if backend == "spatial_graph" else ("balanced_pca_kmeans" if backend == "balanced_pca" else "smopca_with_unimodal_controls")
    method = {"name": method_name,
              "version": BACKEND_VERSION, "code_reference": "spatial_collab.integration:" + BACKEND_VERSION,
              "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "environment": numerical_environment(),
              "parameters": {"components": components, "clusters": clusters, "backend": backend, "feature_selection": feature_selection},
              "seed": seed, "uses_annotations": False}
    if backend == "smopca":
        from .smopca_adapter import provenance
        method["environment"]["smopca"] = provenance()
        method["parameters"]["smopca"] = {"kernel": "matern", "nu": 1.5, "gamma_init": 1,
            "estimate_gamma": False, "iterations_gamma": 1, "iterations_sigma_W": 20, "tol_sigma": 2e-5,
            "sigma_init_list": [1, 1], "sigma_xtol_list": [1e-6, 1e-6], "intercept": False, "omics_weight": False,
            "coordinate_scale": "centered XY divided by median nearest-neighbor distance", "joint_scaling": "per-feature zscore; block balancing undone"}
    elif backend == "spatial_graph":
        method["parameters"]["spatial_graph"] = {"smoothing_alpha": 0.65, "neighbors_k": 6, "kernel": "adaptive_gaussian"}
    return register_result(project, {"schema": SCHEMA,
        "input": input_identity(project, revision_id, assay_id, rna_features, protein_features, ids),
        "fit_scope": "explicit_roi" if observation_ids is not None else "whole_slice", "observation_ids": ids,
        "method": method, "preprocessing": {"rna": "log1p(panel_per_10000)", "protein": protein_transform,
            "protein_library_normalization": "per_10000_selected_protein_counts" if backend == "smopca" else "none",
            "protein_asinh_cofactor": 5, "scaling": "feature z-score then divide each block by sqrt(feature count)",
            "rna_feature_selection": feature_selection if rna_features is None else ("explicit" if len(rna_features) != 128 else "see exact recorded axis"),
            "fitted_rna_scale": rna_scale, "fitted_protein_scale": protein_scale, "missing": "reject"},
        "output_origin": "transformed_measurements", "representations": representations, "domains": domains,
        "fitted_models": models, "elapsed_seconds": time.perf_counter() - start,
        "limitations": ["Single-slice descriptive partitions; no biological ground truth or batch correction.",
                        "Spatial coordinates are incorporated via graph smoothing or spatial kernels when using spatial_graph or smopca.",
                        "Equal block scaling is an explicit baseline choice, not learned modality reliability."]})


def fixed_filter(project, result_id, revision_id):
    parent = objects.get(project, result_id, "integration")
    if parent["method"]["uses_annotations"] and revision_id != parent["input"]["fit_revision"]:
        raise SpatialError("Annotation-trained result requires refitting after a revision.")
    retained = {c["cell_id"] for c in project.cells(revision_id) if c["included"]}
    positions = [i for i, cid in enumerate(parent["observation_ids"]) if cid in retained]
    payload = {k: v for k, v in parent.items() if k not in {"object_id", "object_sha256", "object_kind", "elapsed_seconds"}}
    payload.update(mode="fixed_model_filter", parent_result_id=result_id, view_revision=revision_id,
                   observation_ids=[parent["observation_ids"][i] for i in positions])
    for key in ("representations", "domains", "uncertainty"):
        if key in parent:
            payload[key] = {name: [values[i] for i in positions] for name, values in parent[key].items()}
    return register_result(project, payload)


def compare_results(project, left_id, right_id=None, left_partition="rna_only", right_partition="joint", limit=100, offset=0):
    if type(limit) is not int or not 1 <= limit <= 1000 or type(offset) is not int or offset < 0:
        raise SpatialError("Use limit 1..1000 and nonnegative offset.")
    left, right = get_result(project, left_id), get_result(project, right_id or left_id)
    if left_partition not in left["domains"] or right_partition not in right["domains"]:
        raise SpatialError("Unknown method partition.")
    common = sorted(set(left["observation_ids"]) & set(right["observation_ids"]))
    if len(common) < 2:
        raise SpatialError("At least two common observations required.")
    li, ri = ({cid: i for i, cid in enumerate(r["observation_ids"])} for r in (left, right))
    a = np.array([left["domains"][left_partition][li[c]] for c in common])
    b = np.array([right["domains"][right_partition][ri[c]] for c in common])
    cells = {c["cell_id"]: c for c in project.cells()}
    xy = np.array([[cells[c]["x"], cells[c]["y"]] for c in common])
    _, neighbors = cKDTree(xy).query(xy, k=min(7, len(common)))
    # Remove self by ID (not by distance/order, which is ambiguous at identical XY).
    scores = []
    for i, near in enumerate(neighbors):
        near = [int(j) for j in np.atleast_1d(near) if j != i][:6]
        scores.append(float(np.mean((a[i] == a[near]) != (b[i] == b[near]))))
    pairs = Counter(zip(a.tolist(), b.tolist()))
    def choose2(n):
        return n * (n - 1) / 2
    same = sum(choose2(n) for n in pairs.values())
    ca, cb = sum(choose2(n) for n in Counter(a).values()), sum(choose2(n) for n in Counter(b).values())
    expected = ca * cb / choose2(len(common))
    denominator = (ca + cb) / 2 - expected
    ari = (same - expected) / denominator if denominator else 1.0
    order = sorted(range(len(common)), key=lambda i: (-scores[i], common[i]))
    rows = [{"cell_id": common[i], "x": float(xy[i, 0]), "y": float(xy[i, 1]),
             "left_domain": str(a[i]), "right_domain": str(b[i]), "boundary_disagreement": scores[i],
             "current_label": cells[common[i]]["label"], "included": cells[common[i]]["included"]} for i in order[offset:offset + limit]]
    return {"left_id": left_id, "right_id": right_id or left_id, "common_count": len(common),
            "left_only_count": len(li) - len(common), "right_only_count": len(ri) - len(common),
            "adjusted_rand_index": float(ari), "mean_boundary_disagreement": float(np.mean(scores)),
            "rows": rows, "offset": offset, "next_offset": offset + limit if offset + limit < len(common) else None,
            "label_comparison": "Permutation-invariant co-membership on six nearest common observations; numeric labels are not aligned identities.",
            "scope": "Common output observations only; removed observations are counted separately.",
            "scientific_authorization": "NOT_ESTABLISHED"}


def inspect_result(project, result_id, partition="joint", bounds=None, limit=2000, offset=0):
    result = get_result(project, result_id)
    if partition not in result["domains"] or type(limit) is not int or not 1 <= limit <= 10000 or type(offset) is not int or offset < 0:
        raise SpatialError("Unknown partition or invalid pagination.")
    if bounds is not None and (not isinstance(bounds, list) or len(bounds) != 4 or not np.isfinite(bounds).all() or bounds[0] > bounds[2] or bounds[1] > bounds[3]):
        raise SpatialError("Bounds require finite xmin,ymin,xmax,ymax.")
    cells = {c["cell_id"]: c for c in project.cells()}
    rows = []
    for i, cid in enumerate(result["observation_ids"]):
        c = cells[cid]
        if bounds is None or bounds[0] <= c["x"] <= bounds[2] and bounds[1] <= c["y"] <= bounds[3]:
            rows.append({"cell_id": cid, "x": c["x"], "y": c["y"], "domain": result["domains"][partition][i],
                         "included": c["included"], "label": c["label"]})
    return {"result_id": result_id, "partition": partition, "points": rows[offset:offset + limit],
            "total_in_view": len(rows), "next_offset": offset + limit if offset + limit < len(rows) else None,
            "historical_view": result["historical_view"], "mode": result["mode"], "input": result["input"],
            "method": result["method"], "preprocessing": result["preprocessing"], "display_sampling": False}
