"""Executable official method adapters; scientific outputs stay distinct by task."""

from importlib import metadata, util
import hashlib
import random

import numpy as np
from scipy import sparse
from scipy.linalg import svd
from scipy.optimize import nnls
from scipy.spatial import cKDTree
from threadpoolctl import threadpool_limits

from . import objects
from .integration import numerical_environment
from .spatial_statistics import graph, normalized_counts, svg, bh
from .store import SpatialError, _hash
from .workflow_inputs import load_counts

METHODS = {
    "progeny": (
        "PROGENy / ULM 通路活性",
        "pathway_activity",
        "decoupler",
        "decoupler",
        "https://decoupler.readthedocs.io/en/stable/api/generated/decoupler.mt.ulm.html",
    ),
    "mofa": ("MOFA+ 多模态因子", "multimodal", "mofapy2", "mofapy2", "https://biofam.github.io/MOFA2/"),
    "mefisto": (
        "MEFISTO 空间多模态因子",
        "multimodal",
        "mofapy2",
        "mofapy2",
        "https://biofam.github.io/MOFA2/MEFISTO.html",
    ),
    "nnls": (
        "参考 NNLS 去卷积基线",
        "deconvolution",
        "scipy",
        "scipy",
        "https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.nnls.html",
    ),
    "cell2location": (
        "Cell2location 丰度后验",
        "deconvolution",
        "cell2location",
        "cell2location",
        "https://cell2location.readthedocs.io/",
    ),
    "harmony": (
        "Harmony 表达批次校正",
        "batch_correction",
        "harmonypy",
        "harmonypy",
        "https://github.com/slowkow/harmonypy",
    ),
    "paste": ("PASTE 切片配准", "alignment", "paste", "paste-bio", "https://paste-bio.readthedocs.io/"),
    "spagcn": ("SpaGCN 空间结构域", "domains", "SpaGCN", "SpaGCN", "https://github.com/jianhuupenn/SpaGCN"),
    "spatial_spectral": (
        "稀疏空间图谱聚类基线",
        "domains",
        "sklearn",
        "scikit-learn",
        "https://scikit-learn.org/stable/modules/generated/sklearn.cluster.SpectralClustering.html",
    ),
    "moran_svg": (
        "Moran's I 置换与 FDR",
        "svg",
        "scipy",
        "scipy",
        "https://pysal.org/esda/generated/esda.Moran.html",
    ),
    "spatial_lr": (
        "数据库支持的空间配体—受体检验",
        "communication",
        "liana",
        "liana",
        "https://liana.readthedocs.io/en/stable/tutorials/notebooks/prior_knowledge.html",
    ),
}


def catalog():
    rows = []
    for key, (name, task, module, dist, source) in METHODS.items():
        available = util.find_spec(module) is not None
        try:
            version = metadata.version(dist)
        except metadata.PackageNotFoundError:
            version = None
        rows.append(
            {
                "method": key,
                "name": name,
                "task": task,
                "available": available,
                "package": dist,
                "version": version,
                "source": source,
                "input_roles": (
                    ["spatial", "reference"]
                    if key in {"nnls", "cell2location"}
                    else ["spatial", "spatial"]
                    if key in {"harmony", "paste"}
                    else ["spatial"]
                ),
                "parameter_contract": METHOD_PARAMETERS[key],
            }
        )
    task_order = [
        "multimodal",
        "deconvolution",
        "batch_correction",
        "alignment",
        "domains",
        "svg",
        "communication",
        "pathway_activity",
    ]
    rows.sort(key=lambda row: task_order.index(row["task"]))
    return {
        "methods": rows,
        "availability_meaning": "Package discoverability only; successful runtime receipt is separate.",
    }


METHOD_PARAMETERS = {
    "mofa": {
        "assay_id": "required exact paired protein assay",
        "n_features": "2..2000, default 512",
        "features": "optional explicit measured stable IDs",
        "components": "2..32, bounded by inputs, default 8",
        "iterations": "10..5000, default 200",
    },
    "mefisto": {
        "assay_id": "required exact paired protein assay",
        "n_features": "2..2000, default 512",
        "features": "optional explicit measured stable IDs",
        "components": "2..32, bounded by inputs, default 8",
        "iterations": "10..5000, default 200",
        "frac_inducing": "0.02..0.5, default 0.1; sparse GP",
    },
    "nnls": {
        "n_features": "2..5000, default 512",
        "features": "optional explicit shared stable IDs",
        "min_reference_cells": "2..1000, default 5, per declared type",
    },
    "cell2location": {
        "n_features": "2..5000, default 512",
        "features": "optional explicit shared stable IDs",
        "min_reference_cells": "2..1000, default 5",
        "cells_per_location": "REQUIRED 0.1..1000; histology-supported prior",
        "detection_alpha": "0.01..10000, default 20",
        "max_epochs": "10..30000, default 2000",
        "posterior_samples": "20..1000, default 100",
    },
    "harmony": {
        "n_features": "2..2000, default 512",
        "features": "optional explicit shared stable IDs",
        "components": "2..50, default 10",
        "iterations": "2..100, default 20",
        "sample_conditions": "required nonempty labels in input order",
        "sample_batches": "required batch labels in input order; reject complete confounding",
    },
    "paste": {
        "n_features": "2..2000, default 512",
        "features": "optional explicit shared stable IDs",
        "overlap_assumption": "REQUIRED full_overlap after research review",
        "alpha": "0..1, default 0.1",
        "iterations": "5..1000, default 100",
    },
    "spagcn": {
        "n_features": "2..2000, default 512",
        "features": "optional explicit measured stable IDs",
        "components": "2..50, default 10",
        "clusters": "2..30, default 6",
        "neighbors": "1..64, default 6",
        "length_scale": "optional positive spatial Gaussian length scale",
        "max_epochs": "10..2000, default 200",
    },
    "spatial_spectral": {
        "n_features": "2..2000, default 512",
        "features": "optional explicit measured stable IDs",
        "components": "2..50, default 10",
        "clusters": "2..30, default 6",
        "neighbors": "1..64, default 6",
    },
    "moran_svg": {
        "features": "optional explicit testing family; default all measured genes",
        "neighbors": "1..64, default 6",
        "permutations": "19..9999, default 99",
        "min_detected": "2..N, default 5",
    },
    "spatial_lr": {
        "pairs": "optional unique [ligand,receptor] pairs from species-specific LIANA consensus; default whole resource",
        "neighbors": "1..64, default 6",
        "permutations": "19..9999, default 99",
    },
    "progeny": {
        "top_targets": "10..1000, default 100",
        "min_targets": "3..100, default 5",
        "license": "academic (default), commercial, or nonprofit",
    },
}


def _integer(params, key, default, low, high):
    value = params.get(key, default)
    if type(value) is not int or not low <= value <= high:
        raise SpatialError(f"{key} must be an integer in {low}..{high}.")
    return value


def _number(params, key, default, low, high):
    value = params.get(key, default)
    if type(value) not in (int, float) or not np.isfinite(value) or not low <= value <= high:
        raise SpatialError(f"{key} must be finite in {low}..{high}.")
    return float(value)


def _features(records, matrices, params, maximum=2000):
    common = set(records[0]["features"])
    for r in records[1:]:
        common &= set(r["features"])
    requested = params.get("features")
    if requested is not None:
        if (
            not isinstance(requested, list)
            or not requested
            or len(requested) != len(set(requested))
            or not set(requested) <= common
        ):
            raise SpatialError("Explicit features must be unique and measured in every declared input.")
        chosen = requested
    else:
        limit = _integer(params, "n_features", 512, 2, maximum)
        # Normalized dispersion is a baseline; SVG is a separate tested workflow.
        y = normalized_counts(matrices[0])
        means = np.asarray(y.mean(axis=0)).ravel()
        variance = np.asarray(y.multiply(y).mean(axis=0)).ravel() - means**2
        scores = {g: variance[i] for i, g in enumerate(records[0]["features"]) if g in common}
        chosen = sorted(scores, key=lambda g: (-scores[g], g))[:limit]
    if not 2 <= len(chosen) <= maximum:
        raise SpatialError(f"Need 2..{maximum} shared measured features.")
    blocks = []
    for r, matrix in zip(records, matrices):
        lookup = {g: i for i, g in enumerate(r["features"])}
        blocks.append(matrix[:, [lookup[g] for g in chosen]])
    return chosen, blocks


def _dense(matrix, max_entries=20_000_000):
    if matrix.shape[0] * matrix.shape[1] > max_entries:
        raise SpatialError("Dense method budget exceeded; explicitly restrict features or fit ROI.")
    return matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)


def _zscore(x):
    std = x.std(axis=0)
    std[std < 1e-12] = 1
    return (x - x.mean(axis=0)) / std


def _pca(x, components):
    u, s, _ = svd(_zscore(_dense(x)), full_matrices=False)
    return u[:, :components] * s[:components]


def _anndata(x, record, features):
    import anndata as ad
    import pandas as pd

    return ad.AnnData(
        x.copy(), obs=pd.DataFrame(index=record["observation_ids"]), var=pd.DataFrame(index=features)
    )


def _reference(records, matrices, params):
    if len(records) != 2 or records[0]["kind"] != "spatial" or records[1]["kind"] != "reference":
        raise SpatialError(
            "Deconvolution needs spatial input followed by an annotated single-cell reference."
        )
    features, blocks = _features(records, matrices, params, maximum=5000)
    labels = np.asarray(records[1]["labels"])
    types = sorted(set(labels))
    minimum = _integer(params, "min_reference_cells", 5, 2, 1000)
    sizes = {t: int(np.sum(labels == t)) for t in types}
    if min(sizes.values()) < minimum:
        raise SpatialError(
            "Each reference type needs the declared minimum number of cells; no silent type dropping."
        )
    signature = np.column_stack([np.asarray(blocks[1][labels == t].mean(axis=0)).ravel() for t in types])
    if np.any(signature.sum(axis=0) <= 0) or np.linalg.matrix_rank(signature) < len(types):
        raise SpatialError(
            "Reference signatures are zero or rank deficient on shared features; types are not identifiable."
        )
    return features, blocks[0], types, signature, sizes


def _multimodal(project, record, raw, params, method, seed, checkpoint):
    from .proteomics import get_assay
    from mofapy2.run.entry_point import entry_point

    assay = get_assay(project, params.get("assay_id"))
    if (
        record["provenance"].get("source_sha256") != assay["source_sha256"]
        or record["provenance"].get("project_id") != assay["project_id"]
    ):
        raise SpatialError("Multimodal inputs must share the exact project/source with the protein assay.")
    genes, _ = _features([record], [raw], params)
    lookup = {g: i for i, g in enumerate(record["features"])}
    rna = _dense(normalized_counts(raw)[:, [lookup[g] for g in genes]])
    proteins = [f["feature_id"] for f in assay["features"]]
    protein = np.asarray(
        [assay["values"].get(cid, [None] * len(proteins)) for cid in record["observation_ids"]], dtype=float
    )
    if not np.isfinite(protein).all():
        raise SpatialError(
            "This paired MOFA adapter requires complete protein measurements; no hidden imputation."
        )
    if len(rna) * (len(genes) + len(proteins)) > 20_000_000:
        raise SpatialError("Multimodal matrix budget exceeded.")
    factors = _integer(params, "components", 8, 2, min(32, len(rna) - 1, len(genes), len(proteins)))
    iterations = _integer(params, "iterations", 200, 10, 5000)
    model = entry_point()
    model.set_data_options(scale_views=True)
    model.set_data_matrix(
        [[rna], [np.log1p(protein)]],
        likelihoods=["gaussian", "gaussian"],
        views_names=["RNA", "protein"],
        groups_names=[record["sample_id"]],
        samples_names=[record["observation_ids"]],
        features_names=[genes, proteins],
    )
    model.set_model_options(factors=factors, spikeslab_weights=True, ard_factors=True, ard_weights=True)
    model.set_train_options(iter=iterations, seed=seed, verbose=False, quiet=True, convergence_mode="medium")
    if method == "mefisto":
        xy = np.asarray(record["coordinates"], dtype=float)
        if len(xy) > 10000:
            raise SpatialError("MEFISTO local budget is 10,000 observations; choose an explicit ROI.")
        model.set_covariates([_zscore(xy)], covariates_names=["x", "y"])
        model.set_smooth_options(
            sparseGP=True,
            frac_inducing=_number(params, "frac_inducing", 0.1, 0.02, 0.5),
            start_opt=min(10, iterations // 2),
            opt_freq=10,
        )
    checkpoint()
    model.build()
    model.run()
    checkpoint()
    z = np.asarray(model.model.nodes["Z"].getExpectation())
    weights = [np.asarray(node) for node in model.model.nodes["W"].getExpectation()]
    return {
        "embedding": z.tolist(),
        "features": genes,
        "protein_features": proteins,
        "factor_loadings": {"RNA": weights[0].tolist(), "protein": weights[1].tolist()},
        "controls": {
            "rna_pca": _pca(rna, factors).tolist(),
            "protein_pca": _pca(np.log1p(protein), factors).tolist(),
        },
        "iterations_requested": iterations,
        "training_diagnostics": _clean(model.model.getTrainingStats()),
        "spatial_covariates_used": method == "mefisto",
        "preprocessing": {
            "rna": "full-panel library 10000 then log1p; selected measured features",
            "protein": "log1p measured nonnegative values",
            "views": "official MOFA scale_views=True",
            "coordinates": "per-axis z score for MEFISTO; original coordinates retained"
            if method == "mefisto"
            else "unused",
        },
        "protein_assay_sha256": assay["assay_sha256"],
        "interpretation": "Latent factors/loadings; not cell identities. MOFA does not use spatial coordinates; MEFISTO does.",
    }


def _clean(value):
    if type(value).__module__.startswith("pandas.") and hasattr(value, "to_dict"):
        return _clean(value.to_dict(orient="split") if value.ndim == 2 else value.to_dict())
    if isinstance(value, np.str_):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (tuple, list, np.ndarray)):
        return [_clean(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def _deconvolve(records, matrices, params, method, seed, checkpoint):
    genes, x, types, signature, sizes = _reference(records, matrices, params)
    if method == "nnls":
        signature = signature / signature.sum(axis=0)
        y = _dense(x)
        if np.any(y.sum(axis=1) <= 0):
            raise SpatialError(
                "Some spatial observations have no counts on shared selected features; revise the explicit feature set."
            )
        y = y / y.sum(axis=1)[:, None]
        coef, residual = [], []
        for i, row in enumerate(y):
            if i % 128 == 0:
                checkpoint()
            beta, error = nnls(signature, row)
            coef.append((beta / beta.sum()).tolist() if beta.sum() else [0] * len(types))
            residual.append(float(error))
        return {
            "cell_types": types,
            "features": genes,
            "rna_contribution_fraction": coef,
            "residual_l2": residual,
            "signature_condition_number": float(np.linalg.cond(signature)),
            "reference_cells_per_type": sizes,
            "interpretation": "Nonnegative reference regression. RNA contributions are not absolute cell counts or calibrated cell fractions; absent reference types cannot be discovered.",
        }
    import pandas as pd
    import scvi
    import torch
    from cell2location.models import Cell2location

    torch.set_num_threads(1)
    scvi.settings.seed = seed
    expected_cells = _number(params, "cells_per_location", None, 0.1, 1000)
    alpha = _number(params, "detection_alpha", 20, 0.01, 10000)
    epochs = _integer(params, "max_epochs", 2000, 10, 30000)
    draws = _integer(params, "posterior_samples", 100, 20, 1000)
    data = _anndata(x, records[0], genes)
    Cell2location.setup_anndata(data)
    model = Cell2location(
        data,
        cell_state_df=pd.DataFrame(signature, index=genes, columns=types),
        N_cells_per_location=expected_cells,
        detection_alpha=alpha,
    )
    checkpoint()
    model.train(
        max_epochs=epochs,
        batch_size=None,
        train_size=1,
        accelerator="cpu",
        device=1,
        enable_progress_bar=False,
        logger=False,
    )
    checkpoint()
    data = model.export_posterior(
        data,
        sample_kwargs={
            "num_samples": draws,
            "batch_size": min(1024, len(data)),
            "accelerator": "cpu",
            "device": 1,
        },
    )
    outputs = {q: np.asarray(data.obsm[f"{q}_cell_abundance_w_sf"]).tolist() for q in ("means", "q05", "q95")}
    return {
        "cell_types": types,
        "features": genes,
        "abundance": outputs,
        "reference_cells_per_type": sizes,
        "signature_estimator": "Per-type raw count means; uncertainty in reference signatures is not fitted.",
        "training_history": {k: _clean(np.asarray(v)) for k, v in model.history.items()},
        "max_epochs": epochs,
        "posterior_samples": draws,
        "interpretation": "Model-dependent absolute abundance posterior with an explicitly supplied cells/location prior; training length alone does not establish convergence.",
    }


def _harmony(records, matrices, params, seed):
    import harmonypy
    import pandas as pd

    if len(records) < 2:
        raise SpatialError("Harmony requires two or more declared samples.")
    if len({r["sample_id"] for r in records}) != len(records):
        raise SpatialError("Distinct sample IDs required for batch correction.")
    conditions = params.get("sample_conditions")
    batches = params.get("sample_batches")
    if (
        not isinstance(conditions, list)
        or not isinstance(batches, list)
        or len(conditions) != len(records)
        or len(batches) != len(records)
    ):
        raise SpatialError("Declare sample_conditions and sample_batches for every input, in input order.")
    if any(not isinstance(v, str) or not v.strip() for v in conditions + batches) or len(set(batches)) < 2:
        raise SpatialError("Nonempty design labels and at least two batches required.")
    if len(set(conditions)) > 1 and all(
        len({c for c, b in zip(conditions, batches) if b == batch}) == 1 for batch in set(batches)
    ):
        raise SpatialError(
            "Condition is completely confounded with batch; correction cannot identify condition versus batch effects."
        )
    genes, _ = _features(records, matrices, params)
    normalized = []
    batch_rows = []
    for r, x, batch in zip(records, matrices, batches):
        lookup = {g: i for i, g in enumerate(r["features"])}
        normalized.append(normalized_counts(x)[:, [lookup[g] for g in genes]])
        batch_rows.extend([batch] * x.shape[0])
    components = _integer(params, "components", 10, 2, min(50, len(genes)))
    baseline = _pca(sparse.vstack(normalized), components)
    result = harmonypy.run_harmony(
        baseline,
        pd.DataFrame({"batch": batch_rows}),
        "batch",
        random_state=seed,
        nclust=min(30, max(2, len(baseline) // 30)),
        max_iter_harmony=_integer(params, "iterations", 20, 2, 100),
        verbose=False,
    )
    corrected = np.asarray(result.Z_corr)
    if corrected.shape != baseline.shape and corrected.T.shape == baseline.shape:
        corrected = corrected.T
    if corrected.shape != baseline.shape:
        raise SpatialError("Harmony output axes do not match the input.")

    def mixing(embedding):
        idx = cKDTree(embedding).query(embedding, k=min(16, len(embedding)))[1][:, 1:]
        b = np.asarray(batch_rows)
        return float(np.mean(b[idx] != b[:, None]))

    return {
        "features": genes,
        "embedding": corrected.tolist(),
        "uncorrected_embedding": baseline.tolist(),
        "sample_conditions": conditions,
        "sample_batches": batches,
        "cross_batch_neighbor_fraction": {"before": mixing(baseline), "after": mixing(corrected)},
        "interpretation": "Corrected expression embedding only; no coordinate registration or corrected raw counts. Increased batch mixing alone is not evidence of biological preservation.",
    }


def _paste(records, matrices, params, seed):
    import paste

    if len(records) != 2 or any(r["kind"] != "spatial" for r in records):
        raise SpatialError("PASTE requires two spatial slices.")
    if matrices[0].shape[0] * matrices[1].shape[0] > 4_000_000:
        raise SpatialError("PASTE transport budget is 4M pairs; choose explicit comparable ROIs.")
    units = [r["provenance"].get("units", "unknown") for r in records]
    if units[0] != units[1]:
        raise SpatialError(
            "Cross-slice coordinate units differ; explicitly convert and register before alignment."
        )
    if params.get("overlap_assumption") != "full_overlap":
        raise SpatialError(
            "This PASTE adapter requires an explicit full_overlap assumption; use a partial-overlap method otherwise."
        )
    genes, _ = _features(records, matrices, params)
    slides = []
    for r, x in zip(records, matrices):
        lookup = {g: i for i, g in enumerate(r["features"])}
        y = normalized_counts(x)[:, [lookup[g] for g in genes]]
        data = _anndata(y, r, genes)
        data.obsm["spatial"] = np.asarray(r["coordinates"])
        slides.append(data)
    from .paste_compat import line_search_compatibility

    with line_search_compatibility() as compatibility:
        pi = np.asarray(
            paste.pairwise_align(
                slides[0],
                slides[1],
                alpha=_number(params, "alpha", 0.1, 0, 1),
                numItermax=_integer(params, "iterations", 100, 5, 1000),
                dissimilarity="kl",
                norm=True,
                verbose=False,
            )
        )
    if pi.shape != (len(slides[0]), len(slides[1])) or not np.isfinite(pi).all() or np.any(pi < -1e-10):
        raise SpatialError("PASTE returned invalid transport.")
    # Preserve the complete sparse transport, not a misleading one-to-one assignment.
    transport = sparse.coo_matrix(pi)
    rowmass = pi.sum(axis=1)
    conditional = pi / np.maximum(rowmass[:, None], 1e-300)
    mapped = conditional @ slides[1].obsm["spatial"]
    return {
        "features": genes,
        "transport": {
            "shape": list(pi.shape),
            "row": transport.row.tolist(),
            "col": transport.col.tolist(),
            "mass": transport.data.tolist(),
        },
        "compatibility": compatibility,
        "target_barycentric_coordinates": mapped.tolist(),
        "row_mass": rowmass.tolist(),
        "row_entropy": (-np.sum(conditional * np.log(np.maximum(conditional, 1e-300)), axis=1)).tolist(),
        "marginal_max_error": float(
            max(np.abs(rowmass - 1 / len(slides[0])).max(), np.abs(pi.sum(axis=0) - 1 / len(slides[1])).max())
        ),
        "interpretation": "Optimal-transport computed correspondence under full-overlap assumption; barycentric coordinates are inferred, never measured coordinates or exact same-cell identity.",
    }


def _domains(record, raw, params, method, seed, checkpoint):
    genes, _ = _features([record], [raw], params)
    lookup = {g: i for i, g in enumerate(record["features"])}
    y = normalized_counts(raw)[:, [lookup[g] for g in genes]]
    clusters = _integer(params, "clusters", 6, 2, min(30, raw.shape[0] - 1))
    components = _integer(params, "components", 10, 2, min(50, len(genes), raw.shape[0] - 1))
    xy = np.asarray(record["coordinates"])
    z = _pca(y, components)
    k = _integer(params, "neighbors", 6, 1, min(64, len(xy) - 1))
    w = graph(xy, k)
    if method == "spatial_spectral":
        from sklearn.cluster import SpectralClustering

        # Expression similarity only on real spatial edges. No dense NxN affinity.
        edges = w.tocoo()
        dist2 = np.sum((z[edges.row] - z[edges.col]) ** 2, axis=1)
        scale = np.median(dist2[dist2 > 0]) if np.any(dist2 > 0) else 1
        affinity = sparse.csr_matrix(
            (edges.data * np.exp(-dist2 / (2 * scale)), (edges.row, edges.col)), shape=w.shape
        )
        labels = SpectralClustering(
            n_clusters=clusters, affinity="precomputed", random_state=seed, assign_labels="cluster_qr"
        ).fit_predict(affinity)
        details = {"affinity_nonzero": affinity.nnz, "expression_kernel_scale": float(scale)}
    else:
        import SpaGCN as spg
        import torch
        from scipy.spatial.distance import cdist

        if len(xy) > 6000:
            raise SpatialError(
                "Official SpaGCN dense-adjacency budget is 6000 positions; use the sparse baseline or explicit ROI."
            )
        torch.set_num_threads(1)
        torch.manual_seed(seed)
        random.seed(seed)
        np.random.seed(seed)
        # SpaGCN 1.2.7 calls the removed scipy .A attribute on sparse X.
        # Its official dense input path is supported and stays within our matrix budget.
        data = _anndata(_dense(y), record, genes)
        model = spg.SpaGCN()
        length = _number(
            params, "length_scale", float(np.median(cKDTree(xy).query(xy, k=2)[0][:, 1])), 1e-9, 1e9
        )
        model.set_l(length)
        checkpoint()
        model.train(
            data,
            cdist(xy, xy).astype(np.float32),
            init_spa=False,
            init="kmeans",
            n_clusters=clusters,
            num_pcs=components,
            max_epochs=_integer(params, "max_epochs", 200, 10, 2000),
            lr=0.05,
            tol=1e-3,
        )
        labels, probabilities = model.predict()
        learned, _ = model.model.predict(model.embed, model.adj_exp)
        z = learned.detach().cpu().numpy()
        details = {
            "input_pca_embedding": np.asarray(model.embed).tolist(),
            "assignment_probabilities": np.asarray(probabilities).tolist(),
            "length_scale": length,
            "histology_used": False,
            "interpretation": "Official SpaGCN coordinate+expression mode; probabilities are model assignments, not calibrated biological confidence.",
        }
    return {
        "features": genes,
        "embedding": z.tolist(),
        "domains": [str(v) for v in labels],
        "spatial_graph_edges": w.nnz,
        **details,
    }


def _symbol_lookup(record):
    symbols = record.get("feature_symbols")
    if not symbols:
        raise SpatialError(
            "This database requires explicitly declared gene symbols; import the producer's feature mapping first."
        )
    from collections import Counter

    counts = Counter(symbols)
    # Ambiguous symbols are untestable, never silently combined or chosen arbitrarily.
    return {g: i for i, g in enumerate(symbols) if g is not None and counts[g] == 1}


def _lr(record, raw, params, seed, checkpoint):
    import liana

    resource_name = "mouseconsensus" if record["species"] == "mouse" else "consensus"
    resource = liana.resource.select_resource(resource_name)
    pairs = sorted(set((str(r.ligand), str(r.receptor)) for r in resource.itertuples()))
    resource_digest = _hash(pairs)
    requested = params.get("pairs")
    if requested is not None:
        if not isinstance(requested, list) or any(
            not isinstance(p, list) or len(p) != 2 or not all(isinstance(v, str) for v in p)
            for p in requested
        ):
            raise SpatialError("Pairs must be an explicit list of [ligand, receptor] strings.")
        requested = [tuple(p) for p in requested]
        if not requested or len(requested) != len(set(requested)) or not set(requested) <= set(pairs):
            raise SpatialError(
                "Requested ligand-receptor pairs are absent from the pinned species-specific database."
            )
        pairs = requested
    lookup = _symbol_lookup(record)
    y = normalized_counts(raw)
    w = graph(record["coordinates"], _integer(params, "neighbors", 6, 1, min(64, raw.shape[0] - 1)))
    w = sparse.diags(1 / np.asarray(w.sum(axis=1)).ravel()) @ w
    permutations = _integer(params, "permutations", 99, 19, 9999)
    if len(pairs) * raw.shape[0] * permutations > 2_000_000_000:
        raise SpatialError("LR permutation budget exceeded; explicitly restrict database pairs or ROI.")
    cache = {}

    def component(entity):
        subunits = entity.split("_")
        if not set(subunits) <= set(lookup):
            return None
        if entity not in cache:
            cache[entity] = y[:, [lookup[g] for g in subunits]].toarray().min(axis=1)
        return cache[entity]

    rows = []
    for i, (ligand, receptor) in enumerate(pairs):
        if i % 32 == 0:
            checkpoint()
        a, b = component(ligand), component(receptor)
        if a is None or b is None:
            rows.append(
                {
                    "ligand": ligand,
                    "receptor": receptor,
                    "status": "unmeasured_or_ambiguous_required_subunit",
                    "unresolved_subunits": sorted(
                        (set(ligand.split("_")) | set(receptor.split("_"))) - set(lookup)
                    ),
                    "score": None,
                    "p_value": None,
                    "q_value": None,
                }
            )
            continue
        observed = float(a @ (w @ b) / len(a))
        rng = np.random.default_rng(seed)
        exceed = sum(
            float(a @ (w @ b[rng.permutation(len(b))]) / len(a)) >= observed for _ in range(permutations)
        )
        rows.append(
            {
                "ligand": ligand,
                "receptor": receptor,
                "status": "tested",
                "score": observed,
                "ligand_detection": float(np.mean(a > 0)),
                "receptor_detection": float(np.mean(b > 0)),
                "p_value": (exceed + 1) / (permutations + 1),
                "q_value": None,
            }
        )
    tested = [r for r in rows if r["status"] == "tested"]
    for row, q in zip(tested, bh([r["p_value"] for r in tested])):
        row["q_value"] = float(q)
    rows.sort(
        key=lambda r: (
            r["q_value"] if r["q_value"] is not None else 2,
            -(r["score"] or 0),
            r["ligand"],
            r["receptor"],
        )
    )
    return {
        "rows": rows,
        "resource_name": resource_name,
        "resource_sha256": resource_digest,
        "feature_mapping": "Producer-declared exact symbols; ambiguous symbols are untestable",
        "database_pairs": len(pairs),
        "complex_rule": "minimum expression across every required subunit",
        "tested_pairs": len(tested),
        "null": "Shuffle receptor positions conditional on observed ligand pattern; library-normalized log1p RNA.",
        "multiple_testing": "BH across all measured tested pairs before display filtering",
        "interpretation": "Database-supported spatial RNA coexpression hypotheses, not observed signaling, protein binding, directionality, pathway activation or a CellChat/CellPhoneDB method reproduction.",
    }


def _pathways(record, raw, params, checkpoint, frozen_network=None):
    import decoupler as dc
    import pandas as pd

    lookup = _symbol_lookup(record)
    top = _integer(params, "top_targets", 100, 10, 1000)
    minimum = _integer(params, "min_targets", 5, 3, 100)
    license_scope = params.get("license", "academic")
    if license_scope not in {"academic", "commercial", "nonprofit"}:
        raise SpatialError("Declare a supported PROGENy resource license scope.")
    # Official species-specific footprint weights; record exact network for later review.
    net = (
        dc.op.progeny(organism=record["species"], top=top, license=license_scope)
        if frozen_network is None
        else pd.DataFrame(frozen_network)
    )
    network = net.sort_values(["source", "target"]).to_dict(orient="records")
    coverage = []
    for pathway, group in net.groupby("source"):
        measured = int(group.target.isin(lookup).sum())
        coverage.append(
            {
                "pathway": str(pathway),
                "database_targets": len(group),
                "measured_targets": measured,
                "status": "tested" if measured >= minimum else "insufficient_measured_targets",
            }
        )
    measured_symbols = sorted(set(net.target) & set(lookup))
    if not any(r["status"] == "tested" for r in coverage):
        raise SpatialError("No pathway has enough unambiguous measured targets.")
    # Include the complete unambiguous measured gene background, not only footprint targets.
    symbols = sorted(lookup)
    y = normalized_counts(raw)[:, [lookup[g] for g in symbols]]
    if y.shape[0] * y.shape[1] > 20_000_000:
        raise SpatialError("Pathway model budget exceeded; select an explicit ROI.")
    data = pd.DataFrame(y.toarray(), index=record["observation_ids"], columns=symbols)
    checkpoint()
    scores, pvalues = dc.mt.ulm(data, net, tmin=minimum, verbose=False)
    checkpoint()
    if not np.isfinite(scores.to_numpy()).all() or not np.isfinite(pvalues.to_numpy()).all():
        raise SpatialError("ULM produced nonfinite results; inspect expression/weight variance.")
    p = pvalues.to_numpy()
    return {
        "pathways": scores.columns.tolist(),
        "pathway_activity": scores.to_numpy().tolist(),
        "pathway_pvalues": p.tolist(),
        "pathway_qvalues": bh(p.ravel()).reshape(p.shape).tolist(),
        "coverage": coverage,
        "measured_unique_targets": len(measured_symbols),
        "resource": "PROGENy",
        "resource_license": license_scope,
        "resource_species_mapping": "native human resource"
        if record["species"] == "human"
        else "official decoupler human-to-mouse HCOP translation (default minimum 3 evidence resources, up to 5 orthologs); not a mouse perturbation-trained footprint",
        "resource_network": network,
        "resource_sha256": _hash(network),
        "multiple_testing": "BH across all fitted observation-by-pathway tests",
        "interpretation": "Official decoupler ULM scores with PROGENy footprints; inferred RNA pathway activity, not direct biochemical activation or evidence that a particular ligand caused it. Gene-level model p-values are not donor-level inference.",
    }


@threadpool_limits.wrap(limits=1, user_api="blas")
def run(project, spec, checkpoint=lambda: None, *, frozen_pathway_network=None):
    if not isinstance(spec, dict) or set(spec) - {"method", "input_ids", "parameters", "seed"}:
        raise SpatialError("Workflow spec accepts method, input_ids, parameters and seed only.")
    method = spec.get("method")
    if method not in METHODS:
        raise SpatialError("Unknown executable method.")
    ids = spec.get("input_ids")
    if not isinstance(ids, list) or not 1 <= len(ids) <= 8 or len(set(ids)) != len(ids):
        raise SpatialError("Provide 1..8 distinct immutable input IDs.")
    seed = _integer(spec, "seed", 0, 0, 2**32 - 1)
    params = spec.get("parameters", {})
    if not isinstance(params, dict):
        raise SpatialError("Parameters must be an object.")
    if set(params) - set(METHOD_PARAMETERS[method]):
        raise SpatialError("Unknown method parameters; misspelled options are never silently ignored.")
    records, matrices = zip(*(load_counts(project, oid) for oid in ids))
    if len({r["species"] for r in records}) != 1:
        raise SpatialError(
            "Cross-species inputs require a separately reviewed ortholog mapping, not symbol case conversion."
        )
    if method not in {"nnls", "cell2location", "harmony", "paste"} and len(ids) != 1:
        raise SpatialError("This method requires exactly one spatial input.")
    if records[0]["kind"] != "spatial":
        raise SpatialError("The first input must be a spatial sample.")
    if not util.find_spec(METHODS[method][2]):
        raise SpatialError(
            f"Missing optional method package: {METHODS[method][3]}; install the analysis environment. No fallback executed."
        )
    checkpoint()
    if method in {"mofa", "mefisto"}:
        output = _multimodal(project, records[0], matrices[0], params, method, seed, checkpoint)
    elif method in {"nnls", "cell2location"}:
        output = _deconvolve(records, matrices, params, method, seed, checkpoint)
    elif method == "harmony":
        output = _harmony(records, matrices, params, seed)
    elif method == "paste":
        output = _paste(records, matrices, params, seed)
    elif method in {"spagcn", "spatial_spectral"}:
        output = _domains(records[0], matrices[0], params, method, seed, checkpoint)
    elif method == "moran_svg":
        features = params.get("features", records[0]["features"])
        if (
            not features
            or len(features) != len(set(features))
            or not set(features) <= set(records[0]["features"])
        ):
            raise SpatialError("SVG features must be unique measured genes.")
        # Normalize against full imported panel, then select for the tested family.
        output = svg(
            matrices[0],
            records[0]["coordinates"],
            records[0]["features"],
            k=_integer(params, "neighbors", 6, 1, min(64, matrices[0].shape[0] - 1)),
            permutations=_integer(params, "permutations", 99, 19, 9999),
            seed=seed,
            min_detected=_integer(params, "min_detected", 5, 2, matrices[0].shape[0]),
            checkpoint=checkpoint,
            test_features=features,
        )
    elif method == "progeny":
        output = _pathways(records[0], matrices[0], params, checkpoint, frozen_pathway_network)
    else:
        output = _lr(records[0], matrices[0], params, seed, checkpoint)
    checkpoint()
    result = {
        "schema": "spatial-collab.analysis-result.v1",
        "method": method,
        "task": METHODS[method][1],
        "recipe": spec,
        "inputs": [
            {
                "input_id": r["object_id"],
                "input_sha256": r["object_sha256"],
                "sample_id": r["sample_id"],
                "observation_ids": r["observation_ids"],
            }
            for r in records
        ],
        "output": _clean(output),
        "environment": {**numerical_environment(), METHODS[method][3]: metadata.version(METHODS[method][3])},
        "source": METHODS[method][4],
        "scientific_authorization": "NOT_ESTABLISHED",
    }
    result["environment"]["packages"] = dict(
        sorted((d.metadata["Name"], d.version) for d in metadata.distributions())
    )
    from pathlib import Path

    result["implementation_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    from .workflow_validation import validate_result

    validate_result(project, result)
    return objects.put(project, "analysisresult", result)
