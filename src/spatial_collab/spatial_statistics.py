"""Sparse spatial graphs and explicitly conditional permutation statistics."""

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

from .store import SpatialError


def graph(xy, k=6):
    xy = np.asarray(xy, dtype=float)
    n = len(xy)
    if xy.shape != (n, 2) or not np.isfinite(xy).all() or not 1 <= k < n:
        raise SpatialError("Spatial graph requires finite coordinates and 1 <= k < N.")
    if len(np.unique(xy, axis=0)) != n:
        raise SpatialError("Duplicate positions need an explicit aggregation decision before graph analysis.")
    distances, neighbors = cKDTree(xy).query(xy, k=k + 1)
    distances, neighbors = distances[:, 1:], neighbors[:, 1:]
    scale = np.median(distances, axis=1)
    weights = np.exp(-0.5 * (distances / scale[:, None]) ** 2)
    w = sparse.csr_matrix((weights.ravel(), (np.repeat(np.arange(n), k), neighbors.ravel())), shape=(n, n))
    w = w.maximum(w.T)
    w.setdiag(0)
    w.eliminate_zeros()
    return w


def normalized_counts(x):
    x = sparse.csr_matrix(x, dtype=float)
    library = np.asarray(x.sum(axis=1)).ravel()
    if np.any(library <= 0):
        raise SpatialError("Zero-library rows cannot be normalized.")
    result = sparse.diags(10000 / library) @ x
    result.data = np.log1p(result.data)
    return result.tocsr()


def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p, kind="stable")
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1)
    return out


def moran_scores(matrix, weights):
    x = np.asarray(matrix, dtype=float)
    z = x - x.mean(axis=0)
    denominator = (z * z).sum(axis=0)
    scores = np.full(x.shape[1], np.nan)
    valid = denominator > 1e-12
    scores[valid] = len(x) / weights.sum() * (z * (weights @ z)).sum(axis=0)[valid] / denominator[valid]
    return scores


def svg(
    x,
    xy,
    features,
    *,
    k=6,
    permutations=99,
    seed=0,
    min_detected=5,
    checkpoint=lambda: None,
    test_features=None,
):
    if type(permutations) is not int or not 19 <= permutations <= 9999:
        raise SpatialError("Use 19..9999 permutations; low counts give coarse p-value resolution.")
    w = graph(xy, k)
    y = normalized_counts(x)
    counts = np.asarray((x > 0).sum(axis=0)).ravel()
    if test_features is not None:
        lookup = {g: i for i, g in enumerate(features)}
        columns = [lookup[g] for g in test_features]
        y, counts, features = y[:, columns], counts[columns], test_features
    selected = np.where(counts >= min_detected)[0]
    if len(selected) * len(xy) * permutations > 2_000_000_000:
        raise SpatialError("Permutation budget exceeded; explicitly restrict features/ROI or permutations.")
    rows = [
        {
            "feature_id": g,
            "detected": int(counts[i]),
            "status": "below_detection_filter",
            "morans_i": None,
            "p_value": None,
            "q_value": None,
        }
        for i, g in enumerate(features)
    ]
    for start in range(0, len(selected), 64):
        checkpoint()
        ix = selected[start : start + 64]
        dense = y[:, ix].toarray()
        observed = moran_scores(dense, w)
        exceed = np.zeros(len(ix), dtype=int)
        rng = np.random.default_rng(seed)  # Same permutations across chunks: no chunk-size effect.
        for _ in range(permutations):
            exceed += moran_scores(dense[rng.permutation(len(xy))], w) >= observed
        for j, i in enumerate(ix):
            rows[i].update(
                status="tested" if np.isfinite(observed[j]) else "constant_after_normalization",
                morans_i=float(observed[j]) if np.isfinite(observed[j]) else None,
                p_value=float((exceed[j] + 1) / (permutations + 1)) if np.isfinite(observed[j]) else None,
            )
    tested = [r for r in rows if r["p_value"] is not None]
    for row, q in zip(tested, bh([r["p_value"] for r in tested])):
        row["q_value"] = float(q)
    rows.sort(
        key=lambda r: (
            r["q_value"] if r["q_value"] is not None else 2,
            -(r["morans_i"] or 0),
            r["feature_id"],
        )
    )
    return {
        "rows": rows,
        "tested_features": len(tested),
        "graph_edges": w.nnz,
        "normalization": "full-panel library 10000 then log1p",
        "alternative": "positive spatial autocorrelation",
        "null": "Random spatial placement of expression among observed positions; not subject-level inference.",
        "multiple_testing": "BH across all tested features before display filtering",
        "minimum_p": 1 / (permutations + 1),
    }
