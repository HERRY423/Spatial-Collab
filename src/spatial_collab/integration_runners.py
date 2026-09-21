"""Spatial graph representation runner and Moran's I spatial feature selection.

Integrates spatial topology via adaptive nearest-neighbor graph smoothing,
overcoming non-spatial baseline limits and small spot-count constraints.
"""
from __future__ import annotations

import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.linalg import svd
from scipy.spatial import cKDTree
from scipy import sparse

from .store import SpatialError


def compute_spatial_weights(xy_coords: np.ndarray, k: int = 6, sparse_output: bool = False):
    """Build row-stochastic Gaussian spatial adjacency matrix from 2D coordinates."""
    xy = np.asarray(xy_coords, dtype=float)
    if xy.ndim != 2 or xy.shape[1] != 2 or not np.isfinite(xy).all():
        raise SpatialError("xy_coords must be a finite N x 2 array.")
    n = len(xy)
    if type(k) is not int or k < 1:
        raise SpatialError("k must be a positive integer.")
    if not sparse_output and n > 4096:
        raise SpatialError("Dense graph previews are limited to 4096 positions; use sparse_output=True.")
    if len(np.unique(xy, axis=0)) != n:
        raise SpatialError("Duplicate spatial positions require an explicit aggregation decision.")
    if n < 2:
        return np.eye(n, dtype=float)
    effective_k = min(k + 1, n)
    tree = cKDTree(xy)
    distances, indices = tree.query(xy, k=effective_k)
    # Exclude self (first column)
    if effective_k > 1:
        nbr_indices = indices[:, 1:]
        nbr_dists = distances[:, 1:]
    else:
        nbr_indices = indices
        nbr_dists = distances

    # Gaussian kernel with bandwidth proportional to local median distance
    bandwidth = np.median(nbr_dists, axis=1, keepdims=True)
    bandwidth[bandwidth <= 1e-12] = 1.0

    raw_weights = np.exp(-0.5 * (nbr_dists / bandwidth) ** 2)
    # Row normalization
    row_sums = raw_weights.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    norm_weights = raw_weights / row_sums

    W = sparse.csr_matrix((norm_weights.ravel(), (np.repeat(np.arange(n), nbr_indices.shape[1]), nbr_indices.ravel())), shape=(n, n))
    return W if sparse_output else W.toarray()


def smooth_features(matrix: np.ndarray, W: np.ndarray, alpha: float = 0.65) -> np.ndarray:
    """Smooth feature matrix along spatial graph: Z_smooth = alpha * Z + (1 - alpha) * W @ Z."""
    if not (0.0 < alpha <= 1.0):
        raise SpatialError("alpha must be in (0, 1].")
    smoothed = alpha * matrix + (1.0 - alpha) * (W @ matrix)
    return smoothed


def compute_morans_i(matrix: np.ndarray, W: np.ndarray) -> np.ndarray:
    """Compute Moran's I spatial autocorrelation for each column of matrix.

    I = (N / sum(W)) * (sum_i sum_j W_ij (x_i - mean)(x_j - mean)) / sum_i (x_i - mean)^2
    """
    X = np.asarray(matrix, dtype=float)
    n = X.shape[0]
    if n < 3:
        return np.zeros(X.shape[1], dtype=float)
    means = X.mean(axis=0, keepdims=True)
    centered = X - means
    s0 = (centered**2).sum(axis=0)
    w_sum = W.sum()
    if w_sum <= 1e-12:
        return np.zeros(X.shape[1], dtype=float)

    lag = W @ centered
    numerator = n * (centered * lag).sum(axis=0)
    denominator = w_sum * s0

    moran = np.zeros(X.shape[1], dtype=float)
    valid = denominator > 1e-12
    moran[valid] = numerator[valid] / denominator[valid]
    return moran


def select_spatial_features(cells: list[dict], panel_features: list[str],
                           xy_coords: np.ndarray, max_features: int = 128) -> list[str]:
    """Select top spatially variable genes ranked by Moran's I spatial autocorrelation."""
    if len(panel_features) <= max_features:
        return sorted(panel_features)
    W = compute_spatial_weights(xy_coords, k=6, sparse_output=True)
    moran_scores = np.empty(len(panel_features))
    for start in range(0, len(panel_features), 64):
        features = panel_features[start:start+64]
        rna_matrix = np.array([[c["counts"].get(f, 0.0) for f in features] for c in cells], dtype=float)
        moran_scores[start:start+len(features)] = compute_morans_i(rna_matrix, W)
    indexed = sorted(range(len(panel_features)),
                     key=lambda idx: (-moran_scores[idx], panel_features[idx]))
    return [panel_features[i] for i in indexed[:max_features]]


def fit_spatial_graph(rna: np.ndarray, protein: np.ndarray, xy_coords: np.ndarray,
                      components: int, clusters: int, seed: int, checkpoint=None,
                      alpha: float = 0.65) -> tuple[dict, dict, dict]:
    """Perform spatial graph-smoothed multi-modal representation and domain partitioning."""
    check = checkpoint or (lambda: None)
    W = compute_spatial_weights(xy_coords, k=6, sparse_output=True)
    check()

    rna_smooth = smooth_features(rna, W, alpha=alpha)
    protein_smooth = smooth_features(protein, W, alpha=alpha)
    joint_smooth = np.concatenate([rna_smooth, protein_smooth], axis=1)

    representations, domains, models = {}, {}, {}
    for name, matrix in (("rna_only", rna_smooth),
                         ("protein_only", protein_smooth),
                         ("joint", joint_smooth)):
        check()
        u, s, vt = svd(matrix, full_matrices=False, check_finite=False)
        n = min(components, len(s))
        embedding = u[:, :n] * s[:n]
        check()
        if len(np.unique(embedding, axis=0)) < clusters:
            raise SpatialError(f"{name}: fewer distinct spatial graph representations than requested clusters.")
        centers, labels = kmeans2(embedding, clusters, iter=50, minit="++", seed=seed, missing="raise")
        representations[name] = embedding.tolist()
        domains[name] = [str(int(x)) for x in labels]
        models[name] = {"loadings": vt[:n].tolist(), "centers": centers.tolist(), "singular_values": s[:n].tolist(),
                        "alpha": alpha, "spatial_neighbors": 6}

    return representations, domains, models
