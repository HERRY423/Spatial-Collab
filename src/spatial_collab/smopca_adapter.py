"""Bounded adapter to the unmodified, MIT-licensed official SMOPCA model.

This tests a new runtime, not reproduction of all paper preprocessing/results.
"""
from pathlib import Path
import hashlib

import numpy as np

from .store import SpatialError

COMMIT = "9d59651d78149520661c57bab8a14d42446e6654"


def provenance():
    try:
        import sklearn
    except ImportError as exc:
        raise SpatialError("SMOPCA requires the optional integration environment (scikit-learn).") from exc
    source = Path(__file__).with_name("vendor") / "smopca_model.py"
    if not source.is_file():
        raise SpatialError("Pinned SMOPCA source is missing.")
    return {"upstream_commit": COMMIT, "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "sklearn": sklearn.__version__, "license": "MIT", "runtime_status": "local_reference_adapter"}


def fit_smopca(rna, protein, xy, components, seed, checkpoint):
    provenance()
    if len(xy) > 3000:
        raise SpatialError("SMOPCA dense-kernel budget is 3000 observations; choose an explicit ROI or import external results.")
    if components > min(rna.shape[1], protein.shape[1]):
        raise SpatialError("SMOPCA latent dimensions cannot exceed either selected modality feature count.")
    if len(np.unique(xy, axis=0)) != len(xy):
        raise SpatialError("SMOPCA reference adapter requires distinct spatial positions.")
    from .vendor.smopca_model import SMOPCA
    from scipy.spatial import cKDTree
    # Dimensionless median-neighbor units, recorded in the adapter contract.
    distances = cKDTree(xy).query(xy, k=2)[0][:, 1]
    scale = float(np.median(distances))
    if not np.isfinite(scale) or scale <= 0:
        raise SpatialError("Invalid coordinate scale for spatial kernel.")
    checkpoint()
    # Official model estimates per-feature noise near unit variance. Undo the
    # comparison baseline's block balancing, which is not SMOPCA preprocessing.
    model = SMOPCA([(rna * np.sqrt(rna.shape[1])).T, (protein * np.sqrt(protein.shape[1])).T], (xy - xy.mean(axis=0)) / scale,
                   Z_dim=components, intercept=False, omics_weight=False, kernel_type="matern", nu=1.5)
    model.estimateParams(iterations_gamma=1, iterations_sigma_W=20, tol_sigma=2e-5,
                         sigma_init_list=(1, 1), sigma_xtol_list=(1e-6, 1e-6), gamma_init=1,
                         estimate_gamma=False)
    checkpoint()
    z = np.asarray(model.calculatePosterior())
    if z.shape != (len(xy), components) or not np.isfinite(z).all():
        raise SpatialError("SMOPCA returned invalid or non-finite output; no fallback result published.")
    return z
