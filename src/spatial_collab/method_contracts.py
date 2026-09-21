"""Responsibility is assigned by the implementation, never by caller metadata."""
from .store import _hash

INTERNAL = {"moran_svg", "spatial_lr", "nnls", "region_comparison", "protein_region_comparison"}


def contract(method):
    internal = method in INTERNAL
    return {"schema": "spatial-collab.method-responsibility.v1", "method": method,
            "tier": "internal_numerical" if internal else "external_adapter",
            "label": "内置数值方法" if internal else "外部方法适配器",
            "responsibility": ("Own the specified statistic, preprocessing, permutation/FDR or solver objective, and numerical regression checks."
                               if internal else "Own frozen input identity, invocation, environment manifest, output retention and schema checks; not the external algorithm's scientific correctness."),
            "not_established": "Biological truth, calibrated cell identity, clinical utility or general method validity.",
            "dependency_boundary": ("Uses NumPy/SciPy numerical primitives; spatial_lr uses LIANA only as a versioned prior database, not as a LIANA inference algorithm."
                                    if internal else "Upstream implementation and model assumptions remain explicit external dependencies.")}


def execution_receipt(method, environment, input_hashes, output, implementation):
    return {"method_contract": contract(method), "environment_sha256": _hash(environment),
            "input_hashes": input_hashes, "output_sha256": _hash(output),
            "implementation_sha256": implementation,
            "verification_scope": "Content identity and adapter/numerical contract; environment hash covers the recorded manifest, not every installed binary."}
