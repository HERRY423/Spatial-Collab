"""Validate task-specific axes and distinguish inferred quantities from observations."""

import numpy as np
from . import objects
from .store import SpatialError, _hash


def validate_result(project, record):
    if record.get("schema") != "spatial-collab.analysis-result.v1":
        raise SpatialError("Unknown analysis result schema.")
    if "execution_receipt" in record:
        from .method_contracts import contract
        receipt = record["execution_receipt"]
        if (record.get("method_contract") != contract(record["method"]) or
                receipt.get("method_contract") != record["method_contract"] or
                receipt.get("environment_sha256") != _hash(record["environment"]) or
                receipt.get("implementation_sha256") != record["implementation_sha256"] or
                receipt.get("output_sha256") != _hash(record["output"]) or
                receipt.get("input_hashes") != [r["input_sha256"] for r in record["inputs"]]):
            raise SpatialError("Method responsibility or execution receipt integrity mismatch.")
    n = []
    for link in record["inputs"]:
        source = objects.get(project, link["input_id"], "analysisinput")
        if (
            link["input_sha256"] != source["object_sha256"]
            or link["observation_ids"] != source["observation_ids"]
            or link["sample_id"] != source["sample_id"]
        ):
            raise SpatialError("Analysis result source identity or observation order mismatch.")
        n.append(len(link["observation_ids"]))
    if (
        record["recipe"]["input_ids"] != [r["input_id"] for r in record["inputs"]]
        or record["recipe"]["method"] != record["method"]
    ):
        raise SpatialError("Analysis result recipe differs from fitted axes/method.")
    out = record["output"]
    expected = sum(n) if record["method"] == "harmony" else n[0]
    for key in (
        "embedding",
        "uncorrected_embedding",
        "rna_contribution_fraction",
        "target_barycentric_coordinates",
        "assignment_probabilities",
        "pathway_activity",
        "pathway_pvalues",
        "pathway_qvalues",
    ):
        if key in out:
            a = np.asarray(out[key], dtype=float)
            if a.ndim != 2 or a.shape[0] != expected or a.shape[1] < 1 or not np.isfinite(a).all():
                raise SpatialError("Analysis output has nonfinite values or incorrect row axes.")
            if key in {"rna_contribution_fraction", "assignment_probabilities"} and (
                np.any(a < 0) or not np.allclose(a.sum(axis=1), 1, atol=1e-5)
            ):
                raise SpatialError("Composition/probability rows must be nonnegative and sum to one.")
            if key.startswith("pathway_") and a.shape[1] != len(out["pathways"]):
                raise SpatialError("Pathway output columns mismatch.")
            if key in {"pathway_pvalues", "pathway_qvalues"} and (np.any(a < 0) or np.any(a > 1)):
                raise SpatialError("Pathway probabilities must be in [0, 1].")
    if "domains" in out and (
        len(out["domains"]) != expected or not all(isinstance(x, str) for x in out["domains"])
    ):
        raise SpatialError("Domain output does not match spatial observation IDs.")
    if "factor_loadings" in out:
        k = len(out["embedding"][0])
        for view, size in (("RNA", len(out["features"])), ("protein", len(out["protein_features"]))):
            a = np.asarray(out["factor_loadings"][view], dtype=float)
            if a.shape != (size, k) or not np.isfinite(a).all():
                raise SpatialError("Factor loading axes must match measured features and returned factors.")
    if "abundance" in out:
        for a in out["abundance"].values():
            a = np.asarray(a, dtype=float)
            if a.shape != (expected, len(out["cell_types"])) or not np.isfinite(a).all() or np.any(a < 0):
                raise SpatialError(
                    "Abundance posterior must have finite nonnegative observation x type values."
                )
        if np.any(np.asarray(out["abundance"]["q05"]) > np.asarray(out["abundance"]["q95"])):
            raise SpatialError("Posterior interval endpoints are reversed.")
    if "transport" in out:
        t = out["transport"]
        if t["shape"] != n or len(t["row"]) != len(t["col"]) or len(t["row"]) != len(t["mass"]):
            raise SpatialError("Transport output axes mismatch.")
        mass = np.asarray(t["mass"], dtype=float)
        if any(
            type(i) is not int or not 0 <= i < n[axis]
            for axis, key in enumerate(("row", "col"))
            for i in t[key]
        ):
            raise SpatialError("Transport indices fall outside the declared input axes.")
        if not np.isfinite(mass).all() or np.any(mass < 0) or not np.isclose(mass.sum(), 1):
            raise SpatialError("Transport must preserve probability mass.")
    return True
