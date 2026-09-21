"""Explicit assay relationships, measurement layers and molecular associations."""
from __future__ import annotations

import numpy as np

from . import objects
from .proteomics import get_assay
from .store import SpatialError, _text

RELATIONS = {"same_object", "known_mapping", "spatial_aggregation", "cross_slice_computed", "unpaired"}
LAYER_KINDS = {"measured_raw", "normalized_corrected", "model_inferred", "unmeasured"}


def assay_descriptor(project, assay_id):
    summary = project.summary()
    if assay_id == "rna":
        return {"assay_id": "rna", "modality": "rna", "sample_id": summary["metadata"].get("sample_id", summary["metadata"]["slice_id"]),
                "source_sha256": summary["source_sha256"], "features": summary["metadata"]["panel_genes"],
                "observation_ids": [c["cell_id"] for c in project.cells()], "measurement_type": "rna_count"}
    if assay_id.startswith("protein_"):
        r = get_assay(project, assay_id)
        return {"assay_id": assay_id, "modality": "protein", "sample_id": r["sample_id"],
                "source_sha256": r["assay_sha256"], "features": [f["feature_id"] for f in r["features"]],
                "observation_ids": list(r["values"]), "measurement_type": r["measurement_type"]}
    return objects.get(project, assay_id, "assaydescriptor")


def register_descriptor(project, spec):
    """Describe external modalities without claiming a loaded matrix or pairing."""
    from .integration import _ids
    for key in ("sample_id", "slice_id", "modality", "measurement_type", "coordinate_system", "units", "source_sha256"):
        _text(spec.get(key), key)
    if len(spec["source_sha256"]) != 64 or any(c not in "0123456789abcdef" for c in spec["source_sha256"]):
        raise SpatialError("Assay source needs a SHA256 digest.")
    _ids(spec.get("observation_ids"), "observation_ids")
    _ids(spec.get("features"), "features")
    return objects.put(project, "assaydescriptor", {**spec, "matrix_status": "external_descriptor_only"})


def register_correspondence(project, spec):
    kind = spec.get("relation")
    if kind not in RELATIONS:
        raise SpatialError("Declare an explicit supported correspondence relation.")
    left, right = (assay_descriptor(project, spec.get(k, "")) for k in ("left_assay", "right_assay"))
    _text(spec.get("evidence"), "correspondence evidence")
    _text(spec.get("method_version"), "correspondence method version")
    edges = spec.get("edges", [])
    if not isinstance(edges, list) or len(edges) > 500_000:
        raise SpatialError("Correspondence edges exceed budget.")
    if kind == "unpaired" and edges:
        raise SpatialError("Unpaired assays cannot contain mapping edges.")
    if kind != "unpaired" and not edges:
        raise SpatialError("Correspondence requires explicit edges, never inferred nearest neighbors.")
    lset, rset = set(left["observation_ids"]), set(right["observation_ids"])
    seen, weights = set(), {}
    for edge in edges:
        if not isinstance(edge, dict) or set(edge) != {"left_id", "right_id", "weight"}:
            raise SpatialError("Mapping edge needs left_id, right_id, weight.")
        a, b, w = edge["left_id"], edge["right_id"], edge["weight"]
        if a not in lset or b not in rset or (a, b) in seen or type(w) not in (float, int) or not np.isfinite(w) or not 0 < w <= 1:
            raise SpatialError("Invalid, foreign, repeated or unweighted correspondence edge.")
        seen.add((a, b))
        weights[a] = weights.get(a, 0) + w
    if kind in {"same_object", "known_mapping"} and (len({a for a, b in seen}) != len(seen) or len({b for a, b in seen}) != len(seen) or any(w != 1 for w in weights.values())):
        raise SpatialError("Same-object/known mappings must be one-to-one with unit weights.")
    if kind == "same_object":
        if left["sample_id"] != right["sample_id"] or any(a != b for a, b in seen):
            raise SpatialError("Same-object relation requires same sample and exact qualified IDs.")
    if kind == "spatial_aggregation" and any(not np.isclose(w, 1) for w in weights.values()):
        raise SpatialError("Aggregation weights must sum to one per left object.")
    if kind == "cross_slice_computed":
        for key in ("transform", "error_description"):
            if not spec.get(key):
                raise SpatialError("Cross-slice matching needs a recorded transform and error description.")
    tolerance = spec.get("coordinate_tolerance", 0)
    if type(tolerance) not in (float, int) or not np.isfinite(tolerance) or tolerance < 0:
        raise SpatialError("Coordinate tolerance must be explicit nonnegative finite.")
    if tolerance and not spec.get("tolerance_units"):
        raise SpatialError("Coordinate tolerance requires units; it does not establish biological identity.")
    return objects.put(project, "correspondence", {**spec, "left_source_sha256": left["source_sha256"],
                       "right_source_sha256": right["source_sha256"],
                       "unmapped_left_count": len(lset - {a for a, b in seen}),
                       "unmapped_right_count": len(rset - {b for a, b in seen}),
                       "analysis_permissions": ["paired_descriptive"] if kind == "same_object" else ["inspect_mapping_only"]})


def register_layer(project, spec):
    from .integration import _ids
    kind = spec.get("layer_kind")
    if kind not in LAYER_KINDS:
        raise SpatialError("Unknown measurement layer kind.")
    assay = assay_descriptor(project, spec.get("assay_id", ""))
    _text(spec.get("output_scale"), "output_scale")
    _text(spec.get("method"), "method")
    _text(spec.get("method_version"), "method_version")
    features = _ids(spec.get("features"), "features")
    if not set(features) <= set(assay["features"]):
        raise SpatialError("Layer features do not belong to the source assay.")
    ids = _ids(spec.get("observation_ids"), "observation_ids")
    if not set(ids) <= set(assay["observation_ids"]):
        raise SpatialError("Layer contains foreign observations.")
    required_controls = spec.get("required_controls", [])
    controls = spec.get("controls", {})
    if not isinstance(required_controls, list) or not isinstance(controls, dict) or any(not controls.get(k) for k in required_controls):
        raise SpatialError("Declared preprocessing controls are missing.")
    if spec.get("correction_claim") and not required_controls:
        raise SpatialError("A correction claim requires declared controls and their evidence.")
    values = spec.get("values")
    if kind == "unmeasured":
        if values is not None:
            raise SpatialError("Unmeasured layers cannot contain numeric values or zeros.")
    else:
        if len(ids) * len(features) > 2_000_000:
            raise SpatialError("Inline derived layer exceeds 2M entries.")
        a = np.asarray(values, dtype=float)
        if a.shape != (len(ids), len(features)) or np.isinf(a).any():
            raise SpatialError("Layer must match declared axes, using null for missing.")
        if kind == "measured_raw":
            raise SpatialError("Raw measurements enter through source import/protein registration; derived registration cannot assert new raw truth.")
        if kind == "model_inferred" and not spec.get("model_reference"):
            raise SpatialError("Inferred layers need a model reference.")
    return objects.put(project, "measurementlayer", {**spec, "source_sha256": assay["source_sha256"],
                       "negative_values_allowed": kind == "normalized_corrected", "scientific_authorization": "NOT_ESTABLISHED"})


def register_molecular_relation(project, spec):
    if spec.get("relation") not in {"same_gene_product", "complex_related", "marker_association", "researcher_defined"}:
        raise SpatialError("Declare molecular relation type; no automatic name matching.")
    for side in ("left", "right"):
        assay = assay_descriptor(project, spec.get(side + "_assay", ""))
        if spec.get(side + "_feature") not in assay["features"]:
            raise SpatialError("Molecular relation references an unmeasured feature.")
    _text(spec.get("evidence"), "molecular relationship evidence")
    _text(spec.get("declared_by"), "declared_by")
    return objects.put(project, "molecularrelation", {**spec, "status": "researcher_declared_not_independently_validated"})
