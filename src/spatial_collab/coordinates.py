"""Explicit coordinate capabilities; unknown scale never becomes micrometres."""
UNITS = {"micrometer", "pixel", "array_index", "unknown"}


def coordinate_capabilities(metadata):
    units = metadata.get("units", "unknown")
    physical = units == "micrometer"
    return {"units": units, "coordinate_system": metadata.get("coordinate_system"),
            "physical_radius_analysis": physical, "centroid_exploration": True,
            "expression_review": True, "annotation_review": True, "region_comparison": True,
            "calibration_status": "declared_micrometer_frame" if physical else "physical_scale_unverified",
            "reason": ("Physical radii use the explicitly declared micrometer frame; metadata is not external verification."
                       if physical else "Coordinates remain in their recorded units. Browsing, marker and region review are available; micrometer radius analysis requires evidence-backed calibration.")}
