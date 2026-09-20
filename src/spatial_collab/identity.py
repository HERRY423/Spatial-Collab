"""Explicit observation namespaces and lossless stable-feature/name resolution."""
from __future__ import annotations

from urllib.parse import quote


def _identifier(value, field: str) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > 500:
        raise ValueError(f"{field} must be a nonempty string of at most 500 characters.")
    return value


def qualify_observation_id(sample_id: str, source_cell_id: str) -> str:
    """Encode both original fields without collisions from embedded delimiters."""
    sample = quote(_identifier(sample_id, "sample_id"), safe="")
    observation = quote(_identifier(source_cell_id, "source_cell_id"), safe="")
    result = f"sample:{sample}::cell:{observation}"
    if len(result) > 512:
        raise ValueError("Qualified observation ID exceeds 512 characters; choose a shorter explicit sample namespace.")
    return result


def validate_features(metadata: dict) -> list[dict]:
    """Validate the feature axis, without adding defaults to immutable metadata.

    Legacy panels treat each name as its own ID. New feature tables may repeat a
    symbol but never an ID. A null symbol means an explicitly unavailable name.
    """
    if not isinstance(metadata, dict):
        raise ValueError("Feature metadata must be an object.")
    panel = metadata.get("panel_genes")
    if not isinstance(panel, list) or any(not isinstance(v, str) or not v.strip() for v in panel):
        raise ValueError("panel_genes must contain nonempty feature IDs.")
    if len(panel) != len(set(panel)):
        raise ValueError("panel_genes must contain unique feature IDs.")
    if "features" not in metadata:
        return [{"feature_id": value, "symbol": value} for value in panel]
    features = metadata["features"]
    if not isinstance(features, list) or len(features) != len(panel):
        raise ValueError("features must describe every panel_genes feature exactly once.")
    result, seen = [], set()
    for item in features:
        if not isinstance(item, dict) or set(item) != {"feature_id", "symbol"}:
            raise ValueError("Each feature requires exactly feature_id and symbol.")
        feature_id = _identifier(item["feature_id"], "feature_id")
        symbol = item["symbol"]
        if symbol is not None:
            _identifier(symbol, "feature symbol")
        if feature_id in seen:
            raise ValueError(f"Duplicate stable feature ID: {feature_id}.")
        seen.add(feature_id)
        result.append({"feature_id": feature_id, "symbol": symbol})
    if seen != set(panel):
        raise ValueError("features IDs must match panel_genes exactly.")
    return result


def describe_feature_query(metadata: dict, query: str, *, kind: str = "auto") -> dict:
    """Describe ambiguity instead of guessing, including ID/symbol collisions.

    Prefixes ``feature_id:`` and ``symbol:`` select the namespace. Explicit kind
    overrides prefix parsing so even IDs beginning with a prefix are addressable.
    """
    _identifier(query, "feature query")
    if not isinstance(kind, str) or kind not in {"auto", "feature_id", "symbol"}:
        raise ValueError("Feature query kind must be auto, feature_id or symbol.")
    original, value = query, query
    if kind == "auto":
        for prefix in ("feature_id", "symbol"):
            if query.startswith(prefix + ":"):
                kind, value = prefix, _identifier(query[len(prefix) + 1:], "feature query")
                break
    features = validate_features(metadata)
    matches = sorted({item["feature_id"] for item in features
                      if (kind in {"auto", "feature_id"} and item["feature_id"] == value)
                      or (kind in {"auto", "symbol"} and item["symbol"] == value)})
    return {"status": "unmeasured" if not matches else "resolved" if len(matches) == 1 else "ambiguous",
            "feature_id": matches[0] if len(matches) == 1 else None,
            "candidate_feature_ids": matches, "query": original}


def resolve_feature(metadata: dict, query: str, *, kind: str = "auto") -> str:
    result = describe_feature_query(metadata, query, kind=kind)
    if result["status"] == "ambiguous":
        raise ValueError(f"Ambiguous feature query {query!r}; candidate feature IDs: {result['candidate_feature_ids']}. Use feature_id:<ID>.")
    if result["status"] == "unmeasured":
        raise ValueError(f"Unmeasured feature query {query!r}; absence is not measured zero.")
    return result["feature_id"]
