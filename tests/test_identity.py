"""Cross-sample identity and feature ambiguity have explicit, lossless semantics."""
from copy import deepcopy

import pytest

from spatial_collab.identity import (describe_feature_query, qualify_observation_id,
                                     resolve_feature, validate_features)


def metadata():
    return {"panel_genes": ["E1", "E2", "E3", "E4"],
            "features": [{"feature_id": "E1", "symbol": "Gene"},
                         {"feature_id": "E2", "symbol": "Gene"},
                         {"feature_id": "E3", "symbol": "E1"},
                         {"feature_id": "E4", "symbol": None}]}


def test_observation_namespaces_escape_delimiters_and_unicode():
    values = [qualify_observation_id(sample, source) for sample, source in
              [("s1", "a"), ("s2", "a"), ("s:a", "b"), ("s", "a:b"),
               ("s%3Aa", "b"), ("组织", "条码")]]
    assert len(values) == len(set(values))
    assert qualify_observation_id("s1", "a") == "sample:s1::cell:a"


@pytest.mark.parametrize("sample,source", [(None, "a"), (True, "a"), ("", "a"), ("s", " "),
                                           ("s", 3), ("s" * 499, "c" * 499)])
def test_observation_invalid_namespace_rejected(sample, source):
    with pytest.raises(ValueError):
        qualify_observation_id(sample, source)


def test_symbols_are_not_summed_and_cross_namespace_collision_is_ambiguous():
    meta = metadata()
    assert describe_feature_query(meta, "Gene") == {
        "status": "ambiguous", "feature_id": None, "candidate_feature_ids": ["E1", "E2"], "query": "Gene"}
    assert describe_feature_query(meta, "E1")["candidate_feature_ids"] == ["E1", "E3"]
    with pytest.raises(ValueError, match="Ambiguous.*E1.*E2"):
        resolve_feature(meta, "Gene")
    assert resolve_feature(meta, "feature_id:E1") == "E1"
    assert resolve_feature(meta, "symbol:E1") == "E3"
    assert resolve_feature(meta, "E1", kind="feature_id") == "E1"


def test_legacy_mapping_and_missing_are_nonmutating():
    meta = {"panel_genes": ["Cd3d", "Cd79a"]}
    before = deepcopy(meta)
    assert resolve_feature(meta, "Cd3d") == "Cd3d"
    assert describe_feature_query(meta, "Cd4") == {
        "status": "unmeasured", "feature_id": None, "candidate_feature_ids": [], "query": "Cd4"}
    with pytest.raises(ValueError, match="Unmeasured"):
        resolve_feature(meta, "Cd4")
    assert meta == before


def test_ids_with_namespace_prefix_still_addressable():
    meta = {"panel_genes": ["symbol:foo"]}
    assert resolve_feature(meta, "symbol:foo", kind="feature_id") == "symbol:foo"
    assert resolve_feature(meta, "feature_id:symbol:foo") == "symbol:foo"


@pytest.mark.parametrize("change", ["duplicate_id", "missing_id", "unknown_id", "invalid_symbol", "extra_key"])
def test_feature_axis_validation_rejects_lossy_alignment(change):
    meta = metadata()
    if change == "duplicate_id":
        meta["features"][1]["feature_id"] = "E1"
    elif change == "missing_id":
        meta["features"].pop()
    elif change == "unknown_id":
        meta["features"][1]["feature_id"] = "E99"
    elif change == "invalid_symbol":
        meta["features"][1]["symbol"] = 12
    else:
        meta["features"][1]["inferred"] = True
    with pytest.raises(ValueError):
        validate_features(meta)


@pytest.mark.parametrize("kind", [[], {}, True, "guess"])
def test_query_kind_is_strict(kind):
    with pytest.raises(ValueError):
        describe_feature_query(metadata(), "Gene", kind=kind)
