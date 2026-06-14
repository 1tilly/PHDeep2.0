"""Tests for post-prediction aggregation functions."""
import numpy as np
import pandas as pd
import pytest

from src.post_prediction.aggregation import (
    aggregate_variant_scores,
    build_skat_input,
    score_variants_abs_max,
)


def _make_delta_df(n_variants=10, n_features=4, seed=0):
    rng = np.random.default_rng(seed)
    data = {
        "variant_id": [f"var_{i}" for i in range(n_variants)],
        "gene_symbol": ["BMPR2"] * 5 + ["SOX17"] * 5,
        "chromosome": ["1"] * n_variants,
        "start": list(range(1000, 1000 + n_variants)),
    }
    for i in range(n_features):
        data[f"delta_{i}"] = rng.uniform(-1, 1, n_variants).tolist()
    return pd.DataFrame(data)


def test_aggregate_variant_scores_shape():
    df = _make_delta_df()
    agg = aggregate_variant_scores(df, group_col="gene_symbol")
    # 2 genes
    assert len(agg) == 2
    assert "gene_symbol" in agg.columns
    assert "n_variants" in agg.columns
    # Each feature gets 3 aggregate columns
    assert "max_abs_delta_delta_0" in agg.columns
    assert "mean_abs_delta_delta_0" in agg.columns
    assert "l2_delta_delta_0" in agg.columns


def test_aggregate_variant_scores_counts():
    df = _make_delta_df()
    agg = aggregate_variant_scores(df, group_col="gene_symbol")
    counts = dict(zip(agg["gene_symbol"], agg["n_variants"]))
    assert counts["BMPR2"] == 5
    assert counts["SOX17"] == 5


def test_aggregate_with_feature_names():
    df = _make_delta_df(n_features=2)
    agg = aggregate_variant_scores(df, group_col="gene_symbol", feature_names=["DNase", "H3K27ac"])
    assert "max_abs_delta_DNase" in agg.columns
    assert "l2_delta_H3K27ac" in agg.columns


def test_aggregate_feature_names_length_mismatch():
    df = _make_delta_df(n_features=4)
    with pytest.raises(ValueError, match="feature_names length"):
        aggregate_variant_scores(df, group_col="gene_symbol", feature_names=["only_one"])


def test_aggregate_no_delta_columns_raises():
    df = pd.DataFrame({"gene_symbol": ["A"], "some_col": [1.0]})
    with pytest.raises(ValueError, match="No delta_"):
        aggregate_variant_scores(df, group_col="gene_symbol")


def test_build_skat_input_shape():
    df = _make_delta_df(n_variants=6, n_features=3)
    matrix = build_skat_input(df, variant_id_col="variant_id")
    assert matrix.shape == (6, 3)
    assert matrix.index.tolist() == [f"var_{i}" for i in range(6)]


def test_build_skat_input_feature_names():
    df = _make_delta_df(n_variants=4, n_features=2)
    matrix = build_skat_input(df, variant_id_col="variant_id", feature_names=["F1", "F2"])
    assert list(matrix.columns) == ["F1", "F2"]


def test_score_variants_abs_max():
    df = _make_delta_df(n_variants=5, n_features=3)
    scores = score_variants_abs_max(df)
    assert len(scores) == 5
    # Each score should be >= 0
    assert (scores >= 0).all()
