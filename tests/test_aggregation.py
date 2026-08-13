"""Tests for post-prediction aggregation functions."""
import numpy as np
import pandas as pd
import pytest

from src.post_prediction.aggregation import (
    aggregate_variant_scores,
    score_variants_abs_max,
)


def _pred_df():
    # 3 variants, 2 model features, 2 groups, 2 samples — all values hand-chosen
    return pd.DataFrame({
        "chromosome": ["1", "1", "2"],
        "start":      [100, 200, 300],          # 1-based VCF POS
        "reference":  ["A", "AT", "G"],
        "alternate":  ["T", "A",  "C"],
        "gene_symbol":["BMPR2", "BMPR2", "SOX17"],
        "ref_pred_0": [0.1, 0.5, 0.25], "ref_pred_1": [0.2, 0.5, 0.25],
        "alt_pred_0": [0.4, 0.1, 0.25], "alt_pred_1": [0.2, 0.5, 0.75],
        "delta_0":    [0.3, -0.4, 0.0], "delta_1":    [0.0, 0.0,  0.5],
        "S1": ["0/0", "0|1", "1/1"],
        "S2": ["1/1", "./.", "."],
    })


def _make_delta_df(n_variants=10, n_features=4, seed=0):
    rng = np.random.default_rng(seed)
    data = {
        "variant_id": [f"var_{i}" for i in range(n_variants)],
        "gene_symbol": ["BMPR2"] * (n_variants // 2) + ["SOX17"] * (n_variants - n_variants // 2),
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
    # Each feature gets 3 aggregate columns (default labels are not doubled,
    # e.g. "l2_delta_0" not "l2_delta_delta_0")
    assert "max_abs_delta_0" in agg.columns
    assert "mean_abs_delta_0" in agg.columns
    assert "l2_delta_0" in agg.columns


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


def test_score_variants_abs_max():
    df = _make_delta_df(n_variants=5, n_features=3)
    scores = score_variants_abs_max(df)
    assert len(scores) == 5
    # Each score should be >= 0
    assert (scores >= 0).all()


def test_variant_weights_table_schema_and_values():
    from src.post_prediction.aggregation import build_variant_weights_table
    result = build_variant_weights_table(_pred_df(), group_col="gene_symbol")
    assert list(result["variant_id"]) == ["1:100:A:T", "1:200:AT:A", "2:300:G:C"]
    assert list(result["end"]) == [100, 201, 300]
    row1 = result.iloc[0]
    assert row1["eis_ref"] == pytest.approx(0.3)
    assert row1["eis_alt"] == pytest.approx(0.6)
    assert row1["eis_diff"] == pytest.approx(0.3)
    assert row1["abs_delta_max"] == pytest.approx(0.3)
    assert row1["abs_delta_sum"] == pytest.approx(0.3)
    assert row1["l2_delta"] == pytest.approx(0.3)
    row3 = result.iloc[2]
    assert row3["eis_diff"] == pytest.approx(0.5)
    assert row3["l2_delta"] == pytest.approx(0.5)
    assert "group" in result.columns
    assert list(result["group"]) == ["BMPR2", "BMPR2", "SOX17"]


def test_variant_weights_table_is_canonically_sorted():
    from src.post_prediction.aggregation import build_variant_weights_table
    shuffled = _pred_df().sample(frac=1, random_state=0).reset_index(drop=True)
    result = build_variant_weights_table(shuffled, group_col="gene_symbol")
    expected = build_variant_weights_table(_pred_df(), group_col="gene_symbol")
    assert list(result["variant_id"]) == list(expected["variant_id"])


def test_variant_weights_table_reports_ungrouped_variants():
    from src.post_prediction.aggregation import build_variant_weights_table
    df = _pred_df()
    df.loc[0, "gene_symbol"] = None
    result = build_variant_weights_table(df, group_col="gene_symbol")
    assert len(result) == 3
    ungrouped = result[result["group"].isna()]
    assert len(ungrouped) == 1
    assert ungrouped.iloc[0]["variant_id"] == "1:100:A:T"


def test_genotype_matrix_recoding():
    from src.post_prediction.aggregation import build_genotype_matrix
    result = build_genotype_matrix(_pred_df(), sample_ids=["S1", "S2"])
    assert list(result["variant_id"]) == ["1:100:A:T", "1:200:AT:A", "2:300:G:C"]
    assert list(result["S1"]) == [0, 1, 2]
    assert list(result["S2"]) == [2, 9, 9]
    assert result["S1"].dtype.kind in "iu"  # integer dtype


def test_genotype_matrix_rejects_unknown_call():
    from src.post_prediction.aggregation import build_genotype_matrix
    df = _pred_df()
    df.loc[0, "S1"] = "1/2"
    with pytest.raises(ValueError, match="S1"):
        build_genotype_matrix(df, sample_ids=["S1", "S2"])


def test_genotype_matrix_missing_sample_column_raises():
    from src.post_prediction.aggregation import build_genotype_matrix
    with pytest.raises(ValueError, match="S404"):
        build_genotype_matrix(_pred_df(), sample_ids=["S1", "S404"])
