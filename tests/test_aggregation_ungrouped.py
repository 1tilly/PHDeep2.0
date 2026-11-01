"""Regression test: aggregate_variant_scores must not drop ungrouped variants."""
from __future__ import annotations

import pandas as pd

from src.post_prediction.aggregation import aggregate_variant_scores


def test_ungrouped_variants_are_retained() -> None:
    df = pd.DataFrame(
        {
            "gene_symbol": ["A", "A", None],
            "delta_0": [0.5, -1.0, 2.0],
        }
    )

    result = aggregate_variant_scores(df, "gene_symbol")

    assert int(result["n_variants"].sum()) == 3
    assert len(result) == 2
    null_rows = result[result["gene_symbol"].isna()]
    assert len(null_rows) == 1
    assert int(null_rows["n_variants"].iloc[0]) == 1
    assert float(null_rows["max_abs_delta_0"].iloc[0]) == 2.0
