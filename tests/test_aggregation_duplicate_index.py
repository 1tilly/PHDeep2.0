from __future__ import annotations

import pandas as pd

from src.post_prediction.aggregation import (
    build_genotype_matrix,
    build_variant_weights_table,
)


def _dup_index_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chromosome": ["chr1", "chr1"],
            "start": [200, 100],
            "reference": ["A", "C"],
            "alternate": ["G", "T"],
            "gene": ["g2", "g1"],
            "ref_pred_0": [1.0, 2.0],
            "alt_pred_0": [1.5, 1.0],
            "delta_0": [0.5, -1.0],
            "s1": ["0/1", "1/1"],
        },
        index=[0, 0],
    )


def test_weights_table_duplicate_index() -> None:
    out = build_variant_weights_table(_dup_index_df(), "gene")
    assert len(out) == 2
    assert list(out["variant_id"]) == ["chr1:100:C:T", "chr1:200:A:G"]
    assert list(out["group"]) == ["g1", "g2"]
    assert list(out["eis_ref"]) == [2.0, 1.0]


def test_genotype_matrix_duplicate_index() -> None:
    out = build_genotype_matrix(_dup_index_df(), ["s1"])
    assert len(out) == 2
    assert list(out["variant_id"]) == ["chr1:100:C:T", "chr1:200:A:G"]
    assert list(out["s1"]) == [2, 1]
