"""
Post-prediction aggregation and scoring.

Takes per-variant delta scores (alt_pred - ref_pred) from
VariantEffectPredictor and aggregates them into gene- or
enhancer-level summaries suitable for downstream association testing
(e.g. SKAT-O).

Functions
---------
aggregate_variant_scores
    Group delta scores by a gene/region key and compute summary
    statistics (max absolute delta, mean absolute delta, L2 norm).

build_skat_input
    Pivot the aggregated scores into the variant × feature matrix
    format expected by the SKAT-O runner.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _delta_columns(df: pd.DataFrame) -> list[str]:
    """Return all delta_<i> columns in order."""
    return [c for c in df.columns if c.startswith("delta_")]


def aggregate_variant_scores(
    delta_df: pd.DataFrame,
    group_col: str,
    feature_names: list[str] | None = None,
) -> pd.DataFrame:
    """
    Aggregate per-variant delta scores to gene/region level.

    Parameters
    ----------
    delta_df : pd.DataFrame
        Output from VariantEffectPredictor.predict_variants(). Must contain
        a grouping column (e.g. "gene_symbol") and delta_<i> columns.
    group_col : str
        Column to group by (e.g. "gene_symbol", "genehancer_id").
    feature_names : list[str] | None
        Optional feature labels for the output columns. If None, column names
        remain delta_0, delta_1, …

    Returns
    -------
    pd.DataFrame
        One row per group with columns:
            <group_col>, n_variants,
            max_abs_delta_<feature>, mean_abs_delta_<feature>, l2_delta_<feature>
    """
    delta_cols = _delta_columns(delta_df)
    if not delta_cols:
        raise ValueError("No delta_* columns found in delta_df.")

    if feature_names is not None and len(feature_names) != len(delta_cols):
        raise ValueError(
            f"feature_names length ({len(feature_names)}) must match "
            f"number of delta columns ({len(delta_cols)})."
        )

    label = feature_names or delta_cols
    rows = []

    for group_key, group in delta_df.groupby(group_col):
        deltas = group[delta_cols].values.astype(np.float32)  # (n_variants, n_features)
        row: dict = {group_col: group_key, "n_variants": len(group)}
        for i, col_label in enumerate(label):
            col = deltas[:, i]
            row[f"max_abs_delta_{col_label}"] = float(np.max(np.abs(col)))
            row[f"mean_abs_delta_{col_label}"] = float(np.mean(np.abs(col)))
            row[f"l2_delta_{col_label}"] = float(np.sqrt(np.sum(col ** 2)))
        rows.append(row)

    return pd.DataFrame(rows)


def build_skat_input(
    delta_df: pd.DataFrame,
    variant_id_col: str,
    feature_names: list[str] | None = None,
) -> pd.DataFrame:
    """
    Pivot delta scores into a variant × feature matrix for SKAT-O.

    Parameters
    ----------
    delta_df : pd.DataFrame
        Per-variant delta scores (output of VariantEffectPredictor).
    variant_id_col : str
        Column that uniquely identifies each variant (e.g. "variant_id").
    feature_names : list[str] | None
        Optional column labels for the feature matrix.

    Returns
    -------
    pd.DataFrame, shape (n_variants, n_features)
        Index set to variant_id_col; columns named by feature_names or
        delta_0, delta_1, …
    """
    delta_cols = _delta_columns(delta_df)
    if not delta_cols:
        raise ValueError("No delta_* columns found in delta_df.")

    matrix = delta_df[[variant_id_col] + delta_cols].set_index(variant_id_col)
    if feature_names is not None:
        if len(feature_names) != len(delta_cols):
            raise ValueError(
                f"feature_names length ({len(feature_names)}) must match "
                f"number of delta columns ({len(delta_cols)})."
            )
        matrix.columns = feature_names
    return matrix


def score_variants_abs_max(delta_df: pd.DataFrame) -> pd.Series:
    """
    Compute a single per-variant score as the maximum absolute delta
    across all features.

    Useful as a simple variant-level effect size for ranking.

    Returns
    -------
    pd.Series, same index as delta_df
    """
    delta_cols = _delta_columns(delta_df)
    return delta_df[delta_cols].abs().max(axis=1)
