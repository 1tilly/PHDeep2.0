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

build_variant_weights_table
    Build a one-row-per-variant table of variant weights (effect
    is-score sums and delta-based magnitudes), sorted canonically and
    annotated with a group label. This is the per-variant weight
    handoff SKAT-O needs (as opposed to aggregate_variant_scores,
    which collapses to one row per group).

build_genotype_matrix
    Build a variant × sample genotype matrix (0/1/2/9 dosage coding)
    aligned to the same variant order as build_variant_weights_table.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd


def _delta_columns(df: pd.DataFrame) -> list[str]:
    """Return all delta_<i> columns in order."""
    return [c for c in df.columns if c.startswith("delta_")]


def _ref_pred_columns(df: pd.DataFrame) -> list[str]:
    """Return all ref_pred_<i> columns in order."""
    return [c for c in df.columns if c.startswith("ref_pred_")]


def _alt_pred_columns(df: pd.DataFrame) -> list[str]:
    """Return all alt_pred_<i> columns in order."""
    return [c for c in df.columns if c.startswith("alt_pred_")]


def _variant_id(chromosome: str, start: int, reference: str, alternate: str) -> str:
    """Build a canonical variant id string ``chrom:start:ref:alt``."""
    return f"{chromosome}:{start}:{reference}:{alternate}"


# Recognized homozygous-alt calls
_HOM_ALT = {"1/1", "1|1"}
# Recognized het calls
_HET = {"0/1", "1/0", "0|1", "1|0"}
# Recognized hom-ref calls
_HOM_REF = {"0/0", "0|0"}
# Recognized missing calls
_MISSING_RE = re.compile(r"^\.([/|]\.)?$")


def _recode_genotype(value: object) -> int:
    """Recode a genotype string to a 0/1/2/9 dosage, or raise if unrecognized."""
    text = str(value).strip()
    if text in _HOM_REF:
        return 0
    if text in _HET:
        return 1
    if text in _HOM_ALT:
        return 2
    if _MISSING_RE.match(text):
        return 9
    raise ValueError(f"Unrecognized genotype call {value!r}")


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
        use the delta column's numeric suffix (e.g. "0", "1", …).

    Returns
    -------
    pd.DataFrame
        One row per group with columns:
            <group_col>, n_variants,
            max_abs_delta_<feature>, mean_abs_delta_<feature>, l2_delta_<feature>
        where <feature> defaults to "0", "1", … (from delta_0, delta_1, …)
        or the supplied feature_names.
    """
    delta_cols = _delta_columns(delta_df)
    if not delta_cols:
        raise ValueError("No delta_* columns found in delta_df.")

    if feature_names is not None and len(feature_names) != len(delta_cols):
        raise ValueError(
            f"feature_names length ({len(feature_names)}) must match "
            f"number of delta columns ({len(delta_cols)})."
        )

    # Default labels strip the redundant "delta_" prefix (e.g. "delta_0" -> "0")
    # so output columns read "l2_delta_0" rather than "l2_delta_delta_0".
    label = feature_names or [c.removeprefix("delta_") for c in delta_cols]
    rows = []

    for group_key, group in delta_df.groupby(group_col):
        deltas = group[delta_cols].values.astype(np.float64)  # (n_variants, n_features)
        row: dict = {group_col: group_key, "n_variants": len(group)}
        for i, col_label in enumerate(label):
            col = deltas[:, i]
            row[f"max_abs_delta_{col_label}"] = float(np.max(np.abs(col)))
            row[f"mean_abs_delta_{col_label}"] = float(np.mean(np.abs(col)))
            row[f"l2_delta_{col_label}"] = float(np.sqrt(np.sum(col ** 2)))
        rows.append(row)

    return pd.DataFrame(rows)


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


def _canonical_variant_frame(pred_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build a chromosome/start/end/reference/alternate/variant_id frame,
    indexed like pred_df but sorted canonically by
    (chromosome, start, end, alternate).
    """
    reference = pred_df["reference"]
    frame = pd.DataFrame(
        {
            "chromosome": pred_df["chromosome"],
            "start": pred_df["start"],
            "end": pred_df["start"] + reference.str.len() - 1,
            "reference": reference,
            "alternate": pred_df["alternate"],
        },
        index=pred_df.index,
    )
    frame["variant_id"] = [
        _variant_id(c, s, r, a)
        for c, s, r, a in zip(
            frame["chromosome"], frame["start"], frame["reference"], frame["alternate"]
        )
    ]
    return frame.sort_values(["chromosome", "start", "end", "alternate"])


def build_variant_weights_table(pred_df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    """
    Build a one-row-per-variant weights table for downstream association
    testing (e.g. SKAT-O), which needs a weight PER VARIANT rather than
    the per-group summary produced by aggregate_variant_scores.

    Rows whose group_col value is null/NaN are RETAINED (not dropped)
    with "group" set to NaN; callers that need to report/exclude
    ungrouped variants can check `result["group"].isna()`.

    Parameters
    ----------
    pred_df : pd.DataFrame
        Output from VariantEffectPredictor.predict_variants(). Must
        contain "chromosome", "start", "reference", "alternate",
        group_col, ref_pred_<i>, alt_pred_<i>, and delta_<i> columns.
    group_col : str
        Column identifying the gene/region each variant belongs to.

    Returns
    -------
    pd.DataFrame
        One row per variant, sorted by (chromosome, start, end,
        alternate), with columns:
            chromosome, start, end, reference, alternate, variant_id,
            group, n_features,
            eis_ref, eis_alt, eis_diff,
            abs_delta_max, abs_delta_sum, l2_delta
    """
    delta_cols = _delta_columns(pred_df)
    ref_cols = _ref_pred_columns(pred_df)
    alt_cols = _alt_pred_columns(pred_df)
    if not delta_cols:
        raise ValueError("No delta_* columns found in pred_df.")

    frame = _canonical_variant_frame(pred_df)
    pred_sorted = pred_df.loc[frame.index].reset_index(drop=True)
    out = frame.reset_index(drop=True).copy()

    out["group"] = pred_sorted[group_col]
    out["n_features"] = len(delta_cols)

    ref_vals = pred_sorted[ref_cols].to_numpy(dtype=np.float64)
    alt_vals = pred_sorted[alt_cols].to_numpy(dtype=np.float64)
    delta_vals = pred_sorted[delta_cols].to_numpy(dtype=np.float64)

    out["eis_ref"] = ref_vals.sum(axis=1)
    out["eis_alt"] = alt_vals.sum(axis=1)
    out["eis_diff"] = out["eis_alt"] - out["eis_ref"]
    out["abs_delta_max"] = np.max(np.abs(delta_vals), axis=1)
    out["abs_delta_sum"] = np.sum(np.abs(delta_vals), axis=1)
    out["l2_delta"] = np.sqrt(np.sum(delta_vals**2, axis=1))

    return out


def build_genotype_matrix(pred_df: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """
    Build a variant x sample genotype (dosage) matrix for SKAT-O, in the
    same variant order as build_variant_weights_table.

    Genotype strings are recoded to integer dosages: 0/0 -> 0; any
    single-copy-alt het (0/1, 1/0, and the phased 0|1, 1|0) -> 1; 1/1
    (or phased 1|1) -> 2; missing calls ("." or "./." or ".|.") -> 9.
    Both "/" and "|" separators are accepted. Any other call (e.g.
    multi-allelic "1/2") raises ValueError.

    Parameters
    ----------
    pred_df : pd.DataFrame
        Must contain "chromosome", "start", "reference", "alternate",
        and one column per sample_id holding genotype call strings.
    sample_ids : list[str]
        Sample columns to include, in output-column order.

    Returns
    -------
    pd.DataFrame
        Columns: variant_id, <sample_id> (integer dosage) for each
        sample_id, sorted by (chromosome, start, end, alternate).
    """
    missing_cols = [s for s in sample_ids if s not in pred_df.columns]
    if missing_cols:
        raise ValueError(
            f"Sample column(s) not found in pred_df: {', '.join(missing_cols)}"
        )

    frame = _canonical_variant_frame(pred_df)
    pred_sorted = pred_df.loc[frame.index].reset_index(drop=True)
    frame_sorted = frame.reset_index(drop=True)

    out = pd.DataFrame({"variant_id": frame_sorted["variant_id"]})
    for sample_id in sample_ids:
        codes = np.empty(len(pred_sorted), dtype=np.int64)
        for i, (variant_id, value) in enumerate(
            zip(frame_sorted["variant_id"], pred_sorted[sample_id])
        ):
            try:
                codes[i] = _recode_genotype(value)
            except ValueError as exc:
                raise ValueError(
                    f"Unrecognized genotype call {value!r} for sample "
                    f"{sample_id!r} at variant {variant_id!r}"
                ) from exc
        out[sample_id] = codes

    return out
