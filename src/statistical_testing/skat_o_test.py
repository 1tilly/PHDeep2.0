"""
SKAT-O (Sequence Kernel Association Test, Optimal) wrapper.

Wraps the R `SKAT` package via rpy2 to perform burden/variance-component
tests on aggregated variant effect scores, grouped by gene or enhancer.

Requires the `stats` optional extra:
    pip install phdeep2[stats]   # installs rpy2
And the R packages:
    install.packages(c("SKAT", "Matrix"))

Typical usage
-------------
from src.statistical_testing.skat_o_test import run_skat_o

results = run_skat_o(
    score_matrix=delta_df,   # rows=variants, cols=features (delta scores)
    phenotype=pheno_series,  # pd.Series, 1=case, 0=control
    group_col="gene_symbol", # column in delta_df that defines gene sets
)
# results: DataFrame with columns gene_symbol, p_value, p_value_burden,
#          p_value_skat, n_variants
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

try:
    from rpy2 import robjects
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.packages import importr
    _HAS_RPY2 = True
except ImportError:
    _HAS_RPY2 = False


def _require_rpy2() -> None:
    if not _HAS_RPY2:
        raise ImportError(
            "rpy2 is required for SKAT-O testing. "
            "Install with: pip install phdeep2[stats]\n"
            "Also ensure R and the SKAT package are installed."
        )


def _load_skat_package():
    """Import the R SKAT package, raising a clear error if missing."""
    try:
        return importr("SKAT")
    except Exception as exc:
        raise ImportError(
            "R package 'SKAT' not found. Install in R with:\n"
            "    install.packages('SKAT')"
        ) from exc


def run_skat_o(
    score_matrix: pd.DataFrame,
    phenotype: pd.Series,
    group_col: str,
    covariates: pd.DataFrame | None = None,
    out_type: str = "D",
    kernel: str = "linear.weighted",
) -> pd.DataFrame:
    """
    Run SKAT-O per gene/group on variant effect score matrices.

    Parameters
    ----------
    score_matrix : pd.DataFrame
        DataFrame with variants as rows. Must contain `group_col` and one or
        more feature score columns (e.g. delta_0, delta_1, ...).
    phenotype : pd.Series
        Case/control labels aligned to the sample axis. Values: 1=case, 0=control.
        For variant-level tests this is a per-variant phenotype mapped from sample
        labels. If all samples share one phenotype, pass a constant Series.
    group_col : str
        Column in score_matrix defining gene/enhancer groups (e.g. "gene_symbol").
    covariates : pd.DataFrame | None
        Optional covariate matrix (rows=samples, cols=covariates) for null model.
    out_type : str
        "D" for dichotomous (case/control), "C" for continuous phenotype.
    kernel : str
        SKAT kernel. "linear.weighted" (default) or "IBS", "IBS.weighted", etc.

    Returns
    -------
    pd.DataFrame with columns:
        group, n_variants, p_value, p_value_burden, p_value_skat
    """
    _require_rpy2()
    skat = _load_skat_package()

    numpy2ri.activate()
    pandas2ri.activate()

    delta_cols = [c for c in score_matrix.columns if c.startswith("delta_")]
    if not delta_cols:
        raise ValueError("No delta_* columns found in score_matrix.")

    y = np.array(phenotype.values, dtype=np.float64)

    # Build null model once (uses phenotype + optional covariates)
    if covariates is not None:
        cov_r = pandas2ri.py2rpy(covariates.astype(np.float64))
        null_model = skat.SKAT_Null_Model(
            robjects.Formula("y ~ ."),
            data=robjects.r["data.frame"](
                **{"y": robjects.FloatVector(y), **{c: robjects.FloatVector(covariates[c].values) for c in covariates.columns}}
            ),
            out_type=out_type,
        )
    else:
        null_model = skat.SKAT_Null_Model(
            robjects.Formula("y ~ 1"),
            data=robjects.r["data.frame"](y=robjects.FloatVector(y)),
            out_type=out_type,
        )

    rows = []
    for group_key, group_df in score_matrix.groupby(group_col):
        Z = group_df[delta_cols].values.astype(np.float64)
        if Z.shape[0] < 2:
            logger.warning("Group %s has fewer than 2 variants — skipping.", group_key)
            continue

        try:
            result = skat.SKAT(
                robjects.r.matrix(Z, nrow=Z.shape[0], ncol=Z.shape[1]),
                null_model,
                kernel=kernel,
                method="SKATO",
            )
            p_skato   = float(result.rx2("p.value")[0])
            p_burden  = float(result.rx2("p.value.burden")[0])
            p_skat    = float(result.rx2("p.value.skat")[0])
        except Exception as exc:
            logger.error("SKAT-O failed for group %s: %s", group_key, exc)
            p_skato = p_burden = p_skat = float("nan")

        rows.append({
            group_col: group_key,
            "n_variants": len(group_df),
            "p_value": p_skato,
            "p_value_burden": p_burden,
            "p_value_skat": p_skat,
        })

    numpy2ri.deactivate()
    pandas2ri.deactivate()

    result_df = pd.DataFrame(rows)
    if len(result_df) > 0:
        # Bonferroni correction
        result_df["q_value"] = np.minimum(result_df["p_value"] * len(result_df), 1.0)
        result_df = result_df.sort_values("p_value")
    return result_df


def run_skat_o_from_feather(
    scores_feather: str | Path,
    output_path: str | Path | None = None,
    group_col: str = "gene_symbol",
    phenotype_col: str = "phenotype",
) -> pd.DataFrame:
    """
    Convenience wrapper that reads an aggregated scores feather file and
    runs SKAT-O, optionally saving results.

    The feather file must contain `group_col`, `phenotype_col`, and delta_* columns.
    """
    _require_rpy2()
    df = pd.read_feather(scores_feather)

    if phenotype_col not in df.columns:
        raise ValueError(
            f"Column '{phenotype_col}' not found in {scores_feather}. "
            "The scores feather must contain per-variant phenotype labels."
        )

    phenotype = df[phenotype_col]
    results = run_skat_o(df, phenotype, group_col=group_col)

    if output_path is not None:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        results.to_feather(out)

    return results
