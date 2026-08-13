"""
SKAT-O (Sequence Kernel Association Test, Optimal) wrapper.

Wraps the R `SKAT` package via rpy2 to perform burden/variance-component
tests on aggregated variant effect scores, grouped by gene or enhancer.

Requires the `stats` optional extra:
    pip install phdeep2[stats]   # installs rpy2
And the R packages:
    install.packages(c("SKAT", "Matrix"))

All R interaction is routed through the `SkatBackend` protocol so that the
module (and its pure-Python input-building logic) can be imported and tested
without rpy2 or R installed. The default backend (`RpyBackend`) is only
instantiated lazily, inside `run_skat_o`, when no other backend is supplied.

Typical usage
-------------
from src.statistical_testing.skat_o_test import run_skat_o

results = run_skat_o(
    variant_table=weights_df,   # rows=variants, has group_col + weight columns
    genotypes=genotype_df,      # rows=variants (variant_id col), cols=samples
    phenotype=pheno_series,     # pd.Series indexed by sample id, 1=case, 0=control
    weight_col="eis_diff",
    group_col="group",
)
# results: DataFrame with columns feature_id, n_variants, n_samples, p_value,
#          p_value_burden, p_value_skat, q_value, weight
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

try:
    from rpy2 import robjects
    from rpy2.robjects import numpy2ri
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


@dataclass
class SkatGroupInput:
    """Per-group inputs ready to hand to a SKAT-O backend."""

    group: str
    Z: np.ndarray  # samples x variants
    y: np.ndarray  # per-sample phenotype, in sample_order
    weights: np.ndarray | None  # per-variant, aligned to Z's columns
    variant_ids: list[str]  # ordered, aligned to Z's columns


def build_skat_o_inputs(
    variant_table: pd.DataFrame,
    genotypes: pd.DataFrame,
    phenotype: pd.Series,
    weight_col: str | None = None,
    group_col: str = "group",
    min_variants: int = 1,
) -> list[SkatGroupInput]:
    """
    Build per-group (samples x variants) inputs for SKAT-O from a real
    genotype matrix and a real per-sample phenotype.

    Parameters
    ----------
    variant_table : pd.DataFrame
        Must contain `variant_id`, `group_col`, and (if `weight_col` given)
        a numeric weight column.
    genotypes : pd.DataFrame
        Must contain `variant_id` plus one integer column per sample, values
        in {0, 1, 2, 9} (9 = missing).
    phenotype : pd.Series
        Indexed by sample id, values 1 (case) / 0 (control). Its index order
        defines the canonical sample order for the analysis.
    weight_col : str | None
        Column in `variant_table` to use as per-variant weights. If None,
        no weights are built (`SkatGroupInput.weights` is None).
    group_col : str
        Column in `variant_table` defining testable units (e.g. gene).
    min_variants : int
        Groups with fewer variants than this are skipped.

    Returns
    -------
    list[SkatGroupInput]
    """
    sample_cols = [c for c in genotypes.columns if c != "variant_id"]

    missing_samples = [s for s in phenotype.index if s not in sample_cols]
    if missing_samples:
        raise ValueError(
            f"Sample(s) {missing_samples} present in phenotype but missing "
            "from genotypes columns."
        )

    sample_order = [s for s in phenotype.index if s in sample_cols]
    if not sample_order:
        raise ValueError(
            "No overlap between phenotype index and genotype sample columns."
        )

    valid_genotype_values = {0, 1, 2, 9}
    geno_values = genotypes[sample_order].to_numpy()
    if not np.isin(geno_values, list(valid_genotype_values)).all():
        bad_mask = ~np.isin(geno_values, list(valid_genotype_values))
        bad_row, bad_col = np.argwhere(bad_mask)[0]
        bad_variant = genotypes.iloc[bad_row]["variant_id"]
        bad_sample = sample_order[bad_col]
        raise ValueError(
            f"Invalid genotype value {geno_values[bad_row, bad_col]!r} for "
            f"variant {bad_variant!r}, sample {bad_sample!r}. "
            f"Allowed values: {sorted(valid_genotype_values)}."
        )

    if weight_col is not None:
        weight_values = variant_table[weight_col]
        if weight_values.isna().any():
            nan_variants = variant_table.loc[weight_values.isna(), "variant_id"].tolist()
            raise ValueError(
                f"NaN weight value(s) in column {weight_col!r} for variant(s) "
                f"{nan_variants}."
            )
        if (weight_values < 0).any():
            neg_variants = variant_table.loc[weight_values < 0, "variant_id"].tolist()
            raise ValueError(
                f"Negative weight value(s) in column {weight_col!r} for "
                f"variant(s) {neg_variants}. SKAT squares weights internally, "
                "so signed scores need an explicit transform upstream "
                "(e.g. abs())."
            )

    genotypes_indexed = genotypes.set_index("variant_id")
    if weight_col is not None:
        weights_indexed = variant_table.set_index("variant_id")[weight_col]

    results: list[SkatGroupInput] = []
    for group_key, group_df in variant_table.groupby(group_col):
        group_variant_ids = group_df["variant_id"].tolist()

        missing_variants = [
            v for v in group_variant_ids if v not in genotypes_indexed.index
        ]
        if missing_variants:
            raise ValueError(
                f"Variant id(s) {missing_variants} present in variant_table "
                "but missing from genotypes."
            )

        if len(group_variant_ids) < min_variants:
            continue

        Z = (
            genotypes_indexed.loc[group_variant_ids, sample_order]
            .to_numpy(dtype=float)
            .T
        )
        y = phenotype.loc[sample_order].to_numpy()
        weights = (
            weights_indexed.loc[group_variant_ids].to_numpy()
            if weight_col is not None
            else None
        )

        results.append(
            SkatGroupInput(
                group=str(group_key),
                Z=Z,
                y=y,
                weights=weights,
                variant_ids=group_variant_ids,
            )
        )

    return results


def bh_fdr(p_values: Any) -> np.ndarray:
    """
    Benjamini-Hochberg FDR correction.

    Standard algorithm: sort p-values ascending, compute
    ``adjusted[i] = p[i] * n / rank[i]`` (rank is 1-based), then enforce
    monotonicity by taking the running minimum from the largest p-value down
    to the smallest, and clip to [0, 1]. Returns values in the same order as
    the input (not sorted).
    """
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranks = np.arange(1, n + 1)
    sorted_p = p[order]

    adjusted_sorted = sorted_p * n / ranks
    # running minimum from the largest p-value down to the smallest
    adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]
    adjusted_sorted = np.clip(adjusted_sorted, 0.0, 1.0)

    adjusted = np.empty(n, dtype=float)
    adjusted[order] = adjusted_sorted
    return adjusted


class SkatBackend(Protocol):
    """Backend interface for running SKAT-O's R routines."""

    def run_null_model(self, y: Any, X: Any, out_type: str) -> Any:
        """Fit and return a SKAT null model."""
        ...

    def run_skat(
        self, Z: np.ndarray, null_model: Any, weights: np.ndarray | None, method: str
    ) -> dict:
        """Run SKAT-O for one group. Returns a dict with at least
        `p_value`, `p_value_burden`, `p_value_skat` keys."""
        ...


class RpyBackend:
    """SkatBackend implementation that drives the R `SKAT` package via rpy2."""

    def __init__(self) -> None:
        _require_rpy2()
        self._skat = self._load_skat_package()

    @staticmethod
    def _load_skat_package():
        try:
            return importr("SKAT")
        except Exception as exc:
            raise ImportError(
                "R package 'SKAT' not found. Install in R with:\n"
                "    install.packages('SKAT')"
            ) from exc

    def run_null_model(self, y: Any, X: Any, out_type: str) -> Any:
        y_arr = np.asarray(y, dtype=np.float64)
        with (robjects.default_converter + numpy2ri.converter).context():
            if X is not None:
                data_kwargs = {"y": robjects.FloatVector(y_arr)}
                data_kwargs.update(
                    {c: robjects.FloatVector(X[c].to_numpy(dtype=np.float64)) for c in X.columns}
                )
                null_model = self._skat.SKAT_Null_Model(
                    robjects.Formula("y ~ ."),
                    data=robjects.r["data.frame"](**data_kwargs),
                    out_type=out_type,
                )
            else:
                null_model = self._skat.SKAT_Null_Model(
                    robjects.Formula("y ~ 1"),
                    data=robjects.r["data.frame"](y=robjects.FloatVector(y_arr)),
                    out_type=out_type,
                )
        return null_model

    def run_skat(
        self, Z: np.ndarray, null_model: Any, weights: np.ndarray | None, method: str
    ) -> dict:
        with (robjects.default_converter + numpy2ri.converter).context():
            r_matrix = robjects.r.matrix(Z, nrow=Z.shape[0], ncol=Z.shape[1])
            kwargs: dict[str, Any] = {"method": method}
            if weights is not None:
                kwargs["weights"] = robjects.FloatVector(np.asarray(weights, dtype=np.float64))
            result = self._skat.SKAT(r_matrix, null_model, **kwargs)
            return {
                "p_value": float(result.rx2("p.value")[0]),
                "p_value_burden": float(result.rx2("p.value.burden")[0]),
                "p_value_skat": float(result.rx2("p.value.skat")[0]),
            }


def run_skat_o(
    variant_table: pd.DataFrame,
    genotypes: pd.DataFrame,
    phenotype: pd.Series,
    covariates: pd.DataFrame | None = None,
    weight_col: str | None = None,
    group_col: str = "group",
    out_type: str = "D",
    method: str = "optimal.adj",
    min_variants: int = 1,
    backend: "SkatBackend | None" = None,
) -> pd.DataFrame:
    """
    Run SKAT-O per gene/group on a real genotype matrix and per-sample
    phenotype.

    Parameters
    ----------
    variant_table : pd.DataFrame
        Must contain `variant_id`, `group_col`, and (if `weight_col` given)
        a numeric weight column.
    genotypes : pd.DataFrame
        Must contain `variant_id` plus one integer column per sample, values
        in {0, 1, 2, 9} (9 = missing).
    phenotype : pd.Series
        Indexed by sample id, values 1 (case) / 0 (control). Its index order
        defines the canonical sample order for the analysis.
    covariates : pd.DataFrame | None
        Optional covariate matrix (rows=samples, indexed like `phenotype`)
        for the null model.
    weight_col : str | None
        Column in `variant_table` to use as per-variant weights. If None,
        no weights are used and the result's `weight` column is `"none"`.
    group_col : str
        Column in `variant_table` defining testable units (e.g. gene).
    out_type : str
        "D" for dichotomous (case/control), "C" for continuous phenotype.
    method : str
        SKAT-O method, e.g. "optimal.adj" (default) or "davies".
    min_variants : int
        Groups with fewer variants than this are skipped.
    backend : SkatBackend | None
        Backend used to run the R null model / SKAT calls. Defaults to
        `RpyBackend()` (requires rpy2 + R + SKAT installed), only
        instantiated here so importing this module never requires rpy2.

    Returns
    -------
    pd.DataFrame with columns:
        feature_id, n_variants, n_samples, p_value, p_value_burden,
        p_value_skat, q_value, weight
    Sorted by p_value ascending.
    """
    if backend is None:
        backend = RpyBackend()

    inputs = build_skat_o_inputs(
        variant_table,
        genotypes,
        phenotype,
        weight_col=weight_col,
        group_col=group_col,
        min_variants=min_variants,
    )

    sample_order = [s for s in phenotype.index if s in genotypes.columns]
    y = phenotype.loc[sample_order].to_numpy()
    X = covariates.loc[sample_order] if covariates is not None else None
    null_model = backend.run_null_model(y, X, out_type)

    rows = []
    for group_input in inputs:
        skat_result = backend.run_skat(
            group_input.Z, null_model, group_input.weights, method
        )
        rows.append(
            {
                "feature_id": group_input.group,
                "n_variants": len(group_input.variant_ids),
                "n_samples": group_input.Z.shape[0],
                "p_value": skat_result["p_value"],
                "p_value_burden": skat_result["p_value_burden"],
                "p_value_skat": skat_result["p_value_skat"],
                "weight": weight_col if weight_col is not None else "none",
            }
        )

    result_df = pd.DataFrame(rows)
    if len(result_df) > 0:
        result_df["q_value"] = bh_fdr(result_df["p_value"])
        result_df = result_df.sort_values("p_value").reset_index(drop=True)
    return result_df


def run_skat_o_from_feather(
    variant_table_feather: str | Path,
    genotypes_feather: str | Path,
    phenotype_table: str | Path,
    sample_col: str = "sample_id",
    phenotype_col: str = "phenotype",
    output_path: str | Path | None = None,
    weight_col: str | None = None,
    group_col: str = "group",
    out_type: str = "D",
    method: str = "optimal.adj",
    min_variants: int = 2,
    backend: "SkatBackend | None" = None,
) -> pd.DataFrame:
    """
    Convenience wrapper that loads a variant-weights table, a genotype
    matrix, and a phenotype table from disk, then runs `run_skat_o`.

    Parameters
    ----------
    variant_table_feather : str | Path
        Feather file with `variant_id`, `group_col`, and (if `weight_col`
        given) a numeric weight column. Matches the schema produced by
        `build_variant_weights_table` in `src/post_prediction/aggregation.py`.
    genotypes_feather : str | Path
        Feather file with `variant_id` plus one integer column per sample.
        Matches the schema produced by `build_genotype_matrix` in
        `src/post_prediction/aggregation.py`.
    phenotype_table : str | Path
        Feather (`.feather`) or TSV (any other extension) file with a
        sample-id column (`sample_col`) and a phenotype column
        (`phenotype_col`, values 1/0).
    sample_col, phenotype_col : str
        Column names in `phenotype_table` identifying sample id and
        phenotype value respectively.
    output_path : str | Path | None
        If given, the result DataFrame is written here as feather.
    weight_col, group_col, out_type, method, min_variants, backend :
        Passed through to `run_skat_o`.

    Returns
    -------
    pd.DataFrame — see `run_skat_o`.
    """
    variant_table = pd.read_feather(variant_table_feather)
    genotypes = pd.read_feather(genotypes_feather)

    phenotype_path = Path(phenotype_table)
    if phenotype_path.suffix == ".feather":
        pheno_df = pd.read_feather(phenotype_path)
    else:
        pheno_df = pd.read_csv(phenotype_path, sep="\t")

    phenotype = pheno_df.set_index(sample_col)[phenotype_col]

    results = run_skat_o(
        variant_table,
        genotypes,
        phenotype,
        weight_col=weight_col,
        group_col=group_col,
        out_type=out_type,
        method=method,
        min_variants=min_variants,
        backend=backend,
    )

    if output_path is not None:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        results.to_feather(out)

    return results
