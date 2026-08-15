"""Shared driver code for the golden `aggregate` -> `stats` fixture (PH2-013).

Not a test file (no `test_` prefix) — used by both
`test_golden_aggregate_stats.py` and `scripts/regenerate_golden_fixtures.py`.
Pure pandas, no torch import.
"""
from __future__ import annotations

import contextlib
from pathlib import Path
from typing import Any, Iterator

import numpy as np
import pandas as pd

from config.pipeline_config import (
    AggregateConfig,
    ExecutionConfig,
    PipelineConfig,
    StatsConfig,
)

WEIGHT_COL = "abs_delta_max"
GROUP_COL = "gene_symbol"
MIN_VARIANTS = 2


class DeterministicSkatBackend:
    """Not a statistical test — a pure, input-sensitive stand-in for
    RpyBackend so golden stats output is reproducible without R/rpy2.

    Matches the `SkatBackend` protocol in
    `src/statistical_testing/skat_o_test.py`: `run_null_model(y, X,
    out_type)` returns an opaque null-model object later passed back into
    `run_skat`; `run_skat(Z, null_model, weights, method)` returns a dict
    with `p_value`, `p_value_burden`, `p_value_skat` keys.
    """

    def run_null_model(self, y: Any, X: Any, out_type: str) -> dict:
        return {"sum_y": float(np.sum(np.asarray(y, dtype=float))), "out_type": out_type}

    def run_skat(
        self, Z: np.ndarray, null_model: dict, weights: np.ndarray | None, method: str
    ) -> dict:
        weight_sum = 0.0 if weights is None else float(np.sum(weights))
        base = float(np.sum(Z)) + 10.0 * weight_sum + null_model["sum_y"]
        p = 1.0 / (1.0 + base)
        return {
            "p_value": p,
            "p_value_burden": min(1.0, 1.5 * p),
            "p_value_skat": 0.5 * p,
        }


@contextlib.contextmanager
def deterministic_skat_backend() -> Iterator[None]:
    """Swap `RpyBackend` in `skat_o_test` for `DeterministicSkatBackend` for
    the duration of the context (so `run_skat_o_from_feather`'s default
    `backend=None` path instantiates the deterministic fake instead of
    requiring R/rpy2)."""
    import src.statistical_testing.skat_o_test as skat_mod

    original = skat_mod.RpyBackend
    skat_mod.RpyBackend = DeterministicSkatBackend  # type: ignore[misc,assignment]
    try:
        yield
    finally:
        skat_mod.RpyBackend = original  # type: ignore[misc,assignment]


def load_golden_predictions(fixture_dir: Path, dest_feather: Path) -> Path:
    """Join `predictions.tsv` with `genotypes.tsv` on the variant key and
    write the result as a feather file shaped like a real `predict`
    (mode=variant) output feather with sample genotype columns attached
    (the shape `LocalRunner._run_aggregate` expects on its
    `input_predictions` — see `tests/test_stats_stage_wiring.py`'s
    `_pred_fixture`, which has the same shape by hand).

    `ref_pred_*`/`alt_pred_*`/`delta_*` columns are cast to float32,
    matching predict.py's real output dtype (torch model output).
    """
    predictions = pd.read_csv(fixture_dir / "predictions.tsv", sep="\t", dtype={"start": int})
    genotypes = pd.read_csv(fixture_dir / "genotypes.tsv", sep="\t", dtype=str)

    join_key = (
        predictions["chromosome"].astype(str)
        + ":"
        + predictions["start"].astype(str)
        + ":"
        + predictions["reference"].astype(str)
        + ":"
        + predictions["alternate"].astype(str)
    )
    predictions = predictions.copy()
    predictions["_variant_id"] = join_key

    merged = predictions.merge(genotypes, left_on="_variant_id", right_on="variant_id", how="left")
    if merged["variant_id"].isna().any():
        missing = merged.loc[merged["variant_id"].isna(), "_variant_id"].tolist()
        raise ValueError(f"No genotype row found for variant id(s): {missing}")

    merged = merged.drop(columns=["_variant_id", "variant_id"])

    pred_cols = [c for c in merged.columns if c.startswith(("ref_pred_", "alt_pred_", "delta_"))]
    for col in pred_cols:
        merged[col] = merged[col].astype("float32")

    dest_feather.parent.mkdir(parents=True, exist_ok=True)
    merged.reset_index(drop=True).to_feather(dest_feather)
    return dest_feather


def golden_pipeline_config(fixture_dir: Path, work_dir: Path) -> PipelineConfig:
    """Build a PipelineConfig for stage_order=("aggregate", "stats") pointed
    at the frozen golden fixture, writing intermediate/output artifacts
    under `work_dir`."""
    predictions_feather = work_dir / "predictions.feather"
    load_golden_predictions(fixture_dir, predictions_feather)

    scores_path = work_dir / "aggregate_scores.feather"
    genotypes_path = work_dir / "aggregate_genotypes.feather"
    results_path = work_dir / "stats_results.feather"

    return PipelineConfig(
        version="1",
        execution=ExecutionConfig(backend="local"),
        stage_order=("aggregate", "stats"),
        aggregate=AggregateConfig(
            input_predictions=[predictions_feather],
            output_scores=scores_path,
            group_col=GROUP_COL,
            sample_ids_file=fixture_dir / "sample_ids.txt",
            output_genotypes=genotypes_path,
        ),
        stats=StatsConfig(
            method="skat_o",
            input_scores=scores_path,
            input_genotypes=genotypes_path,
            phenotype_table=fixture_dir / "phenotype.tsv",
            sample_id_col="sample_id",
            phenotype_col="phenotype",
            group_col="group",
            weight_col=WEIGHT_COL,
            min_variants=MIN_VARIANTS,
            output_results=results_path,
        ),
    )


def run_aggregate_stats(fixture_dir: Path, work_dir: Path) -> dict[str, pd.DataFrame]:
    """Run the real `LocalRunner` over the golden fixture's aggregate ->
    stats stages, with the R-backed SKAT-O backend swapped for
    `DeterministicSkatBackend`. Returns the three output tables read back
    from disk."""
    from src.workflow.runners import LocalRunner

    config = golden_pipeline_config(fixture_dir, work_dir)
    with deterministic_skat_backend():
        LocalRunner().run(config)

    assert config.aggregate is not None and config.stats is not None
    return {
        "scores": pd.read_feather(config.aggregate.output_scores),
        "genotypes": pd.read_feather(config.aggregate.output_genotypes),
        "stats": pd.read_feather(config.stats.output_results),
    }


def write_expected_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, lineterminator="\n", na_rep="")


def read_expected_tsv(path: Path, str_cols: list[str] | None = None) -> pd.DataFrame:
    dtype = {c: str for c in (str_cols or [])}
    return pd.read_csv(path, sep="\t", dtype=dtype, keep_default_na=True, na_values=[""])
