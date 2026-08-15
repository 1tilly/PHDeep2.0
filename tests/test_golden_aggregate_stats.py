"""Exact-value regression test for the `aggregate` -> `stats` stages against
the frozen golden fixture (PH2-013).

Pure pandas — no torch, no R/rpy2 (the R-backed `RpyBackend` is swapped for
`golden_utils.DeterministicSkatBackend`). This is the core deliverable: a
snapshot test against real, versioned expected output, catching exactly the
kind of column-mixup/shape bugs that shipped undetected before this fixture
existed (see docs/2026-08-13-session-handover.md).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest
from golden_utils import WEIGHT_COL, read_expected_tsv, run_aggregate_stats

from config.pipeline_config import AGGREGATION_OUTPUT_COLUMNS, STATS_OUTPUT_COLUMNS

FIXTURE_DIR = Path(__file__).resolve().parent / "data" / "golden_pipeline"
EXPECTED_DIR = FIXTURE_DIR / "expected"


@pytest.fixture(scope="module")
def result(tmp_path_factory: pytest.TempPathFactory) -> dict[str, pd.DataFrame]:
    work_dir = tmp_path_factory.mktemp("golden_aggregate_stats")
    return run_aggregate_stats(FIXTURE_DIR, work_dir)


@pytest.fixture(scope="module")
def expected_scores() -> pd.DataFrame:
    return read_expected_tsv(EXPECTED_DIR / "aggregate_scores.tsv", str_cols=["variant_id", "group"])


@pytest.fixture(scope="module")
def expected_stats() -> pd.DataFrame:
    return read_expected_tsv(EXPECTED_DIR / "stats_results.tsv", str_cols=["feature_id", "weight"])


def test_scores_match_golden_exactly(result, expected_scores) -> None:
    actual = result["scores"]
    pdt.assert_frame_equal(
        actual,
        expected_scores,
        check_dtype=True,
        check_like=False,
        check_exact=False,
        rtol=1e-12,
        atol=1e-15,
    )


def test_stats_match_golden_exactly(result, expected_stats) -> None:
    actual = result["stats"]
    pdt.assert_frame_equal(
        actual,
        expected_stats,
        check_dtype=True,
        check_like=False,
        check_exact=False,
        rtol=1e-12,
        atol=1e-15,
    )


def test_scores_column_order_matches_expected(expected_scores) -> None:
    assert list(expected_scores.columns) == list(AGGREGATION_OUTPUT_COLUMNS)


def test_stats_column_order_matches_expected(expected_stats) -> None:
    assert list(expected_stats.columns) == list(STATS_OUTPUT_COLUMNS)


def test_ungrouped_variants_retained_in_scores_with_null_group(result) -> None:
    scores = result["scores"]
    ungrouped = scores[scores["group"].isna()]
    assert len(ungrouped) == 2


def test_singleton_group_absent_from_stats(result) -> None:
    stats = result["stats"]
    assert "GENE_SINGLE" not in set(stats["feature_id"])


def test_ungrouped_variants_absent_from_stats(result) -> None:
    stats = result["stats"]
    # Ungrouped (null-group) variants never form a testable "group" at all.
    assert set(stats["feature_id"]) == {"GENE_A", "GENE_B", "GENE_C"}


def test_weight_column_reports_configured_column(result) -> None:
    stats = result["stats"]
    assert set(stats["weight"]) == {WEIGHT_COL}


def test_p_values_unique_and_ascending(result) -> None:
    p_values = result["stats"]["p_value"].tolist()
    assert len(set(p_values)) == len(p_values)
    assert p_values == sorted(p_values)


def test_q_value_matches_independent_bh_fdr_recomputation(result) -> None:
    from src.statistical_testing.skat_o_test import bh_fdr

    stats = result["stats"]
    expected_q = bh_fdr(stats["p_value"].to_numpy())
    np.testing.assert_allclose(stats["q_value"].to_numpy(), expected_q, rtol=1e-12, atol=1e-15)
