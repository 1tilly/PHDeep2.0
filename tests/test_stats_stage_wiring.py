"""Integration test for the aggregate -> stats stage wiring in
`src.workflow.runners.LocalRunner`.

No rpy2/R is required: the SKAT-O R backend (`RpyBackend`) is swapped out
for a local fake via monkeypatch, following the same fake-backend pattern
as `tests/test_skat_o_contract.py` (redefined here rather than imported,
since that module is test-only code).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from config.pipeline_config import (
    AggregateConfig,
    ExecutionConfig,
    PipelineConfig,
    StatsConfig,
)
from src.workflow.runners import LocalRunner


class FakeSkatBackend:
    """Records every call it receives; returns deterministic p-values."""

    def __init__(self) -> None:
        self.null_model_calls: list[tuple[list, object, str]] = []
        self.skat_calls: list[tuple[np.ndarray, object, np.ndarray | None, str]] = []

    def run_null_model(self, y, X, out_type):
        self.null_model_calls.append((list(y), X, out_type))
        return "FAKE_NULL_MODEL"

    def run_skat(self, Z, null_model, weights, method):
        self.skat_calls.append(
            (Z.copy(), null_model, None if weights is None else weights.copy(), method)
        )
        # Distinct-but-deterministic p-value per call, for sort-order checks.
        p = 0.1 + 0.1 * len(self.skat_calls)
        return {"p_value": p, "p_value_burden": p, "p_value_skat": p}


def _pred_fixture() -> pd.DataFrame:
    """5 variants, 2 model features, 3 groups (GENE_A x2, GENE_B x2,
    GENE_C x1), 4 samples. Deltas are chosen so abs_delta_max and eis_diff
    are numerically distinct per variant, so a weight-column mixup would be
    caught by an exact-value assertion."""
    return pd.DataFrame(
        {
            "chromosome": ["1", "1", "1", "1", "1"],
            "start": [1000, 1010, 2000, 2010, 3000],
            "reference": ["A", "G", "C", "T", "A"],
            "alternate": ["T", "C", "T", "A", "G"],
            "gene_symbol": ["GENE_A", "GENE_A", "GENE_B", "GENE_B", "GENE_C"],
            "ref_pred_0": [0.1, 0.1, 0.1, 0.1, 0.1],
            "ref_pred_1": [0.1, 0.1, 0.1, 0.1, 0.1],
            "alt_pred_0": [0.2, 0.3, 0.4, 0.5, 0.6],
            "alt_pred_1": [0.6, 0.15, 0.2, 0.15, 0.15],
            "delta_0": [0.1, 0.2, 0.3, 0.4, 0.5],
            "delta_1": [0.5, 0.05, 0.1, 0.05, 0.05],
            "S1": ["0/0", "0/1", "1/1", "0/0", "0/1"],
            "S2": ["0/1", "0/1", "0/0", "1/1", "0/0"],
            "S3": ["1/1", "0/0", "0/1", "0/1", "1/1"],
            "S4": ["0/0", "1/1", "0/1", "0/0", "0/1"],
        }
    )


def test_aggregate_then_stats_round_trip(tmp_path, monkeypatch):
    pred_path = tmp_path / "predictions.feather"
    _pred_fixture().to_feather(pred_path)

    sample_ids_path = tmp_path / "sample_ids.txt"
    sample_ids_path.write_text("S1\nS2\nS3\nS4\n")

    phenotype_path = tmp_path / "phenotype.feather"
    pd.DataFrame({"sample_id": ["S1", "S2", "S3", "S4"], "phenotype": [1, 0, 1, 0]}).to_feather(
        phenotype_path
    )

    weights_path = tmp_path / "weights.feather"
    genotypes_path = tmp_path / "genotypes.feather"
    results_path = tmp_path / "skat_results.feather"

    config = PipelineConfig(
        version="1",
        execution=ExecutionConfig(backend="local"),
        stage_order=("aggregate", "stats"),
        aggregate=AggregateConfig(
            input_predictions=[pred_path],
            output_scores=weights_path,
            group_col="gene_symbol",
            sample_ids_file=sample_ids_path,
            output_genotypes=genotypes_path,
        ),
        stats=StatsConfig(
            method="skat_o",
            input_scores=weights_path,
            input_genotypes=genotypes_path,
            phenotype_table=phenotype_path,
            sample_id_col="sample_id",
            phenotype_col="phenotype",
            group_col="group",
            weight_col="abs_delta_max",
            min_variants=2,
            output_results=results_path,
        ),
    )

    fake_backend = FakeSkatBackend()
    monkeypatch.setattr(
        "src.statistical_testing.skat_o_test.RpyBackend", lambda: fake_backend
    )

    summary = LocalRunner().run(config)

    # ---- files exist ----
    assert weights_path.exists()
    assert genotypes_path.exists()
    assert results_path.exists()

    weights_df = pd.read_feather(weights_path)
    genotypes_df = pd.read_feather(genotypes_path)
    results_df = pd.read_feather(results_path)

    # ---- weights/genotypes share the same variant order ----
    assert list(weights_df["variant_id"]) == list(genotypes_df["variant_id"])
    assert list(genotypes_df["variant_id"]) == [
        "1:1000:A:T", "1:1010:G:C", "1:2000:C:T", "1:2010:T:A", "1:3000:A:G",
    ]
    assert list(genotypes_df.columns) == ["variant_id", "S1", "S2", "S3", "S4"]

    # ---- stats output contract ----
    assert set(results_df.columns) >= {
        "feature_id", "n_variants", "n_samples", "p_value",
        "p_value_burden", "p_value_skat", "q_value", "weight",
    }

    # ---- min_variants threaded: GENE_C (1 variant) excluded, others kept ----
    assert set(results_df["feature_id"]) == {"GENE_A", "GENE_B"}

    # ---- weight_col threaded: result reports the configured column ----
    assert set(results_df["weight"]) == {"abs_delta_max"}

    # ---- weight_col threaded correctly (not e.g. eis_diff, which is
    # numerically different for this fixture) — check the fake backend's
    # recorded per-group weight arrays directly ----
    weight_arrays = [
        None if w is None else sorted(w.tolist()) for _, _, w, _ in fake_backend.skat_calls
    ]
    assert sorted([0.5, 0.2]) in weight_arrays  # GENE_A abs_delta_max
    assert sorted([0.3, 0.4]) in weight_arrays  # GENE_B abs_delta_max
    # eis_diff would have been [0.6, 0.25] / [0.4, 0.45] — assert those did NOT reach SKAT
    assert sorted([0.6, 0.25]) not in weight_arrays
    assert sorted([0.4, 0.45]) not in weight_arrays

    # ---- group_col threaded: null model saw all 4 phenotypes once ----
    assert len(fake_backend.null_model_calls) == 1
    assert fake_backend.null_model_calls[0][0] == [1, 0, 1, 0]
    assert len(fake_backend.skat_calls) == 2  # GENE_A, GENE_B

    # ---- run summary is sane ----
    assert summary["aggregate"]["n_variants"] == 5
    assert summary["aggregate"]["n_groups"] == 3
    assert summary["aggregate"]["n_ungrouped"] == 0
    assert summary["aggregate"]["group_col"] == "gene_symbol"
    assert summary["aggregate"]["n_samples"] == 4
    assert summary["stats"]["n_tests"] == 2
