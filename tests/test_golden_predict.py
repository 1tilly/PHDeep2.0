"""Property-based (not exact-value) regression tests for `predict` against
the frozen golden fixture (PH2-013).

Uses an UNTRAINED, seeded checkpoint (`golden_checkpoint`), so values are
not hand-verifiable the way `test_golden_aggregate_stats.py`'s formula-based
predictions.tsv is — these tests check shape/range/schema invariants
instead, and that predict.py's own ref-allele sanity-check warning never
fires (confirming the frozen fixture's ref alleles are genuinely correct,
not just self-consistent with each other).
"""
from __future__ import annotations

import logging

import pytest

pytest.importorskip("torch")

import pandas as pd  # noqa: E402


def _feature_triple_columns(n_features: int) -> list[str]:
    cols = []
    for i in range(n_features):
        cols += [f"ref_pred_{i}", f"alt_pred_{i}", f"delta_{i}"]
    return cols


def _assert_delta_equals_alt_minus_ref(df: pd.DataFrame, n_features: int) -> None:
    for i in range(n_features):
        actual = df[f"delta_{i}"].to_numpy()
        expected = (df[f"alt_pred_{i}"] - df[f"ref_pred_{i}"]).to_numpy()
        assert (actual == expected).all()


def _assert_preds_in_unit_interval(df: pd.DataFrame, pred_cols: list[str]) -> None:
    for col in pred_cols:
        assert df[col].between(0.0, 1.0).all(), f"{col} out of [0, 1]"


def test_variant_effect_predictor_direct_call(
    golden_checkpoint, golden_genome, golden_variants, golden_feature_list, caplog
):
    import torch

    from src.models.registry import build_model
    from src.prediction.predict import VariantEffectPredictor

    n_features = len(golden_feature_list)
    model = build_model("deepsea", sequence_length=1000, n_targets=n_features)
    predictor = VariantEffectPredictor(
        model=model,
        genome_fasta=golden_genome,
        seq_len=1000,
        device=torch.device("cpu"),
        batch_size=8,
        checkpoint_path=golden_checkpoint,
    )

    with caplog.at_level(logging.WARNING, logger="src.prediction.predict"):
        result = predictor.predict_variants(golden_variants)

    assert len(result) == 18

    base_cols = ["chromosome", "start", "reference", "alternate", "gene_symbol"]
    expected_cols = base_cols + _feature_triple_columns(n_features)
    assert list(result.columns) == expected_cols

    pred_cols = [c for c in result.columns if c.startswith(("ref_pred_", "alt_pred_"))]
    _assert_preds_in_unit_interval(result, pred_cols)
    _assert_delta_equals_alt_minus_ref(result, n_features)

    mismatch_warnings = [
        r for r in caplog.records if "reference allele mismatch" in r.getMessage().lower()
    ]
    assert not mismatch_warnings, (
        "Frozen fixture's ref alleles must genuinely match the FASTA -- "
        f"got: {[r.getMessage() for r in mismatch_warnings]}"
    )


def test_local_runner_variant_mode_end_to_end(
    tmp_path, golden_checkpoint, golden_genome, golden_variants_feather, golden_feature_list
):
    from config.pipeline_config import PipelineConfig
    from src.workflow.runners import LocalRunner

    n_features = len(golden_feature_list)
    out_path = tmp_path / "variant_predictions.feather"

    config = PipelineConfig.from_dict(
        {
            "version": "1",
            "execution": {"backend": "local"},
            "stage_order": ["predict"],
            "predict": {
                "mode": "variant",
                "model_name": "deepsea",
                "model_checkpoint": str(golden_checkpoint),
                "sequence_length": 1000,
                "n_targets": n_features,
                "reference_fasta": str(golden_genome),
                "output_predictions": str(out_path),
                "input_variants_feather": str(golden_variants_feather),
                "batch_size": 8,
            },
        }
    )
    result = LocalRunner().run(config)

    assert result["predict"]["n_variants"] == 18
    assert out_path.exists()

    df = pd.read_feather(out_path)
    assert len(df) == 18

    base_cols = ["chromosome", "start", "reference", "alternate", "gene_symbol"]
    expected_cols = base_cols + _feature_triple_columns(n_features)
    assert list(df.columns) == expected_cols


def test_local_runner_reference_mode_on_frozen_bed(
    tmp_path, golden_checkpoint, golden_genome, golden_training_dir, golden_feature_list
):
    from config.pipeline_config import PipelineConfig
    from src.workflow.runners import LocalRunner

    n_features = len(golden_feature_list)
    bed_path = golden_training_dir / "training_regions.bed"
    out_path = tmp_path / "reference_predictions.feather"

    config = PipelineConfig.from_dict(
        {
            "version": "1",
            "execution": {"backend": "local"},
            "stage_order": ["predict"],
            "predict": {
                "mode": "reference",
                "model_name": "deepsea",
                "model_checkpoint": str(golden_checkpoint),
                "sequence_length": 1000,
                "n_targets": n_features,
                "reference_fasta": str(golden_genome),
                "output_predictions": str(out_path),
                "input_regions_bed": str(bed_path),
                "batch_size": 8,
            },
        }
    )
    LocalRunner().run(config)

    df = pd.read_feather(out_path)
    bed_line_count = sum(
        1 for line in bed_path.read_text().splitlines() if line.strip() and not line.startswith("#")
    )
    assert len(df) == bed_line_count

    pred_cols = [f"pred_{i}" for i in range(n_features)]
    assert all(c in df.columns for c in pred_cols)
    _assert_preds_in_unit_interval(df, pred_cols)
