"""Tests that literal-typed config fields are validated at runtime."""

from __future__ import annotations

from typing import Any

import pytest

from config.pipeline_config import PipelineConfig


def _predict(mode: str) -> dict[str, Any]:
    return {
        "mode": mode,
        "model_name": "m",
        "model_checkpoint": "c",
        "sequence_length": 1,
        "n_targets": 1,
        "reference_fasta": "f",
        "output_predictions": "o",
        "input_variants_feather": "v.feather",
    }


def test_invalid_predict_mode_rejected() -> None:
    with pytest.raises(ValueError, match="predict.mode"):
        PipelineConfig.from_dict(
            {"stage_order": ["predict"], "predict": _predict("Variant")}
        )


def test_valid_predict_mode_accepted() -> None:
    cfg = PipelineConfig.from_dict(
        {"stage_order": ["predict"], "predict": _predict("variant")}
    )
    assert cfg.predict is not None
    assert cfg.predict.mode == "variant"


def test_invalid_stats_method_rejected() -> None:
    with pytest.raises(ValueError, match="stats.method"):
        PipelineConfig.from_dict(
            {"stage_order": ["stats"], "stats": {"method": "skat-o"}}
        )


def test_stats_method_none_accepted() -> None:
    cfg = PipelineConfig.from_dict({"stage_order": ["stats"], "stats": {"method": "none"}})
    assert cfg.stats is not None
