from pathlib import Path

import pytest

from config.pipeline_config import (
    AGGREGATION_OUTPUT_COLUMNS,
    BED_TO_TRAINING_OUTPUT_COLUMNS,
    PREDICTION_OUTPUT_COLUMNS,
    STATS_OUTPUT_COLUMNS,
    PipelineConfig,
    load_pipeline_config,
)


def test_load_example_pipeline_config():
    config = load_pipeline_config("config/pipeline.example.json")
    assert config.version == "1"
    assert config.execution.backend == "local"
    assert config.stage_order == ("bed_to_training",)
    assert config.bed_to_training is not None
    assert config.bed_to_training.assembly == "GRCh38"
    assert config.bed_to_training.metadata_tsv == Path("tests/data/encode_dnase_fixture/metadata.tsv")


def test_contract_columns_are_stable():
    assert BED_TO_TRAINING_OUTPUT_COLUMNS == ("chrom", "start", "end", "feature")
    assert PREDICTION_OUTPUT_COLUMNS == ("chromosome", "start", "end", "prediction")
    assert AGGREGATION_OUTPUT_COLUMNS == ("chromosome", "start", "end", "score")
    assert STATS_OUTPUT_COLUMNS == ("feature_id", "p_value", "q_value", "effect_size")


def test_invalid_backend_rejected():
    with pytest.raises(ValueError, match="Unsupported backend"):
        PipelineConfig.from_dict(
            {
                "version": "1",
                "execution": {"backend": "kubernetes"},
                "stage_order": [],
            }
        )


def test_missing_stage_config_rejected():
    with pytest.raises(ValueError, match="Missing stage config"):
        PipelineConfig.from_dict(
            {
                "version": "1",
                "execution": {"backend": "local"},
                "stage_order": ["bed_to_training"],
            }
        )


def test_predict_reference_requires_regions():
    with pytest.raises(ValueError, match="input_regions_bed"):
        PipelineConfig.from_dict(
            {
                "version": "1",
                "execution": {"backend": "local"},
                "stage_order": ["predict"],
                "predict": {
                    "mode": "reference",
                    "model_checkpoint": "models/model.pt",
                    "output_predictions": "results/pred.feather",
                },
            }
        )
