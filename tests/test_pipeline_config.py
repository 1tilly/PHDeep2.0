from pathlib import Path

import pytest

from config.pipeline_config import (
    AGGREGATION_OUTPUT_COLUMNS,
    BED_TO_TRAINING_OUTPUT_COLUMNS,
    GENOTYPE_MATRIX_KEY_COLUMN,
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
    assert PREDICTION_OUTPUT_COLUMNS == ("chromosome", "start", "reference", "alternate")
    assert AGGREGATION_OUTPUT_COLUMNS == (
        "chromosome", "start", "end", "reference", "alternate", "variant_id",
        "group", "n_features", "eis_ref", "eis_alt", "eis_diff",
        "abs_delta_max", "abs_delta_sum", "l2_delta",
    )
    assert GENOTYPE_MATRIX_KEY_COLUMN == "variant_id"
    assert STATS_OUTPUT_COLUMNS == (
        "feature_id", "n_variants", "n_samples", "p_value",
        "p_value_burden", "p_value_skat", "weight", "q_value",
    )


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
                    "model_name": "deepsea",
                    "model_checkpoint": "models/model.pt",
                    "sequence_length": 1000,
                    "n_targets": 919,
                    "reference_fasta": "data/hg38.fa",
                    "output_predictions": "results/pred.feather",
                    # input_regions_bed intentionally omitted → should raise
                },
            }
        )


def test_skat_o_requires_genotypes_and_phenotype():
    with pytest.raises(ValueError, match="phenotype_table"):
        PipelineConfig.from_dict(
            {
                "version": "1",
                "execution": {"backend": "local"},
                "stage_order": ["stats"],
                "stats": {
                    "method": "skat_o",
                    # input_scores, input_genotypes, phenotype_table all
                    # intentionally omitted → should raise
                },
            }
        )


def test_aggregate_and_stats_config_round_trip_new_fields():
    config = PipelineConfig.from_dict(
        {
            "version": "1",
            "execution": {"backend": "local"},
            "stage_order": ["aggregate", "stats"],
            "aggregate": {
                "input_predictions": ["preds/a.feather", "preds/b.feather"],
                "output_scores": "out/weights.feather",
                "group_col": "gene_symbol",
                "sample_ids_file": "out/sample_ids.txt",
                "output_genotypes": "out/genotypes.feather",
                "output_group_summary": "out/group_summary.feather",
            },
            "stats": {
                "method": "skat_o",
                "input_scores": "out/weights.feather",
                "input_genotypes": "out/genotypes.feather",
                "phenotype_table": "cohort/phenotype.tsv",
                "sample_id_col": "subject_id",
                "phenotype_col": "case_control",
                "group_col": "gene_symbol",
                "weight_col": "eis_diff",
                "min_variants": 3,
                "output_results": "out/skat_results.feather",
            },
        }
    )

    assert config.aggregate is not None
    ac = config.aggregate
    assert ac.sample_ids_file == Path("out/sample_ids.txt")
    assert ac.output_genotypes == Path("out/genotypes.feather")
    assert ac.output_group_summary == Path("out/group_summary.feather")

    assert config.stats is not None
    sc = config.stats
    assert sc.input_genotypes == Path("out/genotypes.feather")
    assert sc.phenotype_table == Path("cohort/phenotype.tsv")
    assert sc.sample_id_col == "subject_id"
    assert sc.phenotype_col == "case_control"
    assert sc.group_col == "gene_symbol"
    assert sc.weight_col == "eis_diff"
    assert sc.min_variants == 3


def test_aggregate_and_stats_config_new_field_defaults():
    config = PipelineConfig.from_dict(
        {
            "version": "1",
            "execution": {"backend": "local"},
            "stage_order": ["aggregate", "stats"],
            "aggregate": {
                "input_predictions": ["preds/a.feather"],
                "output_scores": "out/weights.feather",
            },
            "stats": {
                "method": "skat_o",
                "input_scores": "out/weights.feather",
                "input_genotypes": "out/genotypes.feather",
                "phenotype_table": "cohort/phenotype.tsv",
            },
        }
    )

    assert config.aggregate is not None
    ac = config.aggregate
    assert ac.sample_ids_file is None
    assert ac.output_genotypes is None
    assert ac.output_group_summary is None

    assert config.stats is not None
    sc = config.stats
    assert sc.sample_id_col == "sample_id"
    assert sc.phenotype_col == "phenotype"
    assert sc.group_col == "group"
    assert sc.weight_col is None
    assert sc.min_variants == 2
