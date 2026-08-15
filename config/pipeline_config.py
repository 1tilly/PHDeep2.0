"""Typed pipeline configuration and stage contracts."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

BED_TO_TRAINING_OUTPUT_COLUMNS = ("chrom", "start", "end", "feature")
# `predict` (mode="variant") output has these columns plus one
# ref_pred_<i>/alt_pred_<i>/delta_<i> triple per model output feature — a
# fixed-width tuple can't express "one per model feature", so only the
# fixed columns are captured here; see docs/stage_contracts.md.
PREDICTION_OUTPUT_COLUMNS = ("chromosome", "start", "reference", "alternate")
AGGREGATION_OUTPUT_COLUMNS = (
    "chromosome", "start", "end", "reference", "alternate", "variant_id",
    "group", "n_features", "eis_ref", "eis_alt", "eis_diff",
    "abs_delta_max", "abs_delta_sum", "l2_delta",
)
GENOTYPE_MATRIX_KEY_COLUMN = "variant_id"
STATS_OUTPUT_COLUMNS = (
    "feature_id", "n_variants", "n_samples", "p_value",
    "p_value_burden", "p_value_skat", "weight", "q_value",
)


def _to_path(value: str | Path) -> Path:
    return value if isinstance(value, Path) else Path(value)


@dataclass(frozen=True)
class ExecutionConfig:
    backend: Literal["local", "slurm", "aws_batch"] = "local"
    profile: str = "default"


@dataclass(frozen=True)
class BedToTrainingConfig:
    metadata_tsv: Path
    input_bed_dir: Path
    output_bed: Path
    output_features: Path
    output_errors: Path
    assembly: Literal["hg19", "GRCh38"] = "hg19"
    organism: str = "Homo sapiens"


@dataclass(frozen=True)
class TrainConfig:
    model_name: str
    sequence_length: int
    n_targets: int
    train_regions_bed: Path
    reference_fasta: Path
    output_dir: Path
    # Training hyperparameters
    learning_rate: float = 0.01
    batch_size: int = 128
    n_epochs: int = 100
    patience: int = 10
    seed: int = 42
    val_split: float = 0.1
    num_workers: int = 4


@dataclass(frozen=True)
class PredictConfig:
    mode: Literal["reference", "variant"]
    model_name: str
    model_checkpoint: Path
    sequence_length: int
    n_targets: int
    reference_fasta: Path
    output_predictions: Path
    input_regions_bed: Path | None = None
    input_variants_feather: Path | None = None
    batch_size: int = 256


@dataclass(frozen=True)
class AggregateConfig:
    input_predictions: list[Path]
    output_scores: Path
    group_col: str = "gene_symbol"
    feature_names_file: Path | None = None
    sample_ids_file: Path | None = None
    output_genotypes: Path | None = None
    output_group_summary: Path | None = None


@dataclass(frozen=True)
class StatsConfig:
    method: Literal["skat_o", "fisher", "none"] = "skat_o"
    input_scores: Path | None = None
    output_results: Path | None = None
    input_genotypes: Path | None = None
    phenotype_table: Path | None = None
    sample_id_col: str = "sample_id"
    phenotype_col: str = "phenotype"
    group_col: str = "group"
    weight_col: str | None = None
    min_variants: int = 2


@dataclass(frozen=True)
class PipelineConfig:
    version: str
    execution: ExecutionConfig
    stage_order: tuple[str, ...]
    bed_to_training: BedToTrainingConfig | None = None
    train: TrainConfig | None = None
    predict: PredictConfig | None = None
    aggregate: AggregateConfig | None = None
    stats: StatsConfig | None = None

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "PipelineConfig":
        execution_data = data.get("execution", {})
        execution = ExecutionConfig(
            backend=execution_data.get("backend", "local"),
            profile=execution_data.get("profile", "default"),
        )

        stage_order = tuple(data.get("stage_order", []))
        bed_data = data.get("bed_to_training")
        train_data = data.get("train")
        predict_data = data.get("predict")
        aggregate_data = data.get("aggregate")
        stats_data = data.get("stats")

        bed_cfg = None
        if bed_data is not None:
            bed_cfg = BedToTrainingConfig(
                metadata_tsv=_to_path(bed_data["metadata_tsv"]),
                input_bed_dir=_to_path(bed_data["input_bed_dir"]),
                output_bed=_to_path(bed_data["output_bed"]),
                output_features=_to_path(bed_data["output_features"]),
                output_errors=_to_path(bed_data["output_errors"]),
                assembly=bed_data.get("assembly", "hg19"),
                organism=bed_data.get("organism", "Homo sapiens"),
            )

        train_cfg = None
        if train_data is not None:
            train_cfg = TrainConfig(
                model_name=train_data["model_name"],
                sequence_length=int(train_data["sequence_length"]),
                n_targets=int(train_data["n_targets"]),
                train_regions_bed=_to_path(train_data["train_regions_bed"]),
                reference_fasta=_to_path(train_data["reference_fasta"]),
                output_dir=_to_path(train_data["output_dir"]),
                learning_rate=float(train_data.get("learning_rate", 0.01)),
                batch_size=int(train_data.get("batch_size", 128)),
                n_epochs=int(train_data.get("n_epochs", 100)),
                patience=int(train_data.get("patience", 10)),
                seed=int(train_data.get("seed", 42)),
                val_split=float(train_data.get("val_split", 0.1)),
                num_workers=int(train_data.get("num_workers", 4)),
            )

        predict_cfg = None
        if predict_data is not None:
            predict_cfg = PredictConfig(
                mode=predict_data["mode"],
                model_name=predict_data["model_name"],
                model_checkpoint=_to_path(predict_data["model_checkpoint"]),
                sequence_length=int(predict_data["sequence_length"]),
                n_targets=int(predict_data["n_targets"]),
                reference_fasta=_to_path(predict_data["reference_fasta"]),
                output_predictions=_to_path(predict_data["output_predictions"]),
                input_regions_bed=_to_path(predict_data["input_regions_bed"]) if predict_data.get("input_regions_bed") else None,
                input_variants_feather=_to_path(predict_data["input_variants_feather"]) if predict_data.get("input_variants_feather") else None,
                batch_size=int(predict_data.get("batch_size", 256)),
            )

        aggregate_cfg = None
        if aggregate_data is not None:
            aggregate_cfg = AggregateConfig(
                input_predictions=[_to_path(p) for p in aggregate_data["input_predictions"]],
                output_scores=_to_path(aggregate_data["output_scores"]),
                group_col=aggregate_data.get("group_col", "gene_symbol"),
                feature_names_file=_to_path(aggregate_data["feature_names_file"]) if aggregate_data.get("feature_names_file") else None,
                sample_ids_file=_to_path(aggregate_data["sample_ids_file"]) if aggregate_data.get("sample_ids_file") else None,
                output_genotypes=_to_path(aggregate_data["output_genotypes"]) if aggregate_data.get("output_genotypes") else None,
                output_group_summary=_to_path(aggregate_data["output_group_summary"]) if aggregate_data.get("output_group_summary") else None,
            )

        stats_cfg = None
        if stats_data is not None:
            stats_cfg = StatsConfig(
                method=stats_data.get("method", "skat_o"),
                input_scores=_to_path(stats_data["input_scores"]) if stats_data.get("input_scores") else None,
                output_results=_to_path(stats_data["output_results"]) if stats_data.get("output_results") else None,
                input_genotypes=_to_path(stats_data["input_genotypes"]) if stats_data.get("input_genotypes") else None,
                phenotype_table=_to_path(stats_data["phenotype_table"]) if stats_data.get("phenotype_table") else None,
                sample_id_col=stats_data.get("sample_id_col", "sample_id"),
                phenotype_col=stats_data.get("phenotype_col", "phenotype"),
                group_col=stats_data.get("group_col", "group"),
                weight_col=stats_data.get("weight_col"),
                min_variants=int(stats_data.get("min_variants", 2)),
            )

        config = cls(
            version=str(data.get("version", "1")),
            execution=execution,
            stage_order=stage_order,
            bed_to_training=bed_cfg,
            train=train_cfg,
            predict=predict_cfg,
            aggregate=aggregate_cfg,
            stats=stats_cfg,
        )
        config.validate()
        return config

    @classmethod
    def from_json(cls, config_path: str | Path) -> "PipelineConfig":
        with open(config_path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        return cls.from_dict(payload)

    def validate(self) -> None:
        stage_to_config = {
            "bed_to_training": self.bed_to_training,
            "train": self.train,
            "predict": self.predict,
            "aggregate": self.aggregate,
            "stats": self.stats,
        }
        unknown_stages = set(self.stage_order) - set(stage_to_config.keys())
        if unknown_stages:
            raise ValueError(f"Unknown stages in stage_order: {sorted(unknown_stages)}")

        missing = [stage for stage in self.stage_order if stage_to_config.get(stage) is None]
        if missing:
            raise ValueError(f"Missing stage config for stages in stage_order: {missing}")

        if self.execution.backend not in {"local", "slurm", "aws_batch"}:
            raise ValueError(f"Unsupported backend: {self.execution.backend}")

        if self.predict is not None:
            if self.predict.mode == "reference" and self.predict.input_regions_bed is None:
                raise ValueError("predict.input_regions_bed is required when mode='reference'")
            if self.predict.mode == "variant" and self.predict.input_variants_feather is None:
                raise ValueError("predict.input_variants_feather is required when mode='variant'")
            if self.predict.mode == "variant" and self.predict.input_regions_bed is not None:
                raise ValueError("predict.input_regions_bed should not be set when mode='variant'")

        if "stats" in self.stage_order and self.stats is not None and self.stats.method == "skat_o":
            if (
                self.stats.input_scores is None
                or self.stats.input_genotypes is None
                or self.stats.phenotype_table is None
            ):
                raise ValueError(
                    "stats.input_scores, stats.input_genotypes, and stats.phenotype_table "
                    "are all required when method='skat_o' — these are cohort inputs "
                    "(per-sample genotypes and phenotype) that this pipeline does not "
                    "generate."
                )


def load_pipeline_config(config_path: str | Path) -> PipelineConfig:
    return PipelineConfig.from_json(config_path)
