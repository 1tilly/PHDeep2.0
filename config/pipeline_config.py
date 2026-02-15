"""Typed pipeline configuration and stage contracts."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Literal


BED_TO_TRAINING_OUTPUT_COLUMNS = ("chrom", "start", "end", "feature")
PREDICTION_OUTPUT_COLUMNS = ("chromosome", "start", "end", "prediction")
AGGREGATION_OUTPUT_COLUMNS = ("chromosome", "start", "end", "score")
STATS_OUTPUT_COLUMNS = ("feature_id", "p_value", "q_value", "effect_size")


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


@dataclass(frozen=True)
class PredictConfig:
    mode: Literal["reference", "variant"]
    model_checkpoint: Path
    output_predictions: Path
    input_regions_bed: Path | None = None
    input_variants_feather: Path | None = None


@dataclass(frozen=True)
class AggregateConfig:
    input_predictions: list[Path]
    output_scores: Path


@dataclass(frozen=True)
class StatsConfig:
    method: Literal["skat_o", "fisher", "none"] = "skat_o"
    input_scores: Path | None = None
    output_results: Path | None = None


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
            )

        predict_cfg = None
        if predict_data is not None:
            predict_cfg = PredictConfig(
                mode=predict_data["mode"],
                model_checkpoint=_to_path(predict_data["model_checkpoint"]),
                output_predictions=_to_path(predict_data["output_predictions"]),
                input_regions_bed=_to_path(predict_data["input_regions_bed"]) if predict_data.get("input_regions_bed") else None,
                input_variants_feather=_to_path(predict_data["input_variants_feather"]) if predict_data.get("input_variants_feather") else None,
            )

        aggregate_cfg = None
        if aggregate_data is not None:
            aggregate_cfg = AggregateConfig(
                input_predictions=[_to_path(p) for p in aggregate_data["input_predictions"]],
                output_scores=_to_path(aggregate_data["output_scores"]),
            )

        stats_cfg = None
        if stats_data is not None:
            stats_cfg = StatsConfig(
                method=stats_data.get("method", "skat_o"),
                input_scores=_to_path(stats_data["input_scores"]) if stats_data.get("input_scores") else None,
                output_results=_to_path(stats_data["output_results"]) if stats_data.get("output_results") else None,
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


def load_pipeline_config(config_path: str | Path) -> PipelineConfig:
    return PipelineConfig.from_json(config_path)
