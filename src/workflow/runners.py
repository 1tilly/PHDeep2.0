"""Backend-neutral pipeline runners."""

from __future__ import annotations

from typing import Any, Protocol

from config.pipeline_config import PipelineConfig
from src.data_processing.bed_to_training import run_bed_to_training


PipelineRunSummary = dict[str, dict[str, Any]]


class BackendRunner(Protocol):
    backend_name: str

    def run(self, config: PipelineConfig) -> PipelineRunSummary:
        """Run configured stages and return per-stage summaries."""


class LocalRunner:
    """In-process runner for local execution."""

    backend_name = "local"

    def run(self, config: PipelineConfig) -> PipelineRunSummary:
        results: PipelineRunSummary = {}
        for stage in config.stage_order:
            if stage == "bed_to_training":
                if config.bed_to_training is None:
                    raise ValueError("Missing bed_to_training config.")
                bed_cfg = config.bed_to_training
                results[stage] = run_bed_to_training(
                    meta_path=bed_cfg.metadata_tsv,
                    in_bed_path=bed_cfg.input_bed_dir,
                    out_bed_path=bed_cfg.output_bed,
                    error_log_path=bed_cfg.output_errors,
                    feature_out_path=bed_cfg.output_features,
                    assembly=bed_cfg.assembly,
                    organism=bed_cfg.organism,
                )
                continue

            raise NotImplementedError(
                f"Stage '{stage}' is not implemented for local backend."
            )
        return results


def get_runner(backend: str) -> BackendRunner:
    if backend == "local":
        return LocalRunner()
    if backend in {"slurm", "aws_batch"}:
        raise NotImplementedError(
            f"Backend '{backend}' is declared but adapter is not implemented yet."
        )
    raise ValueError(f"Unsupported backend: {backend}")


def run_pipeline(config: PipelineConfig) -> PipelineRunSummary:
    runner = get_runner(config.execution.backend)
    return runner.run(config)
