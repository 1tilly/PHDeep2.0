"""Backend-neutral pipeline runners."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Protocol

from config.pipeline_config import PipelineConfig

logger = logging.getLogger(__name__)

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
            logger.info("Running stage: %s", stage)
            if stage == "bed_to_training":
                results[stage] = self._run_bed_to_training(config)
            elif stage == "train":
                results[stage] = self._run_train(config)
            elif stage == "predict":
                results[stage] = self._run_predict(config)
            elif stage == "aggregate":
                results[stage] = self._run_aggregate(config)
            elif stage == "stats":
                results[stage] = self._run_stats(config)
            else:
                raise NotImplementedError(
                    f"Stage '{stage}' is not implemented for local backend."
                )
        return results

    # ------------------------------------------------------------------
    # Stage implementations
    # ------------------------------------------------------------------

    def _run_bed_to_training(self, config: PipelineConfig) -> dict:
        if config.bed_to_training is None:
            raise ValueError("Missing bed_to_training config.")
        from src.data_processing.bed_to_training import run_bed_to_training
        bed_cfg = config.bed_to_training
        return run_bed_to_training(
            meta_path=bed_cfg.metadata_tsv,
            in_bed_path=bed_cfg.input_bed_dir,
            out_bed_path=bed_cfg.output_bed,
            error_log_path=bed_cfg.output_errors,
            feature_out_path=bed_cfg.output_features,
            assembly=bed_cfg.assembly,
            organism=bed_cfg.organism,
        )

    def _run_train(self, config: PipelineConfig) -> dict:
        if config.train is None:
            raise ValueError("Missing train config.")
        tc = config.train

        # Lazy imports so non-model installs don't fail at import time
        import torch
        from torch.utils.data import DataLoader, random_split

        from src.data_loading.genomics_dataset import GenomicsDataset
        from src.models.registry import build_model
        from src.training.trainer import (
            Trainer,
            get_criterion,
            get_optimizer,
            set_seed,
        )

        set_seed(tc.seed)

        # Load feature list from the features file produced by bed_to_training
        features_file = tc.train_regions_bed.parent / "features.txt"
        if features_file.exists():
            feature_list = [line.strip() for line in features_file.read_text().splitlines() if line.strip()]
        else:
            # Fall back to n_targets integer indices
            feature_list = [str(i) for i in range(tc.n_targets)]

        dataset = GenomicsDataset(
            bed_path=tc.train_regions_bed,
            genome_fasta=tc.reference_fasta,
            feature_list=feature_list,
            seq_len=tc.sequence_length,
        )

        n_val = max(1, int(tc.val_split * len(dataset)))
        n_train = len(dataset) - n_val
        train_ds, val_ds = random_split(dataset, [n_train, n_val])

        train_loader = DataLoader(
            train_ds, batch_size=tc.batch_size, shuffle=True, num_workers=tc.num_workers
        )
        val_loader = DataLoader(
            val_ds, batch_size=tc.batch_size * 2, shuffle=False, num_workers=tc.num_workers
        )

        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model = build_model(tc.model_name, tc.sequence_length, tc.n_targets)
        optimizer = get_optimizer(model, lr=tc.learning_rate)
        criterion = get_criterion()

        checkpoint_dir = Path(tc.output_dir) / "checkpoints"
        trainer = Trainer(model, optimizer, criterion, device, checkpoint_dir, patience=tc.patience)
        history = trainer.fit(train_loader, val_loader, n_epochs=tc.n_epochs)

        best = min(history, key=lambda h: h["val_loss"])
        return {
            "epochs_trained": len(history),
            "best_val_loss": best["val_loss"],
            "best_val_auroc": best["val_auroc"],
            "best_epoch": best["epoch"],
            "checkpoint_dir": str(checkpoint_dir),
        }

    def _run_predict(self, config: PipelineConfig) -> dict:
        if config.predict is None:
            raise ValueError("Missing predict config.")
        pc = config.predict

        import torch

        from src.models.registry import build_model

        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model = build_model(pc.model_name, pc.sequence_length, pc.n_targets)

        out_path = Path(pc.output_predictions)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        if pc.mode == "reference":
            from src.prediction.predict import ReferencePredictor

            if pc.input_regions_bed is None:
                raise ValueError("predict.input_regions_bed is required when mode='reference'")

            ref_predictor = ReferencePredictor(
                model=model,
                genome_fasta=pc.reference_fasta,
                seq_len=pc.sequence_length,
                device=device,
                batch_size=pc.batch_size,
                checkpoint_path=pc.model_checkpoint,
            )
            df = ref_predictor.predict_bed(pc.input_regions_bed)
            df.to_feather(out_path)
            return {"mode": "reference", "n_regions": len(df), "output": str(out_path)}

        else:  # variant
            import pandas as pd

            from src.prediction.predict import VariantEffectPredictor

            if pc.input_variants_feather is None:
                raise ValueError("predict.input_variants_feather is required when mode='variant'")

            var_predictor = VariantEffectPredictor(
                model=model,
                genome_fasta=pc.reference_fasta,
                seq_len=pc.sequence_length,
                device=device,
                batch_size=pc.batch_size,
                checkpoint_path=pc.model_checkpoint,
            )
            var_df = pd.read_feather(pc.input_variants_feather)
            delta_df = var_predictor.predict_variants(var_df)
            delta_df.to_feather(out_path)
            return {"mode": "variant", "n_variants": len(delta_df), "output": str(out_path)}

    def _run_aggregate(self, config: PipelineConfig) -> dict:
        if config.aggregate is None:
            raise ValueError("Missing aggregate config.")
        ac = config.aggregate

        import pandas as pd

        from src.post_prediction.aggregation import (
            aggregate_variant_scores,
            build_genotype_matrix,
            build_variant_weights_table,
        )

        dfs = [pd.read_feather(p) for p in ac.input_predictions]

        # Guard against pd.concat silently NaN-filling mismatched schemas.
        column_sets = [frozenset(df.columns) for df in dfs]
        if len(set(column_sets)) > 1:
            disagreeing = [
                str(p) for p, cols in zip(ac.input_predictions, column_sets)
                if cols != column_sets[0]
            ]
            raise ValueError(
                "Input prediction feathers have mismatched columns; "
                f"the following file(s) disagree with {ac.input_predictions[0]}: "
                f"{', '.join(disagreeing)}"
            )

        combined = pd.concat(dfs, ignore_index=True)

        feature_names = None
        if ac.feature_names_file is not None:
            feature_names = [
                line.strip()
                for line in Path(ac.feature_names_file).read_text().splitlines()
                if line.strip()
            ]

        weights = build_variant_weights_table(combined, group_col=ac.group_col)

        out_path = Path(ac.output_scores)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        weights.to_feather(out_path)

        summary: dict = {
            "n_variants": len(weights),
            "n_groups": weights["group"].dropna().nunique(),
            "n_ungrouped": int(weights["group"].isna().sum()),
            "group_col": ac.group_col,
            "output": str(out_path),
        }

        if ac.sample_ids_file is not None:
            if ac.output_genotypes is None:
                raise ValueError(
                    "aggregate.output_genotypes is required when "
                    "aggregate.sample_ids_file is set."
                )
            sample_ids = [
                line.strip()
                for line in Path(ac.sample_ids_file).read_text().splitlines()
                if line.strip()
            ]
            genotypes = build_genotype_matrix(combined, sample_ids)

            genotypes_path = Path(ac.output_genotypes)
            genotypes_path.parent.mkdir(parents=True, exist_ok=True)
            genotypes.to_feather(genotypes_path)

            summary["n_samples"] = len(sample_ids)
            summary["output_genotypes"] = str(genotypes_path)

        if ac.output_group_summary is not None:
            group_summary = aggregate_variant_scores(
                combined, group_col=ac.group_col, feature_names=feature_names
            )
            group_summary_path = Path(ac.output_group_summary)
            group_summary_path.parent.mkdir(parents=True, exist_ok=True)
            group_summary.to_feather(group_summary_path)
            summary["output_group_summary"] = str(group_summary_path)

        return summary

    def _run_stats(self, config: PipelineConfig) -> dict:
        if config.stats is None:
            raise ValueError("Missing stats config.")
        sc = config.stats

        if sc.method == "none":
            return {"method": "none", "skipped": True}

        if sc.method == "skat_o":
            from src.statistical_testing.skat_o_test import run_skat_o_from_feather

            if sc.input_scores is None or sc.input_genotypes is None or sc.phenotype_table is None:
                raise ValueError(
                    "stats.input_scores, stats.input_genotypes, and "
                    "stats.phenotype_table are all required when method='skat_o'"
                )

            result_df = run_skat_o_from_feather(
                sc.input_scores,
                sc.input_genotypes,
                sc.phenotype_table,
                sample_col=sc.sample_id_col,
                phenotype_col=sc.phenotype_col,
                output_path=sc.output_results,
                weight_col=sc.weight_col,
                group_col=sc.group_col,
                min_variants=sc.min_variants,
            )
            result: dict = {"method": "skat_o", "n_tests": len(result_df)}
            if sc.output_results is not None:
                result["output"] = str(sc.output_results)
            return result

        raise NotImplementedError(f"Stats method '{sc.method}' is not implemented.")


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
