"""Full pipeline integration tests: train → predict (reference).

Uses the frozen PH2-013 golden fixture (`tests/data/golden_pipeline/`, real
GRCh38 chr22 data captured once via `scripts/fetch_golden_genome_fixture.py`
-- no live network call at test time) via conftest.py fixtures. Trains for
2 epochs on real regulatory feature peaks, then runs inference in reference
mode and verifies the feather output.
"""
from __future__ import annotations

from pathlib import Path


def _read_feature_list(training_dir: Path) -> list[str]:
    return [
        line.strip()
        for line in (training_dir / "features.txt").read_text().splitlines()
        if line.strip()
    ]


def test_train_stage(tmp_path, golden_genome, golden_training_dir):
    """Train for 2 epochs on real chr22 data and verify checkpoint is saved."""
    from config.pipeline_config import PipelineConfig
    from src.workflow.runners import LocalRunner

    golden_feature_list = _read_feature_list(golden_training_dir)
    n_targets = len(golden_feature_list)
    out_dir = tmp_path / "artifacts"

    config = PipelineConfig.from_dict(
        {
            "version": "1",
            "execution": {"backend": "local"},
            "stage_order": ["train"],
            "train": {
                "model_name": "deepsea",
                "sequence_length": 1000,
                "n_targets": n_targets,
                "train_regions_bed": str(golden_training_dir / "training_regions.bed"),
                "reference_fasta": str(golden_genome),
                "output_dir": str(out_dir),
                "n_epochs": 2,
                "batch_size": 4,
                "patience": 100,
                "num_workers": 0,
            },
        }
    )
    result = LocalRunner().run(config)

    assert result["train"]["epochs_trained"] == 2
    checkpoint_dir = Path(result["train"]["checkpoint_dir"])
    assert len(list(checkpoint_dir.glob("*.pt"))) > 0


def test_train_then_predict(tmp_path, golden_genome, golden_training_dir):
    """Train → predict (reference mode) on real chr22 data."""
    import pandas as pd

    from config.pipeline_config import PipelineConfig
    from src.workflow.runners import LocalRunner

    golden_feature_list = _read_feature_list(golden_training_dir)
    n_targets = len(golden_feature_list)
    out_dir = tmp_path / "artifacts"
    bed_path = golden_training_dir / "training_regions.bed"

    # ── Train ─────────────────────────────────────────────────────────────────
    train_config = PipelineConfig.from_dict(
        {
            "version": "1",
            "execution": {"backend": "local"},
            "stage_order": ["train"],
            "train": {
                "model_name": "deepsea",
                "sequence_length": 1000,
                "n_targets": n_targets,
                "train_regions_bed": str(bed_path),
                "reference_fasta": str(golden_genome),
                "output_dir": str(out_dir),
                "n_epochs": 2,
                "batch_size": 4,
                "patience": 100,
                "num_workers": 0,
            },
        }
    )
    train_result = LocalRunner().run(train_config)
    checkpoint_dir = Path(train_result["train"]["checkpoint_dir"])
    checkpoints = sorted(checkpoint_dir.glob("*.pt"))
    assert checkpoints, "No checkpoint saved during training"
    best_ckpt = checkpoints[0]

    # ── Predict (reference mode) ───────────────────────────────────────────────
    pred_out = tmp_path / "preds" / "predictions.feather"
    predict_config = PipelineConfig.from_dict(
        {
            "version": "1",
            "execution": {"backend": "local"},
            "stage_order": ["predict"],
            "predict": {
                "mode": "reference",
                "model_name": "deepsea",
                "model_checkpoint": str(best_ckpt),
                "sequence_length": 1000,
                "n_targets": n_targets,
                "reference_fasta": str(golden_genome),
                "output_predictions": str(pred_out),
                "input_regions_bed": str(bed_path),
                "batch_size": 4,
            },
        }
    )
    predict_result = LocalRunner().run(predict_config)

    assert predict_result["predict"]["n_regions"] > 0
    assert pred_out.exists()

    df = pd.read_feather(pred_out)
    assert len(df) > 0
    assert "pred_0" in df.columns
    assert df["pred_0"].between(0.0, 1.0).all(), "Sigmoid output must be in [0, 1]"

    # Now that the region set is frozen and stable, predicted row count is
    # newly pinnable against the BED's line count. `predict_bed` emits one
    # row PER BED LINE (src/prediction/predict.py), unlike `GenomicsDataset`,
    # which dedupes rows sharing the same (chrom, start, end) into a single
    # region (src/data_loading/genomics_dataset.py's `_load_bed`). The
    # golden fixture's training BED happens to have no duplicate regions, so
    # this assertion doesn't exercise that dedup distinction, but it does
    # pin `predict_bed`'s own one-row-per-line contract.
    bed_line_count = sum(
        1 for line in bed_path.read_text().splitlines() if line.strip() and not line.startswith("#")
    )
    assert len(df) == bed_line_count
