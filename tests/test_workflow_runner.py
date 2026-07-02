from pathlib import Path

import pytest

pytest.importorskip("pybedtools", reason="pybedtools not installed")

from config.pipeline_config import PipelineConfig
from src.workflow.runners import LocalRunner, get_runner, run_pipeline


def test_run_pipeline_local_bed_to_training(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    fixture_dir = repo_root / "tests" / "data" / "encode_dnase_fixture"

    out_bed = tmp_path / "out" / "training_regions.bed"
    out_error = tmp_path / "out" / "read_errors.txt"
    out_features = tmp_path / "out" / "features.txt"

    config = PipelineConfig.from_dict(
        {
            "version": "1",
            "execution": {"backend": "local", "profile": "test"},
            "stage_order": ["bed_to_training"],
            "bed_to_training": {
                "metadata_tsv": str(fixture_dir / "metadata.tsv"),
                "input_bed_dir": str(fixture_dir / "input_bed_files"),
                "output_bed": str(out_bed),
                "output_features": str(out_features),
                "output_errors": str(out_error),
                "assembly": "GRCh38",
                "organism": "Homo sapiens",
            },
        }
    )

    results = run_pipeline(config)

    assert "bed_to_training" in results
    assert out_bed.exists()
    assert out_bed.stat().st_size > 0
    assert out_features.exists()
    assert out_error.exists()
    assert results["bed_to_training"]["written_features"] > 0


def test_pipeline_config_rejects_unknown_stage():
    with pytest.raises(ValueError, match="Unknown stages"):
        PipelineConfig.from_dict(
            {
                "version": "1",
                "execution": {"backend": "local"},
                "stage_order": ["nonexistent_stage"],
            }
        )


def test_get_runner_rejects_unimplemented_backend():
    with pytest.raises(NotImplementedError, match="slurm"):
        get_runner("slurm")
