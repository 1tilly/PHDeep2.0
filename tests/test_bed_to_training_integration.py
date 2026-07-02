from pathlib import Path
import subprocess
import sys

import pytest

pytest.importorskip("pybedtools", reason="pybedtools not installed")


def test_bed_to_training_runs_on_encode_fixture(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    fixture_dir = repo_root / "tests" / "data" / "encode_dnase_fixture"
    meta = fixture_dir / "metadata.tsv"
    input_bed_dir = fixture_dir / "input_bed_files"

    assert meta.exists(), "Missing fixture metadata.tsv"
    assert input_bed_dir.exists(), "Missing fixture input bed directory"

    out_bed = tmp_path / "training_regions.bed"
    out_error = tmp_path / "read_errors.txt"
    out_features = tmp_path / "features.txt"

    cmd = [
        sys.executable,
        str(repo_root / "src" / "data_processing" / "bed_to_training.py"),
        "-m",
        str(meta),
        "-i",
        str(input_bed_dir),
        "-o",
        str(out_bed),
        "-e",
        str(out_error),
        "-f",
        str(out_features),
        "-a",
        "GRCh38",
    ]

    completed = subprocess.run(cmd, cwd=repo_root, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert out_features.exists()
    assert out_bed.exists()
    assert out_bed.stat().st_size > 0
