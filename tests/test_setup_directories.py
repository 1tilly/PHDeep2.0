import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "setup_directories.sh"


@pytest.mark.skipif(shutil.which("bash") is None, reason="bash not available")
def test_setup_directories_creates_all_files_without_errors(tmp_path: Path) -> None:
    result = subprocess.run(
        ["bash", str(SCRIPT)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.stderr == ""
    assert (tmp_path / "src/models/cnn_model1/architecture.py").is_file()
    assert (tmp_path / "src/models/cnn_model1/train.py").is_file()


@pytest.mark.skipif(shutil.which("bash") is None, reason="bash not available")
def test_setup_directories_does_not_clobber_existing(tmp_path: Path) -> None:
    readme = tmp_path / "README.md"
    readme.write_text("keep")
    subprocess.run(["bash", str(SCRIPT)], cwd=tmp_path, check=True)
    assert readme.read_text() == "keep"
