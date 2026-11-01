import logging
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

pytest.importorskip("pybedtools")

from src.data_loading.bed_loader import BEDParser  # noqa: E402


def test_handlers_removed_after_calls(tmp_path: Path) -> None:
    parser = BEDParser()
    meta_df = pd.DataFrame({"Experiment accession": []})
    parser.get_peaks_per_feature(meta_df, [], {}, str(tmp_path))
    parser.get_peaks_per_feature(meta_df, [], {}, str(tmp_path))
    assert len(logging.getLogger("read_error_log").handlers) == 0


def test_error_logged_once_on_second_call(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    parser = BEDParser()

    def raise_fnf(*args: Any, **kwargs: Any) -> None:
        raise FileNotFoundError

    monkeypatch.setattr(parser, "read_bed_file", raise_fnf)
    meta_df = pd.DataFrame({"Experiment accession": ["EXP1"], "File accession": ["F1"]})
    feature = "cell|DNase|None"
    parser.get_peaks_per_feature(meta_df, [feature], {feature: ["EXP1"]}, str(tmp_path))
    parser.get_peaks_per_feature(meta_df, [feature], {feature: ["EXP1"]}, str(tmp_path))
    content = (tmp_path / "peak_size_read_errors.txt").read_text()
    assert content.count(f"FNF: F1:{feature}") == 1
    assert len(logging.getLogger("read_error_log").handlers) == 0
