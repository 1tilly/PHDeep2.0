import os
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

pytest.importorskip("pybedtools")

from src.data_loading.bed_loader import BEDParser  # noqa: E402


def _fake_read_bed_file(self: Any, path: str, file: Any, assembly: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [10, 100, 5],
            "end": [20, 150, 25],
            "score": [1, 2, 3],
        }
    )


def test_rerun_does_not_duplicate_peaks(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(BEDParser, "read_bed_file", _fake_read_bed_file)
    parser = BEDParser()
    meta_df = pd.DataFrame(
        {
            "Experiment accession": ["E1"],
            "File accession": ["F1"],
        }
    )
    feature = "cell|DNase|None"
    feature_exp_dict = {feature: ["E1"]}
    out = str(tmp_path)

    parser.get_peaks_per_feature(meta_df, [feature], feature_exp_dict, out)
    tsv = os.path.join(out, "analyses", "cell-DNase-None", "peak_sizes.tsv")
    first = len(pd.read_csv(tsv, sep="\t", header=None))
    assert first == 3

    parser.get_peaks_per_feature(meta_df, [feature], feature_exp_dict, out)
    second = len(pd.read_csv(tsv, sep="\t", header=None))
    assert second == first
