from pathlib import Path
from typing import Any

import pandas as pd
import pytest

pytest.importorskip("pybedtools")

from src.data_loading.bed_loader import BEDParser  # noqa: E402


def _fake_read_bed_file(self: Any, path: str, file: Any, assembly: str) -> pd.DataFrame:
    return pd.DataFrame({"chrom": ["chr1"], "start": [10], "end": [20]})


def test_get_peaks_per_feature_combines_multiple_files(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(BEDParser, "read_bed_file", _fake_read_bed_file)
    parser = BEDParser()
    meta_df = pd.DataFrame(
        {
            "Experiment accession": ["E1", "E1"],
            "File accession": ["F1", "F2"],
        }
    )
    feature = "cell|DNase|None"
    feature_exp_dict = {feature: ["E1"]}

    parser.get_peaks_per_feature(
        meta_df, [feature], feature_exp_dict, str(tmp_path), assembly="hg19"
    )

    out = tmp_path / "analyses" / "cell-DNase-None" / "peak_sizes.tsv"
    assert out.exists()
    result = pd.read_csv(out, sep="\t", header=None)
    assert len(result) == 2
    assert list(result[3]) == [10, 10]
