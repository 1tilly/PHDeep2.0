from pathlib import Path
from typing import Any

import pandas as pd
import pytest

pytest.importorskip("pybedtools")

from src.data_processing.bed_to_training import BedTrainingSet  # noqa: E402


class _FakeParser:
    def read_bed_file(self, in_path: str, row: Any, assembly: str) -> pd.DataFrame:
        return pd.DataFrame({"chrom": ["chr1"], "start": [1], "end": [2]})

    def choose_assay_version(self, df: pd.DataFrame, assembly: str) -> pd.DataFrame:
        return df.copy()


def _setup(n: int) -> tuple[pd.DataFrame, list[str], dict[str, list[str]]]:
    assay_df = pd.DataFrame(
        {
            "Experiment accession": [f"EXP{k}" for k in range(n)],
            "File accession": [f"F{k}" for k in range(n)],
        }
    )
    features = [f"feat{k}" for k in range(n)]
    feature_exp_dict = {f"feat{k}": [f"EXP{k}"] for k in range(n)}
    return assay_df, features, feature_exp_dict


def test_default_has_no_limit(tmp_path: Path) -> None:
    assay_df, features, fed = _setup(25)
    builder = BedTrainingSet(bed_parser=_FakeParser())
    written = builder.write_bed_by_feature(
        assay_df, features, fed, str(tmp_path), str(tmp_path / "out.bed"), str(tmp_path / "err.txt")
    )
    assert written == features


def test_limit_counts_files(tmp_path: Path) -> None:
    assay_df, features, fed = _setup(5)
    builder = BedTrainingSet(bed_parser=_FakeParser())
    written = builder.write_bed_by_feature(
        assay_df,
        features,
        fed,
        str(tmp_path),
        str(tmp_path / "out.bed"),
        str(tmp_path / "err.txt"),
        assembly="hg19",
        exp_count=2,
    )
    assert written == features[:2]
