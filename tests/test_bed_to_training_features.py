from pathlib import Path
from typing import Any

import pandas as pd
import pytest

pytest.importorskip("pybedtools")

import src.data_processing.bed_to_training as btt  # noqa: E402

COLUMNS = [
    "Biosample organism",
    "Assay",
    "Experiment accession",
    "File format",
    "Experiment target",
    "File accession",
]


class FakeBEDParser:
    assays = ["TF ChIP-seq"]
    columns_of_interest = COLUMNS

    def clean_meta_df(self, df: pd.DataFrame, assembly: str = "hg19") -> pd.DataFrame:
        return df

    def generate_feature_exp_dict(self, df: pd.DataFrame) -> dict[str, Any]:
        return {"A": ["E1"], "B": None}

    def generate_feature_list(self, feature_exp_dict: dict[str, Any]) -> list[str]:
        return ["A", "B"]

    def read_bed_file(self, in_path: str, row: Any, assembly: str) -> pd.DataFrame:
        return pd.DataFrame({"chrom": ["chr1"], "start": [10], "end": [20]})

    def choose_assay_version(self, df: pd.DataFrame, assembly: str) -> pd.DataFrame:
        return df.copy()


def test_feature_file_lists_only_written_features(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(btt, "BEDParser", FakeBEDParser)
    meta = pd.DataFrame(
        [
            {
                "Biosample organism": "Homo sapiens",
                "Assay": "TF ChIP-seq",
                "Experiment accession": "E1",
                "File format": "bed narrowPeak",
                "Experiment target": "A-human",
                "File accession": "F1",
            }
        ]
    )
    meta_path = tmp_path / "meta.tsv"
    meta.to_csv(meta_path, sep="\t", index=False)
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    feature_out = tmp_path / "out" / "features.txt"

    result = btt.run_bed_to_training(
        meta_path=meta_path,
        in_bed_path=in_dir,
        out_bed_path=tmp_path / "out" / "train.bed",
        error_log_path=tmp_path / "out" / "errors.txt",
        feature_out_path=feature_out,
    )

    assert result["planned_features"] == 2
    assert result["written_features"] == 1
    assert feature_out.read_text(encoding="utf-8").splitlines() == ["A"]
