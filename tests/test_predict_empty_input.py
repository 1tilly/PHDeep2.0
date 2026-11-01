from __future__ import annotations

from pathlib import Path
from typing import TypeVar

import pandas as pd
import pytest

torch = pytest.importorskip("torch")
nn = torch.nn

from src.prediction.predict import ReferencePredictor, VariantEffectPredictor  # noqa: E402

T = TypeVar("T", ReferencePredictor, VariantEffectPredictor)


def _make(cls: type[T]) -> T:
    # Bypass __init__ so no FASTA is needed; the genome is never touched for empty input.
    p = cls.__new__(cls)
    p.model = nn.Identity()
    p.seq_len = 10
    p.device = torch.device("cpu")
    p.batch_size = 4
    return p


def test_predict_bed_empty_file(tmp_path: Path) -> None:
    bed = tmp_path / "empty.bed"
    bed.write_text("# only a comment\n\n")
    df = _make(ReferencePredictor).predict_bed(bed)
    assert len(df) == 0
    assert list(df.columns) == ["chrom", "start", "end"]


def test_predict_regions_empty() -> None:
    preds = _make(ReferencePredictor).predict_regions([])
    assert preds.shape[0] == 0


def test_predict_variants_empty_df() -> None:
    var_df = pd.DataFrame(columns=["chromosome", "start", "reference", "alternate"])
    result = _make(VariantEffectPredictor).predict_variants(var_df)
    assert len(result) == 0
    assert list(result.columns) == ["chromosome", "start", "reference", "alternate"]
