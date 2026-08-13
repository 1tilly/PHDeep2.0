"""Tests for VariantEffectPredictor's variant-window coordinate handling."""
import textwrap
from pathlib import Path

import numpy as np
import pytest
import torch

from src.prediction.predict import ReferencePredictor, VariantEffectPredictor, _variant_window


@pytest.fixture
def tiny_fasta(tmp_path) -> Path:
    # "ACGT" repeating: base at 0-based index i is "ACGT"[i % 4], so a
    # one-base misalignment shows up as a different letter.
    seq = "ACGT" * 500  # 2000 bp
    p = tmp_path / "tiny.fa"
    p.write_text(">chr1\n" + "\n".join(textwrap.wrap(seq, 60)) + "\n")
    return p


@pytest.fixture
def marker_fasta(tmp_path) -> Path:
    # 2000 bp of 'A' with a single 'G' marker at 0-based index 999
    # (i.e. 1-based VCF position 1000).
    seq = list("A" * 2000)
    seq[999] = "G"
    seq = "".join(seq)
    p = tmp_path / "marker.fa"
    p.write_text(">chr1\n" + "\n".join(textwrap.wrap(seq, 60)) + "\n")
    return p


def _make_vep(fasta, seq_len) -> VariantEffectPredictor:
    return VariantEffectPredictor(
        model=torch.nn.Module(), genome_fasta=fasta,
        seq_len=seq_len, device=torch.device("cpu"), batch_size=4,
    )


def _capture_predict_batch(vep, monkeypatch):
    calls = []

    def fake_predict_batch(seqs):
        calls.append(list(seqs))
        return np.zeros((len(seqs), 2), dtype=np.float32)

    monkeypatch.setattr(vep, "_predict_batch", fake_predict_batch)
    return calls


def test_variant_window_aligns_to_true_1based_position(tiny_fasta, monkeypatch):
    vep = _make_vep(tiny_fasta, seq_len=10)
    calls = _capture_predict_batch(vep, monkeypatch)
    var_df = __import__("pandas").DataFrame([{
        "chromosome": "chr1", "start": 21, "reference": "A", "alternate": "T",
    }])

    vep.predict_variants(var_df)

    ref_seqs, alt_seqs = calls[0], calls[1]
    assert ref_seqs[0] == "TACGTACGTA"
    assert alt_seqs[0] == "TACGTTCGTA"
    assert ref_seqs[0][5] == "A"   # base at 1-based position 21 matches VCF `reference`
    assert alt_seqs[0][5] == "T"   # alt spliced at the same offset


def test_alt_replaces_the_marker_base_not_its_neighbour(marker_fasta, monkeypatch):
    vep = _make_vep(marker_fasta, seq_len=100)
    calls = _capture_predict_batch(vep, monkeypatch)
    var_df = __import__("pandas").DataFrame([{
        "chromosome": "chr1", "start": 1000, "reference": "G", "alternate": "C",
    }])

    vep.predict_variants(var_df)

    ref_seqs, alt_seqs = calls[0], calls[1]
    assert ref_seqs[0].count("G") == 1
    assert alt_seqs[0].count("G") == 0
    assert alt_seqs[0].count("C") == 1
    assert alt_seqs[0].index("C") == ref_seqs[0].index("G")


def test_variant_window_arithmetic():
    assert _variant_window(1000, 1, 500) == (999, 499, 1500)
    assert _variant_window(1, 1, 500) == (0, 0, 501)
    assert _variant_window(10, 1, 500) == (9, 0, 510)
    assert _variant_window(1000, 5, 500) == (999, 499, 1504)


def test_variant_window_region_start_never_negative():
    for pos in (1, 2, 10, 100, 499, 500, 501):
        pos0, region_start, region_end = _variant_window(pos, 1, 500)
        assert region_start >= 0
        assert pos0 - region_start >= 0


def test_reference_predictor_fetch_is_0based_unaffected(tiny_fasta):
    rp = ReferencePredictor(
        model=torch.nn.Module(), genome_fasta=tiny_fasta,
        seq_len=4, device=torch.device("cpu"), batch_size=4,
    )
    assert rp._fetch("chr1", 0, 4) == "ACGT"
    assert rp._fetch("chr1", 3, 7) == "TACG"


def test_variant_effect_predictor_fetch_region_is_0based_unaffected(tiny_fasta):
    vep = _make_vep(tiny_fasta, seq_len=4)
    assert vep._fetch_region("chr1", 0, 4) == "ACGT"
    assert vep._fetch_region("chr1", 3, 7) == "TACG"
