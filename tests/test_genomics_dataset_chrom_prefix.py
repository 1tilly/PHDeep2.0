from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("torch")

from src.data_loading.genomics_dataset import GenomicsDataset  # noqa: E402


def _make_dataset(tmp_path: Path, bed_chrom: str) -> GenomicsDataset:
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">1\nAAAA\n>hr1\nCCCC\n")
    bed = tmp_path / "regions.bed"
    bed.write_text(f"{bed_chrom}\t0\t4\tfeatA\n")
    return GenomicsDataset(bed, fasta, ["featA"], seq_len=4)


def test_chr_prefix_removed_literally(tmp_path: Path) -> None:
    ds = _make_dataset(tmp_path, "chrhr1")
    # "chrhr1" minus the literal "chr" prefix is "hr1", not "1"
    assert ds._fetch_sequence("chrhr1", 0, 4) == "CCCC"


def test_plain_chr_prefix_still_resolves(tmp_path: Path) -> None:
    ds = _make_dataset(tmp_path, "chr1")
    assert ds._fetch_sequence("chr1", 0, 4) == "AAAA"


def test_double_chr_prefix_not_collapsed(tmp_path: Path) -> None:
    ds = _make_dataset(tmp_path, "chrchr1")
    # "chrchr1" -> "chr1", which is absent from the FASTA; must not silently map to "1"
    with pytest.raises(KeyError):
        ds._fetch_sequence("chrchr1", 0, 4)
