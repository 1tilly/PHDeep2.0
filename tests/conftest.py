"""Shared test fixtures for PHDeep2.0.

Real genomic data is downloaded once per test session via public REST APIs:
  - FASTA : GRCh38 chr22:20 Mb–20.1 Mb (Ensembl REST, plain text)
  - BED   : Ensembl Regulatory Build features in that region
            (enhancer / CTCF_binding_site / promoter)

BED coordinates are stored relative to the start of the downloaded FASTA
slice (i.e. position 0 in the FASTA == chr22:20,000,000 in GRCh38) so
that pyfaidx can resolve them correctly against the local file.
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest
import requests

# 100 kb window — gene-dense region of chr22
_CHROM        = "chr22"
_REGION_START = 20_000_000   # 0-based (BED / pyfaidx convention)
_REGION_END   = 20_100_000   # exclusive

_ENSEMBL      = "https://rest.ensembl.org"
_TIMEOUT      = 60           # seconds


# ── Real FASTA ───────────────────────────────────────────────────────────────

@pytest.fixture(scope="session")
def real_fasta(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """
    GRCh38 chr22:20,000,000–20,100,000 from Ensembl REST API.

    Written as a single-sequence FASTA whose position 0 corresponds to
    chr22:20,000,000 in the reference genome.
    """
    dest = tmp_path_factory.mktemp("genome") / "hg38_chr22_20m.fa"
    # Ensembl uses 1-based, inclusive coords
    url = (
        f"{_ENSEMBL}/sequence/region/human/"
        f"22:{_REGION_START + 1}..{_REGION_END}"
        "?content-type=text/plain"
    )
    r = requests.get(url, timeout=_TIMEOUT)
    r.raise_for_status()
    seq = r.text.strip().upper()
    # Wrap at 60 chars so pyfaidx can build a valid .fai index
    wrapped = "\n".join(textwrap.wrap(seq, 60))
    dest.write_text(f">{_CHROM}\n{wrapped}\n")
    return dest


# ── Real training BED (Ensembl Regulatory Build) ─────────────────────────────

@pytest.fixture(scope="session")
def real_training_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """
    Directory with training_regions.bed + features.txt from Ensembl
    Regulatory Build on chr22:20M–20.1M.

    BED coordinates are offset by _REGION_START so they index correctly
    into the real_fasta fixture file.
    """
    d = tmp_path_factory.mktemp("training")

    url = (
        f"{_ENSEMBL}/overlap/region/human/"
        f"22:{_REGION_START + 1}..{_REGION_END}"
        "?feature=regulatory"
    )
    r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=_TIMEOUT)
    r.raise_for_status()

    rows: list[tuple[str, int, int, str]] = []
    for feat in r.json():
        label = feat.get("description", "")
        if not label:
            continue
        # Convert Ensembl 1-based → 0-based, then subtract region start
        start = int(feat["start"]) - 1 - _REGION_START
        end   = int(feat["end"])       - _REGION_START
        if start < 0 or end > (_REGION_END - _REGION_START) or start >= end:
            continue
        rows.append((_CHROM, start, end, label))

    if len(rows) < 5:
        pytest.skip(
            f"Only {len(rows)} Ensembl regulatory features found — "
            "check network connectivity or widen the region."
        )

    # Preserve insertion order for reproducible feature indices
    feature_types = list(dict.fromkeys(f for _, _, _, f in rows))
    (d / "features.txt").write_text("\n".join(feature_types) + "\n")
    bed_lines = [f"{c}\t{s}\t{e}\t{f}" for c, s, e, f in rows]
    (d / "training_regions.bed").write_text("\n".join(bed_lines) + "\n")
    return d


@pytest.fixture(scope="session")
def feature_list(real_training_dir: Path) -> list[str]:
    """Ordered feature labels from the real training directory."""
    return [
        line.strip()
        for line in (real_training_dir / "features.txt").read_text().splitlines()
        if line.strip()
    ]
