"""Integrity checks for the frozen golden pipeline fixture
(`tests/data/golden_pipeline/`, see PH2-013).

Pure pandas + pyfaidx — no torch import, so this file (and
`test_golden_aggregate_stats.py`) can be collected and run in an
environment with zero torch installed.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pyfaidx
import pytest

from src.post_prediction.aggregation import _recode_genotype

FIXTURE_DIR = Path(__file__).resolve().parent / "data" / "golden_pipeline"

# Matches the seq_len/flank used when the fixture's variants were placed
# (scripts/build_golden_variants.py) and the flank formula used by
# VariantEffectPredictor._variant_window (src/prediction/predict.py),
# reimplemented here without importing that (torch-dependent) module.
SEQ_LEN = 1000
FLANK = SEQ_LEN // 2


def _variant_window_bounds(pos_1based: int, ref_len: int, flank: int = FLANK) -> tuple[int, int]:
    pos0 = pos_1based - 1
    region_start = max(0, pos0 - flank)
    region_end = pos0 + flank + ref_len
    return region_start, region_end


@pytest.fixture(scope="module")
def genome() -> pyfaidx.Fasta:
    return pyfaidx.Fasta(str(FIXTURE_DIR / "genome.fa"), as_raw=True)


@pytest.fixture(scope="module")
def genome_length(genome: pyfaidx.Fasta) -> int:
    return len(genome["chr22"][:])


@pytest.fixture(scope="module")
def training_bed_rows() -> list[tuple[str, int, int, str]]:
    rows = []
    with (FIXTURE_DIR / "training" / "training_regions.bed").open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            chrom, start, end, feature = line.split("\t")
            rows.append((chrom, int(start), int(end), feature))
    return rows


@pytest.fixture(scope="module")
def features() -> list[str]:
    return [
        line.strip()
        for line in (FIXTURE_DIR / "training" / "features.txt").read_text().splitlines()
        if line.strip()
    ]


@pytest.fixture(scope="module")
def variants_df() -> pd.DataFrame:
    return pd.read_csv(FIXTURE_DIR / "variants.tsv", sep="\t", dtype={"start": int})


@pytest.fixture(scope="module")
def predictions_df() -> pd.DataFrame:
    return pd.read_csv(FIXTURE_DIR / "predictions.tsv", sep="\t", dtype={"start": int})


@pytest.fixture(scope="module")
def genotypes_df() -> pd.DataFrame:
    return pd.read_csv(FIXTURE_DIR / "genotypes.tsv", sep="\t", dtype=str)


@pytest.fixture(scope="module")
def phenotype_df() -> pd.DataFrame:
    return pd.read_csv(FIXTURE_DIR / "phenotype.tsv", sep="\t")


@pytest.fixture(scope="module")
def sample_ids() -> list[str]:
    return [
        line.strip()
        for line in (FIXTURE_DIR / "sample_ids.txt").read_text().splitlines()
        if line.strip()
    ]


# ── genome ────────────────────────────────────────────────────────────────


def test_genome_is_one_100kb_record(genome: pyfaidx.Fasta, genome_length: int) -> None:
    assert list(genome.keys()) == ["chr22"]
    assert genome_length == 100_000


# ── training BED ─────────────────────────────────────────────────────────


def test_bed_rows_in_bounds_and_features_known(
    training_bed_rows: list[tuple[str, int, int, str]],
    features: list[str],
    genome_length: int,
) -> None:
    assert len(training_bed_rows) > 0
    feature_set = set(features)
    for chrom, start, end, feature in training_bed_rows:
        assert chrom == "chr22"
        assert 0 <= start < end <= genome_length
        assert feature in feature_set


# ── variants.tsv ─────────────────────────────────────────────────────────


def test_variant_ref_alleles_match_frozen_fasta(
    variants_df: pd.DataFrame, genome: pyfaidx.Fasta
) -> None:
    seq = str(genome["chr22"][:])
    for _, row in variants_df.iterrows():
        pos0 = int(row["start"]) - 1
        ref_len = len(row["reference"])
        assert seq[pos0 : pos0 + ref_len] == row["reference"]


def test_variant_alternate_differs_from_reference(variants_df: pd.DataFrame) -> None:
    assert (variants_df["alternate"] != variants_df["reference"]).all()


def test_variant_keys_unique(variants_df: pd.DataFrame) -> None:
    keys = list(
        zip(
            variants_df["chromosome"],
            variants_df["start"],
            variants_df["reference"],
            variants_df["alternate"],
        )
    )
    assert len(set(keys)) == len(keys)


def test_variant_windows_fully_in_bounds(
    variants_df: pd.DataFrame, genome_length: int
) -> None:
    for _, row in variants_df.iterrows():
        region_start, region_end = _variant_window_bounds(
            int(row["start"]), len(row["reference"])
        )
        assert region_start >= 0
        assert region_end <= genome_length


# ── genotypes.tsv / phenotype.tsv ───────────────────────────────────────


def test_genotypes_variant_id_order_matches_variants_tsv(
    variants_df: pd.DataFrame, genotypes_df: pd.DataFrame
) -> None:
    expected_ids = [
        f"{c}:{s}:{r}:{a}"
        for c, s, r, a in zip(
            variants_df["chromosome"],
            variants_df["start"],
            variants_df["reference"],
            variants_df["alternate"],
        )
    ]
    assert list(genotypes_df["variant_id"]) == expected_ids


def test_genotypes_sample_columns_match_sample_ids_order(
    genotypes_df: pd.DataFrame, sample_ids: list[str]
) -> None:
    assert list(genotypes_df.columns) == ["variant_id"] + sample_ids


def test_genotype_calls_all_recodable(genotypes_df: pd.DataFrame, sample_ids: list[str]) -> None:
    for sample in sample_ids:
        for value in genotypes_df[sample]:
            _recode_genotype(value)  # raises ValueError if unrecognized


def test_phenotype_samples_are_subset_with_binary_values(
    phenotype_df: pd.DataFrame, genotypes_df: pd.DataFrame
) -> None:
    genotype_samples = set(genotypes_df.columns) - {"variant_id"}
    assert set(phenotype_df["sample_id"]).issubset(genotype_samples)
    assert set(phenotype_df["phenotype"].unique()).issubset({0, 1})


# ── predictions.tsv ──────────────────────────────────────────────────────


def test_predictions_key_columns_match_variants(
    variants_df: pd.DataFrame, predictions_df: pd.DataFrame
) -> None:
    key_cols = ["chromosome", "start", "reference", "alternate", "gene_symbol"]
    pd.testing.assert_frame_equal(
        predictions_df[key_cols].reset_index(drop=True),
        variants_df[key_cols].reset_index(drop=True),
    )


def test_predictions_delta_equals_alt_minus_ref(predictions_df: pd.DataFrame) -> None:
    delta_cols = [c for c in predictions_df.columns if c.startswith("delta_")]
    assert delta_cols
    for col in delta_cols:
        i = col.removeprefix("delta_")
        expected = predictions_df[f"alt_pred_{i}"] - predictions_df[f"ref_pred_{i}"]
        actual = predictions_df[col]
        assert (actual - expected).abs().max() < 1e-9
