"""Shared test fixtures for PHDeep2.0.

Genomic data used across the test suite is a FROZEN golden fixture (PH2-013)
under `tests/data/golden_pipeline/`:
  - FASTA : GRCh38 chr22:20 Mb-20.1 Mb (fetched once via
            `scripts/fetch_golden_genome_fixture.py`, committed to the repo)
  - BED   : Ensembl Regulatory Build features in that region
            (enhancer / CTCF_binding_site / promoter), also frozen

BED/variant coordinates are stored relative to the start of the fixture's
FASTA (i.e. position 0 in the FASTA == chr22:20,000,000 in GRCh38) so that
pyfaidx can resolve them correctly against the local file. See
`tests/data/golden_pipeline/README.md` for full provenance and regeneration
instructions. Earlier versions of these fixtures fetched this data live from
Ensembl REST on every test session — that was replaced (PH2-013) because it
made the suite depend on a live external API for every run.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

GOLDEN_DIR = Path(__file__).resolve().parent / "data" / "golden_pipeline"


@pytest.fixture(scope="session")
def golden_dir() -> Path:
    """Root of the frozen golden pipeline fixture."""
    return GOLDEN_DIR


@pytest.fixture(scope="session")
def golden_genome() -> Path:
    """Path to the frozen chr22:20M-20.1M FASTA."""
    return GOLDEN_DIR / "genome.fa"


@pytest.fixture(scope="session")
def golden_training_dir() -> Path:
    """Directory with training_regions.bed + features.txt."""
    return GOLDEN_DIR / "training"


@pytest.fixture(scope="session")
def golden_feature_list(golden_training_dir: Path) -> list[str]:
    """Ordered feature labels from the golden training directory."""
    return [
        line.strip()
        for line in (golden_training_dir / "features.txt").read_text().splitlines()
        if line.strip()
    ]


@pytest.fixture(scope="session")
def golden_variants() -> pd.DataFrame:
    """The 18-variant golden fixture as a DataFrame."""
    return pd.read_csv(GOLDEN_DIR / "variants.tsv", sep="\t", dtype={"start": int})


@pytest.fixture(scope="session")
def golden_variants_feather(tmp_path_factory: pytest.TempPathFactory, golden_variants: pd.DataFrame) -> Path:
    """golden_variants written to feather (PredictConfig.input_variants_feather
    requires feather, not TSV)."""
    dest = tmp_path_factory.mktemp("golden_variants") / "variants.feather"
    golden_variants.reset_index(drop=True).to_feather(dest)
    return dest


@pytest.fixture(scope="session")
def golden_checkpoint(tmp_path_factory: pytest.TempPathFactory, golden_feature_list: list[str]) -> Path:
    """An UNTRAINED, seeded `deepsea` model checkpoint.

    A real training run isn't needed to test predict's I/O plumbing —
    just a deterministic, loadable checkpoint.
    """
    import torch

    from src.models.registry import build_model

    torch.manual_seed(1234)
    model = build_model("deepsea", sequence_length=1000, n_targets=len(golden_feature_list))

    dest = tmp_path_factory.mktemp("golden_checkpoint") / "deepsea_untrained.pt"
    torch.save({"model_state_dict": model.state_dict()}, dest)
    return dest


@pytest.fixture(scope="session")
def golden_predictions_feather(
    tmp_path_factory: pytest.TempPathFactory, golden_dir: Path
) -> Path:
    """predictions.tsv joined with genotypes.tsv, written to feather (see
    `tests/golden_utils.load_golden_predictions`)."""
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from golden_utils import load_golden_predictions

    dest = tmp_path_factory.mktemp("golden_predictions") / "predictions.feather"
    return load_golden_predictions(golden_dir, dest)
