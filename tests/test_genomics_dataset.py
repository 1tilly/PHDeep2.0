"""Tests for GenomicsDataset and one_hot_encode.

Fixtures `real_fasta`, `real_training_dir`, and `feature_list` are defined
in conftest.py and backed by real GRCh38 data (Ensembl REST API).
"""
from __future__ import annotations

import numpy as np
import pytest
import torch

# ── one_hot_encode unit tests (no fixtures needed) ───────────────────────────

def test_one_hot_encode_correctness():
    """one_hot_encode produces correct (4, L) arrays for known sequences."""
    from src.data_loading.genomics_dataset import one_hot_encode

    arr = one_hot_encode("ACGTN")
    assert arr.shape == (4, 5)
    expected = np.array(
        [
            [1, 0, 0, 0, 0],  # A channel
            [0, 1, 0, 0, 0],  # C channel
            [0, 0, 1, 0, 0],  # G channel
            [0, 0, 0, 1, 0],  # T channel
            # N → all-zero column (5th column is [0,0,0,0])
        ],
        dtype=np.float32,
    )
    np.testing.assert_array_equal(arr, expected)


def test_one_hot_encode_lowercase():
    """Lowercase nucleotides encode identically to uppercase."""
    from src.data_loading.genomics_dataset import one_hot_encode

    np.testing.assert_array_equal(one_hot_encode("ACGT"), one_hot_encode("acgt"))


# ── Dataset tests (require real genome fixtures) ─────────────────────────────

def test_dataset_loads(real_fasta, real_training_dir, feature_list):
    from src.data_loading.genomics_dataset import GenomicsDataset

    ds = GenomicsDataset(
        bed_path=real_training_dir / "training_regions.bed",
        genome_fasta=real_fasta,
        feature_list=feature_list,
        seq_len=1000,
    )
    assert len(ds) > 0


def test_dataset_item_shapes(real_fasta, real_training_dir, feature_list):
    from src.data_loading.genomics_dataset import GenomicsDataset

    n_targets = len(feature_list)
    ds = GenomicsDataset(
        bed_path=real_training_dir / "training_regions.bed",
        genome_fasta=real_fasta,
        feature_list=feature_list,
        seq_len=1000,
    )
    x, y = ds[0]
    assert x.shape == (4, 1000)
    assert y.shape == (n_targets,)
    assert x.dtype == torch.float32
    assert y.dtype == torch.float32


def test_dataset_label_assigned(real_fasta, real_training_dir, feature_list):
    """At least one label should be active for each item in the dataset."""
    from src.data_loading.genomics_dataset import GenomicsDataset

    ds = GenomicsDataset(
        bed_path=real_training_dir / "training_regions.bed",
        genome_fasta=real_fasta,
        feature_list=feature_list,
        seq_len=1000,
    )
    for i in range(min(len(ds), 5)):
        _, y = ds[i]
        assert y.sum().item() >= 1.0, f"Item {i} has no active label"


def test_dataset_one_hot_valid(real_fasta, real_training_dir, feature_list):
    """Each sequence position is a unit vector (known base) or zero vector (N)."""
    from src.data_loading.genomics_dataset import GenomicsDataset

    ds = GenomicsDataset(
        bed_path=real_training_dir / "training_regions.bed",
        genome_fasta=real_fasta,
        feature_list=feature_list,
        seq_len=1000,
    )
    x, _ = ds[0]
    col_sums = x.sum(dim=0)
    assert torch.all((col_sums == 0) | (col_sums == 1))


def test_dataset_centre_crop(real_fasta, real_training_dir, feature_list):
    """seq_len shorter than a BED window triggers centre-cropping."""
    from src.data_loading.genomics_dataset import GenomicsDataset

    ds = GenomicsDataset(
        bed_path=real_training_dir / "training_regions.bed",
        genome_fasta=real_fasta,
        feature_list=feature_list,
        seq_len=100,   # most regulatory features are > 100 bp
    )
    # Find an item whose BED window is longer than 100 bp
    bed = (real_training_dir / "training_regions.bed").read_text().splitlines()
    for line in bed:
        parts = line.split("\t")
        window_len = int(parts[2]) - int(parts[1])
        if window_len > 100:
            break
    else:
        pytest.skip("No BED window longer than 100 bp to test cropping")

    # Verify the tensor is always exactly seq_len regardless of window size
    for i in range(min(len(ds), 5)):
        x, _ = ds[i]
        assert x.shape == (4, 100)


def test_dataset_right_pad(real_fasta, real_training_dir, feature_list):
    """seq_len longer than a BED window triggers N-padding (zero columns)."""
    from src.data_loading.genomics_dataset import GenomicsDataset

    ds = GenomicsDataset(
        bed_path=real_training_dir / "training_regions.bed",
        genome_fasta=real_fasta,
        feature_list=feature_list,
        seq_len=5000,   # longer than any regulatory feature
    )
    x, _ = ds[0]
    assert x.shape == (4, 5000)
    # Some suffix columns must be zero-padded N
    bed = (real_training_dir / "training_regions.bed").read_text().splitlines()
    parts = bed[0].split("\t")
    window_len = int(parts[2]) - int(parts[1])
    if window_len < 5000:
        assert x[:, window_len:].sum().item() == 0.0
