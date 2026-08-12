"""
PyTorch Dataset for multi-label genomic sequence classification.

Loads genomic windows from a labelled BED file and a reference FASTA,
one-hot encodes each sequence, and returns tensors suitable for training
DeepSEA-family models.

Expected BED format (tab-separated, no header):
    chrom  start  end  feature_name

The feature list must be provided so label vectors are reproducibly ordered.
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pyfaidx
import torch
from torch.utils.data import Dataset


_ENCODE = {
    "A": 1, "a": 1,
    "C": 2, "c": 2,
    "G": 3, "g": 3,
    "T": 4, "t": 4,
    # N / unknown → all-zero row (index 0 in the lookup table)
}
# One-hot lookup: index 0 = [0,0,0,0] (N), 1=A, 2=C, 3=G, 4=T
_OH_TABLE = np.vstack([np.zeros(4, dtype=np.float32), np.eye(4, dtype=np.float32)])


def one_hot_encode(sequence: str) -> np.ndarray:
    """Return a (4, L) float32 array for a nucleotide string of length L."""
    indices = np.array([_ENCODE.get(nt, 0) for nt in sequence], dtype=np.int8)
    # shape: (L, 4) → transpose to (4, L) as expected by Conv1d
    return _OH_TABLE[indices].T


class GenomicsDataset(Dataset):
    """
    Maps labelled BED regions to one-hot sequences and multi-label targets.

    Parameters
    ----------
    bed_path : str | Path
        Path to a tab-separated BED file with columns:
        chrom, start, end, feature_name.
    genome_fasta : str | Path
        Path to a bgzipped or plain FASTA reference genome.
    feature_list : Sequence[str]
        Ordered list of feature names that defines the label vector.
        Features not in this list are ignored.
    seq_len : int
        Fixed sequence length to extract. Windows shorter than seq_len
        are right-padded with Ns; longer windows are centre-cropped.
    """

    def __init__(
        self,
        bed_path: str | Path,
        genome_fasta: str | Path,
        feature_list: Sequence[str],
        seq_len: int = 1000,
    ) -> None:
        self.seq_len = seq_len
        self.feature_index = {f: i for i, f in enumerate(feature_list)}
        self.n_targets = len(feature_list)

        self._genome = pyfaidx.Fasta(str(genome_fasta), as_raw=True)
        self._regions, self._labels = self._load_bed(Path(bed_path))

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _load_bed(self, bed_path: Path):
        """
        Parse the BED file and aggregate feature labels per unique region.

        Returns
        -------
        regions : list[tuple[str, int, int]]
            Sorted list of (chrom, start, end) tuples.
        labels : np.ndarray, shape (n_regions, n_targets), dtype float32
        """
        region_features: dict[tuple[str, int, int], list[str]] = {}

        with bed_path.open() as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                chrom, start, end, feature = parts[0], int(parts[1]), int(parts[2]), parts[3]
                key = (chrom, start, end)
                region_features.setdefault(key, []).append(feature)

        regions = sorted(region_features.keys())
        labels = np.zeros((len(regions), self.n_targets), dtype=np.float32)
        for i, region in enumerate(regions):
            for feat in region_features[region]:
                if feat in self.feature_index:
                    labels[i, self.feature_index[feat]] = 1.0

        return regions, labels

    def _fetch_sequence(self, chrom: str, start: int, end: int) -> str:
        """
        Fetch a sequence from the FASTA, normalising to seq_len via
        centre-crop or right-pad with N.
        """
        # pyfaidx uses 0-based half-open coordinates
        chrom_key = chrom if chrom in self._genome else chrom.lstrip("chr")
        seq = str(self._genome[chrom_key][start:end])

        length = len(seq)
        if length == self.seq_len:
            return seq
        if length > self.seq_len:
            # centre crop
            trim = (length - self.seq_len) // 2
            return seq[trim: trim + self.seq_len]
        # right-pad
        return seq + "N" * (self.seq_len - length)

    # ------------------------------------------------------------------
    # Dataset interface
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return len(self._regions)

    def __getitem__(self, idx: int):
        chrom, start, end = self._regions[idx]
        seq = self._fetch_sequence(chrom, start, end)
        x = torch.from_numpy(one_hot_encode(seq))          # (4, seq_len)
        y = torch.from_numpy(self._labels[idx])             # (n_targets,)
        return x, y
