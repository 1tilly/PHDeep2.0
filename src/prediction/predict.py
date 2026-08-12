"""
Prediction pipelines for trained DeepSEA-family models.

ReferencePredictor
    Slides a fixed-length window across a BED file of regions and
    produces model predictions for each window against the reference
    genome.

VariantEffectPredictor
    For each variant in a variant dataframe (produced by BCFParser),
    computes reference and alternate predictions and returns delta
    scores (alt - ref) representing predicted functional impact.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

from src.data_loading.genomics_dataset import one_hot_encode

logger = logging.getLogger(__name__)


def _load_model(model: nn.Module, checkpoint_path: str | Path, device: torch.device) -> nn.Module:
    """Load model weights from a checkpoint file."""
    ckpt = torch.load(checkpoint_path, map_location=device)
    state = ckpt.get("model_state_dict", ckpt)
    model.load_state_dict(state)
    model.to(device)
    model.eval()
    return model


class ReferencePredictor:
    """
    Predict genomic feature activity across reference genome windows.

    Parameters
    ----------
    model : nn.Module
        Trained DeepSEA-family model (expects input shape [B, 4, seq_len]).
    genome_fasta : str | Path
        Path to reference genome FASTA.
    seq_len : int
        Sequence window length expected by the model.
    device : torch.device | None
        Defaults to CUDA if available, else CPU.
    batch_size : int
        Number of windows to predict at once.
    checkpoint_path : str | Path | None
        If provided, load model weights from this checkpoint before predicting.
    """

    def __init__(
        self,
        model: nn.Module,
        genome_fasta: str | Path,
        seq_len: int = 1000,
        device: torch.device | None = None,
        batch_size: int = 256,
        checkpoint_path: str | Path | None = None,
    ) -> None:
        self.model = model
        self.seq_len = seq_len
        self.device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.batch_size = batch_size
        if checkpoint_path is not None:
            self.model = _load_model(self.model, checkpoint_path, self.device)
        else:
            self.model.to(self.device).eval()

        import pyfaidx
        self._genome = pyfaidx.Fasta(str(genome_fasta), as_raw=True)

    def _fetch(self, chrom: str, start: int, end: int) -> str:
        chrom_key = chrom if chrom in self._genome else chrom.lstrip("chr")
        seq = str(self._genome[chrom_key][start:end])
        length = len(seq)
        if length == self.seq_len:
            return seq
        if length > self.seq_len:
            trim = (length - self.seq_len) // 2
            return seq[trim: trim + self.seq_len]
        return seq + "N" * (self.seq_len - length)

    @torch.no_grad()
    def predict_regions(
        self, regions: Sequence[tuple[str, int, int]]
    ) -> np.ndarray:
        """
        Predict over a list of (chrom, start, end) tuples.

        Returns
        -------
        predictions : np.ndarray, shape (n_regions, n_targets)
        """
        encoded = [one_hot_encode(self._fetch(c, s, e)) for c, s, e in regions]
        results = []
        for i in range(0, len(encoded), self.batch_size):
            batch = torch.from_numpy(np.stack(encoded[i: i + self.batch_size])).to(self.device)
            results.append(self.model(batch).cpu().numpy())
        return np.concatenate(results, axis=0)

    def predict_bed(self, bed_path: str | Path) -> pd.DataFrame:
        """
        Predict for every region in a BED file.

        Parameters
        ----------
        bed_path : str | Path
            Tab-separated BED file (chrom, start, end, ...).

        Returns
        -------
        pd.DataFrame with columns: chrom, start, end, pred_0 … pred_N-1
        """
        records = []
        with Path(bed_path).open() as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                records.append((parts[0], int(parts[1]), int(parts[2])))

        regions = [(c, s, e) for c, s, e in records]
        preds = self.predict_regions(regions)

        df = pd.DataFrame(records, columns=["chrom", "start", "end"])
        pred_cols = {f"pred_{i}": preds[:, i] for i in range(preds.shape[1])}
        return pd.concat([df, pd.DataFrame(pred_cols)], axis=1)


class VariantEffectPredictor:
    """
    Compute variant effect scores (Δ = alt_pred − ref_pred) for a set of
    variants.

    Uses SequenceParser to fetch reference sequences from a local FASTA and
    VariantParser to introduce each alternate allele, then runs both through
    the model.

    Parameters
    ----------
    model : nn.Module
        Trained DeepSEA-family model.
    genome_fasta : str | Path
        Path to reference genome FASTA.
    seq_len : int
        Sequence window length expected by the model.
    device : torch.device | None
    batch_size : int
    checkpoint_path : str | Path | None
        If provided, load model weights from this checkpoint before predicting.
    """

    def __init__(
        self,
        model: nn.Module,
        genome_fasta: str | Path,
        seq_len: int = 1000,
        device: torch.device | None = None,
        batch_size: int = 256,
        checkpoint_path: str | Path | None = None,
    ) -> None:
        self.model = model
        self.seq_len = seq_len
        self.device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.batch_size = batch_size
        if checkpoint_path is not None:
            self.model = _load_model(self.model, checkpoint_path, self.device)
        else:
            self.model.to(self.device).eval()

        import pyfaidx
        self._genome = pyfaidx.Fasta(str(genome_fasta), as_raw=True)

    def _fetch_region(self, chrom: str, start: int, end: int) -> str:
        chrom_key = chrom if chrom in self._genome else chrom.lstrip("chr")
        return str(self._genome[chrom_key][int(start): int(end)])

    @torch.no_grad()
    def _predict_batch(self, sequences: list[str]) -> np.ndarray:
        encoded = [one_hot_encode(s) for s in sequences]
        results = []
        for i in range(0, len(encoded), self.batch_size):
            batch = torch.from_numpy(np.stack(encoded[i: i + self.batch_size])).to(self.device)
            results.append(self.model(batch).cpu().numpy())
        return np.concatenate(results, axis=0)

    def predict_variants(self, var_df: pd.DataFrame, frame_length: int = 500) -> pd.DataFrame:
        """
        Compute ref and alt predictions for every row in `var_df`.

        Expected columns in var_df: chromosome, start, reference, alternate.

        Returns
        -------
        pd.DataFrame with original variant metadata plus:
            ref_pred_<i>, alt_pred_<i>, delta_<i> for each output feature i.
        """
        from src.data_processing.vcf_processing import VariantParser

        ref_seqs: list[str] = []
        alt_seqs: list[str] = []

        for _, var in var_df.iterrows():
            chrom = str(var["chromosome"])
            pos = int(var["start"])
            ref_allele = str(var["reference"])
            alt_allele = str(var["alternate"])

            # Fetch a padded reference window centred on the variant
            flank = self.seq_len // 2
            region_start = max(0, pos - flank)
            region_end = pos + flank + len(ref_allele)
            ref_region = self._fetch_region(chrom, region_start, region_end)

            alt_seq, ref_window = VariantParser.find_variant_in_reference(
                (chrom, pos, ref_allele, alt_allele),
                ref_region,
                region_start,
                self.seq_len,
            )
            ref_seqs.append(ref_window)
            alt_seqs.append(alt_seq)

        ref_preds = self._predict_batch(ref_seqs)
        alt_preds = self._predict_batch(alt_seqs)
        delta = alt_preds - ref_preds

        n_features = ref_preds.shape[1]
        result = var_df.copy().reset_index(drop=True)
        for i in range(n_features):
            result[f"ref_pred_{i}"] = ref_preds[:, i]
            result[f"alt_pred_{i}"] = alt_preds[:, i]
            result[f"delta_{i}"] = delta[:, i]

        return result
