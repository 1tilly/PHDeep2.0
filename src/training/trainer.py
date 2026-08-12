"""
Model trainer for DeepSEA-family multi-label classifiers.

Features
--------
- Train / validation loop with configurable epochs and batch size.
- Per-epoch mean AUROC and AUPR across all genomic features.
- Early stopping on validation loss.
- Checkpoint saving (best model) and loading.
- Deterministic seed for reproducibility.
"""
from __future__ import annotations

import logging
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import average_precision_score, roc_auc_score
from torch.utils.data import DataLoader

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Loss / optimiser factories
# ---------------------------------------------------------------------------

def get_criterion() -> nn.BCELoss:
    """Binary cross-entropy loss for multi-label classification."""
    return nn.BCELoss()


def get_optimizer(
    model: nn.Module,
    lr: float = 0.01,
    weight_decay: float = 1e-6,
    momentum: float = 0.9,
) -> torch.optim.SGD:
    """SGD with momentum, matching the original DeepSEA training protocol."""
    return torch.optim.SGD(
        model.parameters(),
        lr=lr,
        weight_decay=weight_decay,
        momentum=momentum,
    )


# ---------------------------------------------------------------------------
# Evaluation helpers
# ---------------------------------------------------------------------------

@dataclass
class EpochMetrics:
    loss: float
    mean_auroc: float
    mean_aupr: float
    per_feature_auroc: np.ndarray = field(repr=False)
    per_feature_aupr: np.ndarray = field(repr=False)

    def __str__(self) -> str:
        return (
            f"loss={self.loss:.4f}  "
            f"AUROC={self.mean_auroc:.4f}  "
            f"AUPR={self.mean_aupr:.4f}"
        )


def compute_metrics(labels: np.ndarray, preds: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute per-feature AUROC and AUPR.

    Features with only one class present in `labels` receive NaN.

    Parameters
    ----------
    labels : np.ndarray, shape (N, n_targets)
    preds  : np.ndarray, shape (N, n_targets)

    Returns
    -------
    auroc : np.ndarray, shape (n_targets,)
    aupr  : np.ndarray, shape (n_targets,)
    """
    n_targets = labels.shape[1]
    auroc = np.full(n_targets, np.nan)
    aupr = np.full(n_targets, np.nan)
    for i in range(n_targets):
        if labels[:, i].sum() == 0 or labels[:, i].sum() == len(labels):
            continue
        try:
            auroc[i] = roc_auc_score(labels[:, i], preds[:, i])
            aupr[i] = average_precision_score(labels[:, i], preds[:, i])
        except ValueError:
            pass
    return auroc, aupr


# ---------------------------------------------------------------------------
# Trainer
# ---------------------------------------------------------------------------

class Trainer:
    """
    Manages the training loop for a DeepSEA-family model.

    Parameters
    ----------
    model : nn.Module
    optimizer : torch.optim.Optimizer
    criterion : nn.Module
    device : torch.device
    checkpoint_dir : str | Path
    patience : int
        Number of epochs without validation-loss improvement before stopping.
    """

    def __init__(
        self,
        model: nn.Module,
        optimizer: torch.optim.Optimizer,
        criterion: nn.Module,
        device: torch.device,
        checkpoint_dir: str | Path,
        patience: int = 10,
    ) -> None:
        self.model = model.to(device)
        self.optimizer = optimizer
        self.criterion = criterion
        self.device = device
        self.checkpoint_dir = Path(checkpoint_dir)
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        self.patience = patience
        self._best_val_loss = float("inf")
        self._epochs_no_improve = 0

    # ------------------------------------------------------------------
    # Single epoch
    # ------------------------------------------------------------------

    def train_epoch(self, loader: DataLoader) -> float:
        """Run one training epoch and return mean batch loss."""
        self.model.train()
        total_loss = 0.0
        for x, y in loader:
            x = x.to(self.device)
            y = y.to(self.device)
            self.optimizer.zero_grad()
            pred = self.model(x)
            loss = self.criterion(pred, y)
            loss.backward()
            self.optimizer.step()
            total_loss += loss.item()
        return total_loss / len(loader)

    @torch.no_grad()
    def validate(self, loader: DataLoader) -> EpochMetrics:
        """Evaluate on a validation DataLoader and return EpochMetrics."""
        self.model.eval()
        total_loss = 0.0
        all_labels: list[np.ndarray] = []
        all_preds: list[np.ndarray] = []

        for x, y in loader:
            x = x.to(self.device)
            y = y.to(self.device)
            pred = self.model(x)
            total_loss += self.criterion(pred, y).item()
            all_labels.append(y.cpu().numpy())
            all_preds.append(pred.cpu().numpy())

        labels = np.concatenate(all_labels, axis=0)
        preds = np.concatenate(all_preds, axis=0)
        auroc, aupr = compute_metrics(labels, preds)

        return EpochMetrics(
            loss=total_loss / len(loader),
            mean_auroc=float(np.nanmean(auroc)),
            mean_aupr=float(np.nanmean(aupr)),
            per_feature_auroc=auroc,
            per_feature_aupr=aupr,
        )

    # ------------------------------------------------------------------
    # Full training run
    # ------------------------------------------------------------------

    def fit(
        self,
        train_loader: DataLoader,
        val_loader: DataLoader,
        n_epochs: int,
    ) -> list[dict]:
        """
        Train for up to `n_epochs` with early stopping.

        Returns
        -------
        history : list[dict]
            One entry per epoch with keys train_loss, val_loss, val_auroc, val_aupr.
        """
        history = []
        for epoch in range(1, n_epochs + 1):
            train_loss = self.train_epoch(train_loader)
            val_metrics = self.validate(val_loader)

            logger.info(
                "Epoch %d/%d — train_loss=%.4f  val: %s",
                epoch, n_epochs, train_loss, val_metrics,
            )
            history.append({
                "epoch": epoch,
                "train_loss": train_loss,
                "val_loss": val_metrics.loss,
                "val_auroc": val_metrics.mean_auroc,
                "val_aupr": val_metrics.mean_aupr,
            })

            if val_metrics.loss < self._best_val_loss:
                self._best_val_loss = val_metrics.loss
                self._epochs_no_improve = 0
                self.save_checkpoint(epoch, val_metrics)
                logger.info("  ✓ New best — checkpoint saved.")
            else:
                self._epochs_no_improve += 1
                if self._epochs_no_improve >= self.patience:
                    logger.info("Early stopping after %d epochs without improvement.", epoch)
                    break

        return history

    # ------------------------------------------------------------------
    # Checkpointing
    # ------------------------------------------------------------------

    def save_checkpoint(self, epoch: int, metrics: EpochMetrics) -> Path:
        path = self.checkpoint_dir / f"best_model_epoch{epoch:04d}_val{metrics.loss:.4f}.pt"
        torch.save(
            {
                "epoch": epoch,
                "model_state_dict": self.model.state_dict(),
                "optimizer_state_dict": self.optimizer.state_dict(),
                "val_loss": metrics.loss,
                "val_auroc": metrics.mean_auroc,
                "val_aupr": metrics.mean_aupr,
            },
            path,
        )
        return path

    def load_checkpoint(self, path: str | Path) -> dict:
        """Load a checkpoint into the model and optimizer. Returns the checkpoint dict."""
        ckpt = torch.load(path, map_location=self.device)
        self.model.load_state_dict(ckpt["model_state_dict"])
        self.optimizer.load_state_dict(ckpt["optimizer_state_dict"])
        logger.info("Loaded checkpoint from %s (epoch %d)", path, ckpt.get("epoch", "?"))
        return ckpt


# ---------------------------------------------------------------------------
# Seed helper
# ---------------------------------------------------------------------------

def set_seed(seed: int) -> None:
    """Set seeds for Python, NumPy and PyTorch for deterministic training."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
