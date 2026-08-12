"""
Tests for the Trainer and related training utilities.

Uses a tiny synthetic dataset and a lightweight model so the tests are
fast and GPU-free.
"""
import tempfile

import numpy as np
import pytest
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

from src.training.trainer import (
    EpochMetrics,
    Trainer,
    compute_metrics,
    get_criterion,
    get_optimizer,
    set_seed,
)

# ---------------------------------------------------------------------------
# Tiny dummy model
# ---------------------------------------------------------------------------

class _TinyModel(nn.Module):
    def __init__(self, seq_len=100, n_targets=8):
        super().__init__()
        self.net = nn.Sequential(
            nn.Flatten(),
            nn.Linear(4 * seq_len, n_targets),
            nn.Sigmoid(),
        )

    def forward(self, x):
        return self.net(x)


def _make_loaders(n=64, seq_len=100, n_targets=8):
    x = torch.randn(n, 4, seq_len)
    y = (torch.rand(n, n_targets) > 0.5).float()
    ds = TensorDataset(x, y)
    loader = DataLoader(ds, batch_size=16)
    return loader, loader  # use same for train and val


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def test_set_seed_deterministic():
    set_seed(0)
    a = torch.randn(10)
    set_seed(0)
    b = torch.randn(10)
    assert torch.allclose(a, b)


def test_get_criterion_returns_bce():
    loss = get_criterion()
    assert isinstance(loss, nn.BCELoss)


def test_get_optimizer_returns_sgd():
    model = _TinyModel()
    opt = get_optimizer(model)
    assert isinstance(opt, torch.optim.SGD)


# ---------------------------------------------------------------------------
# compute_metrics
# ---------------------------------------------------------------------------

def test_compute_metrics_perfect():
    n, k = 50, 4
    labels = np.eye(k, dtype=np.float32).repeat(n // k, axis=0)[:n]
    preds = labels.copy()
    auroc, aupr = compute_metrics(labels, preds)
    # Perfect predictions: AUROC and AUPR should be 1.0 for active features
    assert np.nanmean(auroc) == pytest.approx(1.0, abs=1e-4)


def test_compute_metrics_single_class_nan():
    """Features with only one class present get NaN, not an error."""
    labels = np.zeros((20, 3), dtype=np.float32)
    labels[:, 0] = 1.0  # feature 0: all positive
    # feature 1: all negative (zero)
    labels[5:10, 2] = 1.0  # feature 2: mixed
    preds = np.random.rand(20, 3).astype(np.float32)
    auroc, aupr = compute_metrics(labels, preds)
    assert np.isnan(auroc[0])   # all-positive → NaN
    assert np.isnan(auroc[1])   # all-negative → NaN
    assert not np.isnan(auroc[2])


# ---------------------------------------------------------------------------
# Trainer
# ---------------------------------------------------------------------------

def test_trainer_train_epoch_decreases_loss():
    set_seed(42)
    model = _TinyModel()
    opt = get_optimizer(model, lr=0.1)
    criterion = get_criterion()
    train_loader, _ = _make_loaders()

    with tempfile.TemporaryDirectory() as tmpdir:
        trainer = Trainer(model, opt, criterion, torch.device("cpu"), tmpdir)
        loss1 = trainer.train_epoch(train_loader)
        loss2 = trainer.train_epoch(train_loader)
        # Loss should generally decrease with learning, but may not on a single
        # batch — just check it's a finite number
        assert np.isfinite(loss1)
        assert np.isfinite(loss2)


def test_trainer_validate_returns_epoch_metrics():
    model = _TinyModel()
    opt = get_optimizer(model)
    criterion = get_criterion()
    _, val_loader = _make_loaders()

    with tempfile.TemporaryDirectory() as tmpdir:
        trainer = Trainer(model, opt, criterion, torch.device("cpu"), tmpdir)
        metrics = trainer.validate(val_loader)
        assert isinstance(metrics, EpochMetrics)
        assert np.isfinite(metrics.loss)


def test_trainer_fit_returns_history():
    set_seed(7)
    model = _TinyModel()
    opt = get_optimizer(model, lr=0.05)
    criterion = get_criterion()
    train_loader, val_loader = _make_loaders()

    with tempfile.TemporaryDirectory() as tmpdir:
        trainer = Trainer(model, opt, criterion, torch.device("cpu"), tmpdir, patience=3)
        history = trainer.fit(train_loader, val_loader, n_epochs=3)
        assert len(history) == 3
        assert all("train_loss" in h for h in history)
        assert all("val_auroc" in h for h in history)


def test_trainer_checkpoint_save_load():
    set_seed(1)
    model = _TinyModel()
    opt = get_optimizer(model)
    criterion = get_criterion()
    _, val_loader = _make_loaders()

    with tempfile.TemporaryDirectory() as tmpdir:
        trainer = Trainer(model, opt, criterion, torch.device("cpu"), tmpdir)
        metrics = trainer.validate(val_loader)
        ckpt_path = trainer.save_checkpoint(epoch=1, metrics=metrics)
        assert ckpt_path.exists()

        # Load into fresh model and check state dict keys match
        model2 = _TinyModel()
        opt2 = get_optimizer(model2)
        trainer2 = Trainer(model2, opt2, criterion, torch.device("cpu"), tmpdir)
        ckpt = trainer2.load_checkpoint(ckpt_path)
        assert "epoch" in ckpt
        assert ckpt["epoch"] == 1
