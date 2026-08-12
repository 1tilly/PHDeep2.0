"""
Forward-pass shape tests for all DeepSEA-family models.

These tests do not require a GPU or a real dataset — they use random
tensors to verify that each model produces the correct output shape.
"""
import pytest
import torch

SEQ_LEN = 1000
N_TARGETS = 64
BATCH = 4


def _random_batch():
    return torch.randn(BATCH, 4, SEQ_LEN)


# ---------------------------------------------------------------------------
# AbstractCNN utility tests
# ---------------------------------------------------------------------------

def test_output_size_1d_no_padding():
    from src.models.base_model.architecture import AbstractCNN
    # Conv1d with kernel=8, stride=1, no padding on input length 1000
    out = AbstractCNN.output_size_1d_pt(1000, padding=0, dilation=1, kernel_size=8, stride=1)
    assert out == 993  # 1000 - 8 + 1 = 993


def test_output_size_1d_maxpool():
    from src.models.base_model.architecture import AbstractCNN
    out = AbstractCNN.output_size_1d_pt(993, padding=0, dilation=1, kernel_size=4, stride=4)
    assert out == 248  # floor((993 - 4) / 4) + 1 = 248


# ---------------------------------------------------------------------------
# DeepSEA
# ---------------------------------------------------------------------------

def test_deepsea_forward_shape():
    from src.models.deepsea.architecture import DeepSEA
    model = DeepSEA(SEQ_LEN, N_TARGETS)
    out = model(_random_batch())
    assert out.shape == (BATCH, N_TARGETS)


def test_deepsea_output_range():
    from src.models.deepsea.architecture import DeepSEA
    model = DeepSEA(SEQ_LEN, N_TARGETS)
    model.eval()
    with torch.no_grad():
        out = model(_random_batch())
    assert out.min() >= 0.0 and out.max() <= 1.0, "Sigmoid output must be in [0, 1]"


# ---------------------------------------------------------------------------
# DeeperDeepSEA
# ---------------------------------------------------------------------------

def test_deeper_deepsea_forward_shape():
    from src.models.deeper_deepsea.architecture import DeeperDeepSEA
    model = DeeperDeepSEA(SEQ_LEN, N_TARGETS)
    out = model(_random_batch())
    assert out.shape == (BATCH, N_TARGETS)


# ---------------------------------------------------------------------------
# BlainvilleDeepSEA (existing model)
# ---------------------------------------------------------------------------

def test_blainville_forward_shape():
    from src.models.blainville.architecture import BlainvilleDeepSEA
    model = BlainvilleDeepSEA(SEQ_LEN, N_TARGETS)
    out = model(_random_batch())
    assert out.shape == (BATCH, N_TARGETS)


# ---------------------------------------------------------------------------
# JellyFishDeepSEA (existing model)
# ---------------------------------------------------------------------------

def test_jellyfish_forward_shape():
    from src.models.jellyfish.architecture import JellyFishDeepSEA
    model = JellyFishDeepSEA(SEQ_LEN, N_TARGETS)
    out = model(_random_batch())
    assert out.shape == (BATCH, N_TARGETS)


def test_jellyfish_custom_dims():
    from src.models.jellyfish.architecture import JellyFishDeepSEA
    model = JellyFishDeepSEA(SEQ_LEN, N_TARGETS, n_convModules=6, layer_dimensions=(128, 512))
    out = model(_random_batch())
    assert out.shape == (BATCH, N_TARGETS)


# ---------------------------------------------------------------------------
# BottleNoseDeepSEA
# ---------------------------------------------------------------------------

def test_bottlenose_forward_shape():
    from src.models.bottlenose.architecture import BottleNoseDeepSEA
    model = BottleNoseDeepSEA(SEQ_LEN, N_TARGETS)
    out = model(_random_batch())
    assert out.shape == (BATCH, N_TARGETS)


def test_bottlenose_output_range():
    from src.models.bottlenose.architecture import BottleNoseDeepSEA
    model = BottleNoseDeepSEA(SEQ_LEN, N_TARGETS)
    model.eval()
    with torch.no_grad():
        out = model(_random_batch())
    assert out.min() >= 0.0 and out.max() <= 1.0
