import pytest

from src.models.base_model.architecture import AbstractCNN


def _layer(kernel_size: int, stride: int = 1, padding: int = 0, dilation: int = 1) -> dict:
    return {
        'padding': padding,
        'dilation': dilation,
        'kernel_size': kernel_size,
        'stride': stride,
    }


def test_negative_length_raises() -> None:
    with pytest.raises(ValueError, match="after layer 0"):
        AbstractCNN.calculate_output_channels(4, [_layer(8)])


def test_zero_length_raises() -> None:
    with pytest.raises(ValueError):
        AbstractCNN.calculate_output_channels(7, [_layer(8)])


def test_failure_reports_layer_index() -> None:
    with pytest.raises(ValueError, match="after layer 1"):
        AbstractCNN.calculate_output_channels(10, [_layer(8), _layer(4)])


def test_valid_configuration_returns_length() -> None:
    assert AbstractCNN.calculate_output_channels(100, [_layer(8), _layer(4, stride=4)]) == 23
