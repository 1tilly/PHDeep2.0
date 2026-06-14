"""
DeepSEA architecture (Zhou & Troyanskaya, 2015).

Original: https://doi.org/10.1038/nmeth.3547
Selene reference implementation: https://github.com/FunctionLab/selene (commit ef18ab2)

Ported to PHDeep2.0 with modern PyTorch conventions and AbstractCNN base class.
"""
import torch.nn as nn

from src.models.base_model.architecture import AbstractCNN


class DeepSEA(AbstractCNN):
    """
    Canonical DeepSEA model.

    Three convolutional blocks (320 → 480 → 960 channels) followed by two
    fully-connected layers and sigmoid output for multi-label classification.

    Parameters
    ----------
    sequence_length : int
        Length of the one-hot encoded input sequence.
    n_targets : int
        Number of genomic features (output classes).
    """

    def __init__(self, sequence_length: int, n_targets: int) -> None:
        super().__init__(sequence_length, n_targets)

        self.conv_net = nn.Sequential(
            nn.Conv1d(4, 320, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(kernel_size=self.pool_kernel_size, stride=self.pool_kernel_size),
            nn.Dropout(p=0.2),

            nn.Conv1d(320, 480, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(kernel_size=self.pool_kernel_size, stride=self.pool_kernel_size),
            nn.Dropout(p=0.2),

            nn.Conv1d(480, 960, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Dropout(p=0.5),
        )

        layer_specs = [
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.pool_kernel_size, "stride": self.pool_kernel_size},
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.pool_kernel_size, "stride": self.pool_kernel_size},
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
        ]
        self._n_channels = self.calculate_output_channels(sequence_length, layer_specs)

        self.classifier = nn.Sequential(
            nn.Linear(960 * self._n_channels, n_targets),
            nn.ReLU(inplace=True),
            nn.Linear(n_targets, n_targets),
            nn.Sigmoid(),
        )

    def forward(self, x):
        out = self.conv_net(x)
        out = out.view(out.size(0), 960 * self._n_channels)
        return self.classifier(out)
