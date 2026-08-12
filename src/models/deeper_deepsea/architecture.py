"""
DeeperDeepSEA architecture.

An extension of DeepSEA (Zhou & Troyanskaya, 2015) with two convolutional
layers per block instead of one, plus BatchNorm after each pooling step.

Based on the Selene example model (FunctionLab/selene, commit ef18ab2).
Ported to PHDeep2.0 with AbstractCNN base class and typed output-size
calculation.
"""
import torch.nn as nn

from src.models.base_model.architecture import AbstractCNN


class DeeperDeepSEA(AbstractCNN):
    """
    DeeperDeepSEA: DeepSEA with doubled convolutional depth per block.

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
            nn.Conv1d(320, 320, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(kernel_size=self.pool_kernel_size, stride=self.pool_kernel_size),
            nn.BatchNorm1d(320),

            nn.Conv1d(320, 480, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(480, 480, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(kernel_size=self.pool_kernel_size, stride=self.pool_kernel_size),
            nn.BatchNorm1d(480),
            nn.Dropout(p=0.2),

            nn.Conv1d(480, 960, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(960, 960, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.BatchNorm1d(960),
            nn.Dropout(p=0.2),
        )

        layer_specs = [
            # Block 1: 2 conv + 1 pool
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.pool_kernel_size, "stride": self.pool_kernel_size},
            # Block 2: 2 conv + 1 pool
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.pool_kernel_size, "stride": self.pool_kernel_size},
            # Block 3: 2 conv
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
        ]
        self._n_channels = self.calculate_output_channels(sequence_length, layer_specs)

        self.classifier = nn.Sequential(
            nn.Linear(960 * self._n_channels, n_targets),
            nn.ReLU(inplace=True),
            nn.BatchNorm1d(n_targets),
            nn.Linear(n_targets, n_targets),
            nn.Sigmoid(),
        )

    def forward(self, x):
        out = self.conv_net(x)
        out = out.view(out.size(0), 960 * self._n_channels)
        return self.classifier(out)
