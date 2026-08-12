"""
BottleNoseDeepSEA architecture.

A wider, deeper variant of DeepSEA that progressively widens channels
(160 → 320 → 720 → 1080) with an additional "bottleneck" expansion in
each block.

Based on the PhDeep deepsea_bottlenose model. Bug fixed: the original
forward pass incorrectly used 960 * n_channels; this version uses the
correct final channel count (1080).
"""
import torch.nn as nn

from src.models.base_model.architecture import AbstractCNN


class BottleNoseDeepSEA(AbstractCNN):
    """
    BottleNoseDeepSEA: wider DeepSEA variant with 1080 output channels.

    Parameters
    ----------
    sequence_length : int
        Length of the one-hot encoded input sequence.
    n_targets : int
        Number of genomic features (output classes).
    """

    _FINAL_CHANNELS = 1080

    def __init__(self, sequence_length: int, n_targets: int) -> None:
        super().__init__(sequence_length, n_targets)

        self.conv_net = nn.Sequential(
            # Block 1: 160 channels (warmup)
            nn.Conv1d(4, 160, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(160, 160, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),

            # Block 2: expand to 320 (Blainville block)
            nn.Conv1d(160, 320, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(320, 320, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(320, 320, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(kernel_size=self.pool_kernel_size, stride=self.pool_kernel_size),
            nn.BatchNorm1d(320),

            # Block 3: expand to 720 (bottlenose expansion)
            nn.Conv1d(320, 480, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(480, 480, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(480, 720, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(720, 720, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(720, 720, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(kernel_size=self.pool_kernel_size, stride=self.pool_kernel_size),
            nn.BatchNorm1d(720),
            nn.Dropout(p=0.2),

            # Block 4: expand to 1080
            nn.Conv1d(720, 1080, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(1080, 1080, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(1080, 1080, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.BatchNorm1d(1080),
            nn.Dropout(p=0.2),
        )

        layer_specs = [
            # Block 1
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            # Block 2 (3 conv + pool)
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.pool_kernel_size, "stride": self.pool_kernel_size},
            # Block 3 (5 conv + pool)
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.pool_kernel_size, "stride": self.pool_kernel_size},
            # Block 4 (3 conv)
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
            {"padding": 0, "dilation": 1, "kernel_size": self.conv_kernel_size, "stride": 1},
        ]
        self._n_channels = self.calculate_output_channels(sequence_length, layer_specs)

        self.classifier = nn.Sequential(
            nn.Linear(self._FINAL_CHANNELS * self._n_channels, n_targets),
            nn.ReLU(inplace=True),
            nn.BatchNorm1d(n_targets),
            nn.Linear(n_targets, n_targets),
            nn.Sigmoid(),
        )

    def forward(self, x):
        out = self.conv_net(x)
        out = out.view(out.size(0), self._FINAL_CHANNELS * self._n_channels)
        return self.classifier(out)
