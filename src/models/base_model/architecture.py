import torch
import torch.nn as nn
from abc import ABC, abstractmethod

class AbstractCNN(ABC, nn.Module):
    def __init__(self, sequence_length: int, n_targets: int, conv_kernel_size: int = 8, pool_kernel_size: int = 4):
        super(AbstractCNN, self).__init__()
        self.sequence_length = sequence_length
        self.n_targets = n_targets
        self.conv_kernel_size = conv_kernel_size
        self.pool_kernel_size = pool_kernel_size
    
    @classmethod
    @abstractmethod
    def output_size_1d_pt(cls, seq_len: int, padding: int, dilation: int, kernel_size: int, stride: int) -> int:
        """
        Calculates the output size of a 1D convolutional layer in PyTorch.

        Args:
            seq_len (int): Length of the input sequence.
            padding (int): Padding size.
            dilation (int): Dilation size.
            kernel_size (int): Kernel size.
            stride (int): Stride size.

        Returns:
            int: Output size of the convolutional layer.
        """
        return (seq_len + 2 * padding - dilation * (kernel_size - 1) - 1) // stride + 1

    @classmethod
    def calculate_output_channels(cls, sequence_length:int, layers:list) -> int:
        """
        Calculate the output channels given a list of layer specifications.
        """
        for layer in reversed(layers):
            sequence_length = cls.output_size_1d_pt(
                sequence_length,
                padding=layer['padding'],
                dilation=layer['dilation'],
                kernel_size=layer['kernel_size'],
                stride=layer['stride']
            )
        return sequence_length


    @abstractmethod
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Stub for the forward function.
        """
        pass