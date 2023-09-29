import torch.nn as nn
from base_model import AbstractCNN

class BlainvilleDeepSEA(AbstractCNN):
    def __init__(self, sequence_length, n_targets):
        super(BlainvilleDeepSEA, self).__init__(sequence_length, n_targets)
        
        # 3 Conv1d 1 MaxPool+BatchNorm, 
        # 3 Conv1d 1 MaxPool+BatchNorm+Dropout, 
        # 2 Conv1d 1 BatchNorm+Dropout
        self.conv_net = nn.Sequential(
            nn.Conv1d(4, 320, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(320, 320, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(320, 320, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(kernel_size=self.pool_kernel_size, stride=self.pool_kernel_size),
            nn.BatchNorm1d(320),
            nn.Conv1d(320, 480, kernel_size=self.conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(480, 480, kernel_size=self.conv_kernel_size),
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
            nn.Dropout(p=0.2)
        )
        
        # Calculate the number of channels after the convolutional layers
        layers = [
            # First set of conv layers and pooling
            {'padding': 0, 'dilation': 1, 'kernel_size': self.conv_kernel_size, 'stride': 1},
            {'padding': 0, 'dilation': 1, 'kernel_size': self.conv_kernel_size, 'stride': 1},
            {'padding': 0, 'dilation': 1, 'kernel_size': self.conv_kernel_size, 'stride': 1},
            {'padding': 0, 'dilation': 1, 'kernel_size': self.pool_kernel_size, 'stride': self.pool_kernel_size},
            # Second set of conv layers and pooling
            {'padding': 0, 'dilation': 1, 'kernel_size': self.conv_kernel_size, 'stride': 1},
            {'padding': 0, 'dilation': 1, 'kernel_size': self.conv_kernel_size, 'stride': 1},
            {'padding': 0, 'dilation': 1, 'kernel_size': self.conv_kernel_size, 'stride': 1},
            {'padding': 0, 'dilation': 1, 'kernel_size': self.pool_kernel_size, 'stride': self.pool_kernel_size},
            # Remaining conv layers
            {'padding': 0, 'dilation': 1, 'kernel_size': self.conv_kernel_size, 'stride': 1},
            {'padding': 0, 'dilation': 1, 'kernel_size': self.conv_kernel_size, 'stride': 1},
        ]

        self._n_channels = self.calculate_output_channels(self.sequence_length, layers)
        
        self.classifier = nn.Sequential(
            nn.Linear(960 * self._n_channels, self.n_targets),
            nn.ReLU(inplace=True),
            nn.BatchNorm1d(self.n_targets),
            nn.Linear(self.n_targets, self.n_targets),
            nn.Sigmoid()
        )
        
    def forward(self, x):
        out = self.conv_net(x)
        reshape_out = out.view(out.size(0), 960 * self._n_channels)
        predict = self.classifier(reshape_out)
        return predict