from functools import partial

import torch.nn as nn

from src.models.base_model.architecture import AbstractCNN


class JellyFishDeepSEA(AbstractCNN):
    def __init__(self, sequence_length, n_targets, dropout=0.2, dropout_freq=6, n_convModules=10, layer_dimensions=None, batch_pool_freq=5):
        super(JellyFishDeepSEA, self).__init__(sequence_length, n_targets)
        if layer_dimensions is None:
            layer_dimensions = (320, 960)

        self.conv_net = self.create_conv_net(dropout, dropout_freq, n_convModules, layer_dimensions, batch_pool_freq)

        self.classifier = nn.Sequential(
            nn.Linear(self._max_dim * self._n_channels, n_targets),
            nn.ReLU(inplace=True),
            nn.BatchNorm1d(n_targets),
            nn.Linear(n_targets, n_targets),
            nn.Sigmoid()
        )

    def forward(self, x):
        out = self.conv_net(x)
        reshape_out = out.view(out.size(0), self._max_dim * self._n_channels)
        predict = self.classifier(reshape_out)
        return predict

    def create_conv_net(self, dropout, dropout_freq, n_convModules, layer_dimensions, batch_pool_freq):
        layer_list = []
        skip_one_layer = False
        conv1d_out = partial(self.output_size_1d_pt, padding=0, dilation=1, kernel_size=self.conv_kernel_size, stride=1)
        maxpool1d_out = partial(
            self.output_size_1d_pt,
            padding=0,
            dilation=1,
            kernel_size=self.pool_kernel_size,
            stride=self.pool_kernel_size,
        )
        channel_value = self.sequence_length
        min_dim = layer_dimensions[0]
        max_dim = layer_dimensions[1]
        self._max_dim = max_dim
        layer_list.append(nn.Conv1d(4, min_dim, kernel_size=self.conv_kernel_size))
        channel_value = conv1d_out(channel_value)

        n_convModules = max(2, int(n_convModules))
        n_convModules -= 1 # needed, as we are already adding one convModule to the list for the input
        dim_step = int((max_dim - min_dim) / n_convModules * 2) # factor of 2, as only every second layer has a dim increase
        if dim_step <= 0:
            dim_step = 1
        if n_convModules%2 == 0:
            skip_one_layer = False # Remember, we already created a layer and subtracted 1 from the number
        else:
            skip_one_layer = True
        in_out = [(x,x+dim_step) for x in range(min_dim, max_dim, dim_step)]
        if not in_out:
            in_out = [(min_dim, max_dim)]
        if in_out[-1][-1] != max_dim:
            in_out[-1] = (in_out[-1][0],max_dim)


        layer_list.append(nn.ReLU(inplace=True)) #second part of the input convolution, the ReLU

        counter = 1
        for (i,out) in in_out:
            if counter == 1 and skip_one_layer:
                pass
            else:
                layer_list.append(nn.Conv1d(i, i, kernel_size=self.conv_kernel_size))
                layer_list.append(nn.ReLU(inplace=True))
                channel_value = conv1d_out(channel_value)

                counter += 1
                if counter%batch_pool_freq == 0:
                    layer_list.append(nn.MaxPool1d(kernel_size=self.pool_kernel_size, stride=self.pool_kernel_size))
                    layer_list.append(nn.BatchNorm1d(i))
                    channel_value = maxpool1d_out(channel_value)

                if counter%dropout_freq == 0:
                    layer_list.append(nn.Dropout(p=dropout))
            layer_list.append(nn.Conv1d(i, out, kernel_size=self.conv_kernel_size))
            channel_value = conv1d_out(channel_value)

            counter += 1
            layer_list.append(nn.ReLU(inplace=True))
            if counter%batch_pool_freq == 0:
                layer_list.append(nn.MaxPool1d(kernel_size=self.pool_kernel_size, stride=self.pool_kernel_size))
                layer_list.append(nn.BatchNorm1d(out))
                channel_value = maxpool1d_out(channel_value)
            if counter%dropout_freq == 0:
                layer_list.append(nn.Dropout(p=dropout))
        layer_list.append(nn.BatchNorm1d(out))
        layer_list.append(nn.Dropout(p=dropout))
        self._n_channels = channel_value
        return nn.Sequential(*layer_list)
