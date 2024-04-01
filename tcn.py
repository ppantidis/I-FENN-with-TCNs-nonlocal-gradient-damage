import torch
import torch.nn as nn
from torch.nn.utils import weight_norm
import gc

gc.collect()

class Chomp1d(nn.Module):
    def __init__(self, chomp_size):
        super(Chomp1d, self).__init__()
        self.chomp_size = chomp_size

    def forward(self, x):
        return x[:, :, :-self.chomp_size].contiguous()


class TemporalBlock(nn.Module):
    def __init__(self, n_inputs, n_outputs, act_func, kernel_size, stride, dilation, padding, dropout=0.2):
        super(TemporalBlock, self).__init__()

        if act_func == 'tanh':
            self.activ1 = nn.Tanh()
            self.activ2 = nn.Tanh()
            self.activ  = nn.Tanh()
        elif act_func == 'relu':
            self.activ1 = nn.ReLU()
            self.activ2 = nn.ReLU()
            self.activ  = nn.ReLU()
        elif act_func == 'silu':
            self.activ1 = nn.SiLU()
            self.activ2 = nn.SiLU()
            self.activ  = nn.SiLU()

        self.conv1 = weight_norm(nn.Conv1d(n_inputs, n_outputs, kernel_size,
                                           stride=stride, padding=padding, dilation=dilation))
        self.chomp1 = Chomp1d(padding)

        self.dropout1 = nn.Dropout(dropout)

        self.conv2 = weight_norm(nn.Conv1d(n_outputs, n_outputs, kernel_size,
                                           stride=stride, padding=padding, dilation=dilation))
        self.chomp2 = Chomp1d(padding)
        self.dropout2 = nn.Dropout(dropout)

        self.net = nn.Sequential(self.conv1, self.chomp1, self.activ1, self.dropout1,
                                 self.conv2, self.chomp2, self.activ2, self.dropout2)
        self.downsample = nn.Conv1d(n_inputs, n_outputs, 1) if n_inputs != n_outputs else None
        self.init_weights()

    def init_weights(self):
        # self.conv1.weight.data.normal_(0, 0.1)
        # self.conv2.weight.data.normal_(0, 0.1)
        torch.nn.init.xavier_normal_(self.conv1.weight)
        torch.nn.init.xavier_normal_(self.conv2.weight)
        if self.downsample is not None:
            self.downsample.weight.data.normal_(0, 0.1) # Change that also

    def forward(self, x):
        out = self.net(x)
        res = x if self.downsample is None else self.downsample(x)
        return self.activ(out + res)


class TemporalConvNet(nn.Module):
    def __init__(self, num_inputs, num_channels, act_func, kernel_size=2, dropout=0.2):
        super(TemporalConvNet, self).__init__()
        layers = []
        num_levels = len(num_channels)
        for i in range(num_levels):
            dilation_size = 2 ** i
            in_channels = num_inputs if i == 0 else num_channels[i-1]
            out_channels = num_channels[i]
            layers += [TemporalBlock(in_channels, out_channels, act_func, kernel_size, stride=1, dilation=dilation_size,
                                     padding=(kernel_size-1) * dilation_size, dropout=dropout)]

        self.network = nn.Sequential(*layers)

    def forward(self, x):
        return self.network(x)
    