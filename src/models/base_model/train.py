import torch
from torch import nn


def criterion(loss=nn.BCELoss):
    """
    Specify the appropriate loss function (criterion) for this
    model.

    Returns
    -------
    torch.nn._Loss
    """
    return loss()

def get_optimizer(optimiser=torch.optim.SGD, lr=0.01, weight_decay=1e-6, momentum=0.9):
    """
    Specify an optimizer and its parameters.

    Returns
    -------
    tuple(torch.optim.Optimizer, dict)
        The optimizer class and the dictionary of kwargs that should
        be passed in to the optimizer constructor.

    """
    return (optimiser,
            {"lr": lr, "weight_decay": weight_decay, "momentum": momentum})
