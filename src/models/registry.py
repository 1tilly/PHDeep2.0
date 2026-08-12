"""
Model registry — maps name strings to model classes.

Usage
-----
from src.models.registry import build_model

model = build_model("deepsea", sequence_length=1000, n_targets=919)
model = build_model("jellyfish", sequence_length=1000, n_targets=919,
                    n_convModules=12, layer_dimensions=(320, 960))
"""
from __future__ import annotations

import torch.nn as nn

from src.models.blainville.architecture import BlainvilleDeepSEA
from src.models.bottlenose.architecture import BottleNoseDeepSEA
from src.models.deeper_deepsea.architecture import DeeperDeepSEA
from src.models.deepsea.architecture import DeepSEA
from src.models.jellyfish.architecture import JellyFishDeepSEA

MODEL_REGISTRY: dict[str, type[nn.Module]] = {
    "deepsea": DeepSEA,
    "deeper_deepsea": DeeperDeepSEA,
    "blainville": BlainvilleDeepSEA,
    "bottlenose": BottleNoseDeepSEA,
    "jellyfish": JellyFishDeepSEA,
}


def build_model(name: str, sequence_length: int, n_targets: int, **kwargs) -> nn.Module:
    """
    Instantiate a model by name.

    Parameters
    ----------
    name : str
        One of: deepsea, deeper_deepsea, blainville, bottlenose, jellyfish.
    sequence_length : int
    n_targets : int
    **kwargs
        Extra keyword arguments forwarded to the model constructor
        (e.g. n_convModules and layer_dimensions for JellyFish).

    Raises
    ------
    ValueError
        If name is not in the registry.
    """
    if name not in MODEL_REGISTRY:
        available = sorted(MODEL_REGISTRY.keys())
        raise ValueError(
            f"Unknown model '{name}'. Available models: {available}"
        )
    return MODEL_REGISTRY[name](sequence_length, n_targets, **kwargs)
