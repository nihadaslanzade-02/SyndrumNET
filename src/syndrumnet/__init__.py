"""
SyndrumNET: Network-based prediction of synergistic drug combinations.

This package implements the complete SyndrumNET pipeline from Iida et al. (2024)
and extends it per Aslanzadə's thesis.
"""

__version__ = "1.0.0"
__author__ = "Nihad Aslanzade"

from syndrumnet.utils.config import load_config
from syndrumnet.utils.logging import setup_logger
from syndrumnet.utils.seeds import set_random_seed

__all__ = ["load_config", "setup_logger", "set_random_seed", "__version__"]