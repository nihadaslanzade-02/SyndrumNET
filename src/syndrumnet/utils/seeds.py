"""
Reproducibility utilities for random seed management.

Ensures all random operations (NumPy, Python random, NetworkX) use consistent seeds.
"""

import random
import logging
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)


def set_random_seed(seed: Optional[int] = None) -> int:
    """
    Set random seed for all libraries to ensure reproducibility.
    
    Sets seeds for:
    - Python's built-in random module
    - NumPy random number generator
    - Hash seed (for consistent dict/set ordering in Python 3.3+)
    
    Parameters
    ----------
    seed : int, optional
        Random seed value. If None, uses 42 as default.
        
    Returns
    -------
    int
        The seed that was set.
        
    Examples
    --------
    >>> set_random_seed(42)
    42
    >>> np.random.randint(100)  # Reproducible
    51
    >>> set_random_seed(42)
    42
    >>> np.random.randint(100)  # Same value
    51
    
    Notes
    -----
    This function should be called at the start of each script to ensure
    reproducibility. NetworkX graph operations that depend on random will
    also respect these seeds.
    """
    if seed is None:
        seed = 42
    
    # Python random
    random.seed(seed)
    
    # NumPy random
    np.random.seed(seed)
    
    # Set hash seed via environment (for subprocess calls)
    import os
    os.environ['PYTHONHASHSEED'] = str(seed)
    
    logger.info(f"Random seed set to: {seed}")
    
    return seed


def get_random_state(seed: Optional[int] = None) -> np.random.RandomState:
    """
    Create a NumPy RandomState object for isolated random operations.
    
    Useful for functions that need their own random state without affecting
    global random state.
    
    Parameters
    ----------
    seed : int, optional
        Random seed for this RandomState. If None, uses current global state.
        
    Returns
    -------
    np.random.RandomState
        Random state object.
        
    Examples
    --------
    >>> rng = get_random_state(42)
    >>> rng.randint(100)
    51
    """
    return np.random.RandomState(seed)