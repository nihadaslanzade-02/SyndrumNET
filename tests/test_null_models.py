"""Tests for null model randomization."""

import pytest
import networkx as nx

from syndrumnet.metrics.null_models import (
    degree_preserving_randomization,
    compute_zscore,
)


def test_degree_preserving_randomization():
    """Test degree-preserving randomization."""
    G = nx.karate_club_graph()
    
    module = {0, 1, 2, 3, 4}
    
    random_modules = degree_preserving_randomization(G, module, n_random=100, seed=42)
    
    assert len(random_modules) == 100
    assert all(len(rm) == len(module) for rm in random_modules)


def test_compute_zscore():
    """Test z-score computation."""
    observed = 5.0
    null = [10.0, 11.0, 9.0, 10.5, 10.2]
    
    z = compute_zscore(observed, null)
    
    # Observed is lower than null mean
    assert z < 0
