"""Tests for PRINCE propagation."""

import pytest
import networkx as nx
import numpy as np

from syndrumnet.propagation.prince import PRINCE


def test_prince_convergence():
    """Test PRINCE converges."""
    G = nx.karate_club_graph()
    
    prince = PRINCE(G, alpha=0.5, tolerance=1e-6, max_iterations=100)
    
    seed_nodes = {0, 1, 2}
    scores = prince.propagate(seed_nodes)
    
    assert len(scores) == G.number_of_nodes()
    assert all(s >= 0 for s in scores.values())
    assert sum(scores.values()) > 0


def test_prince_seed_weights():
    """Test PRINCE with custom seed weights."""
    G = nx.path_graph(5)
    
    prince = PRINCE(G, alpha=0.5)
    
    seed_nodes = {0, 4}
    seed_weights = {0: 1.0, 4: 0.5}
    
    scores = prince.propagate(seed_nodes, seed_weights)
    
    # Node 0 should have higher score than node 4
    assert scores[0] > scores[4]