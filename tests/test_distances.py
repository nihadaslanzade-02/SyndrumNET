"""Tests for distance and proximity calculations."""

import pytest
import networkx as nx

from syndrumnet.metrics.distances import (
    shortest_path_distance,
    module_proximity,
    separation_score,
)


def test_shortest_path_distance():
    """Test shortest path distance calculation."""
    # Create simple graph
    G = nx.Graph()
    G.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
    
    source = {'A'}
    target = {'D'}
    
    dist = shortest_path_distance(G, source, target)
    
    assert dist == 3.0  # A->B->C->D


def test_module_proximity():
    """Test bidirectional module proximity."""
    G = nx.karate_club_graph()
    
    module_a = {0, 1, 2}
    module_b = {30, 31, 32}
    
    d_ab, d_ba = module_proximity(G, module_a, module_b)
    
    assert d_ab > 0
    assert d_ba > 0


def test_separation_score():
    """Test separation score calculation."""
    # Create graph with two separated communities
    G = nx.Graph()
    
    # Community 1
    G.add_edges_from([('A1', 'A2'), ('A2', 'A3'), ('A3', 'A1')])
    
    # Community 2
    G.add_edges_from([('B1', 'B2'), ('B2', 'B3'), ('B3', 'B1')])
    
    # Bridge
    G.add_edge('A1', 'B1')
    
    module_a = {'A1', 'A2', 'A3'}
    module_b = {'B1', 'B2', 'B3'}
    
    s_ab = separation_score(G, module_a, module_b)
    
    # Separated communities should have positive separation
    assert s_ab > 0
