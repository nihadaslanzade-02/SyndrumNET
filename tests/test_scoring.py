"""Tests for scoring functions."""

import pytest
import networkx as nx

from syndrumnet.scoring.tqab import compute_tqab
from syndrumnet.scoring.pqab import compute_pqab
from syndrumnet.scoring.cqab import compute_cqab


def test_tqab_complementary():
    """Test TQAB identifies complementary drugs."""
    # Create graph
    G = nx.Graph()
    
    # Disease module in center
    disease = {'D1', 'D2', 'D3'}
    
    # Drug A on left
    drug_a = {'A1', 'A2'}
    
    # Drug B on right
    drug_b = {'B1', 'B2'}
    
    # Connect
    G.add_edges_from([
        ('A1', 'A2'),
        ('A2', 'D1'),
        ('D1', 'D2'), ('D2', 'D3'),
        ('D3', 'B1'),
        ('B1', 'B2'),
    ])
    
    tqab, topo_class = compute_tqab(G, disease, drug_a, drug_b)
    
    # Should be complementary or intermediate
    assert topo_class in ['complementary', 'intermediate']


def test_cqab_inverse_correlation():
    """Test CQAB computes inverse correlation."""
    disease_sig = {'G1': 2.0, 'G2': 1.5, 'G3': -1.0}
    
    # Drug A reverses disease (up genes are down in disease)
    drug_a_up = {'G3'}  # Was down in disease
    drug_a_down = {'G1', 'G2'}  # Were up in disease
    
    # Drug B same as disease (redundant)
    drug_b_up = {'G1', 'G2'}
    drug_b_down = {'G3'}
    
    cqab, cqa, cqb = compute_cqab(
        disease_sig,
        drug_a_up, drug_a_down,
        drug_b_up, drug_b_down,
    )
    
    # Drug A should have positive correlation (reverses disease)
    # Drug B should have negative correlation (same as disease)
    assert cqa > cqb