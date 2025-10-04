"""
Network distance and proximity calculations.

Implements shortest-path based proximity between gene sets as defined
in Iida et al. (2024) and Guney et al. (2016).
"""

import logging
from typing import Dict, Optional, Set, Tuple

import networkx as nx
import numpy as np

logger = logging.getLogger(__name__)


def shortest_path_distance(
    G: nx.Graph,
    source_set: Set[str],
    target_set: Set[str],
    infinity_value: float = 1000.0,
) -> float:
    """
    Compute average shortest path distance from source set to target set.
    
    Implements: d(S,T) = (1/|S|) * sum_{s in S} min_{t in T} dist(s,t)
    
    Parameters
    ----------
    G : nx.Graph
        Network graph.
    source_set : set
        Source gene set.
    target_set : set
        Target gene set.
    infinity_value : float
        Value to use for disconnected nodes.
        
    Returns
    -------
    float
        Average minimum distance.
        
    Notes
    -----
    This is the proximity measure used in Guney et al. (2016) and
    adopted by Iida et al. (2024) for disease-drug proximity.
    """
    # Filter to genes in network
    source_set = source_set & set(G.nodes())
    target_set = target_set & set(G.nodes())
    
    if len(source_set) == 0 or len(target_set) == 0:
        logger.warning("Empty source or target set after filtering to network")
        return infinity_value
    
    total_distance = 0.0
    
    for source in source_set:
        # Find minimum distance to any target
        min_dist = infinity_value
        
        for target in target_set:
            try:
                dist = nx.shortest_path_length(G, source, target)
                min_dist = min(min_dist, dist)
            except nx.NetworkXNoPath:
                # Nodes are disconnected
                pass
        
        total_distance += min_dist
    
    avg_distance = total_distance / len(source_set)
    
    return avg_distance


def module_proximity(
    G: nx.Graph,
    module_a: Set[str],
    module_b: Set[str],
) -> Tuple[float, float]:
    """
    Compute bidirectional proximity between two modules.
    
    Returns both d(A,B) and d(B,A) as they may differ.
    
    Parameters
    ----------
    G : nx.Graph
        Network graph.
    module_a : set
        First gene module.
    module_b : set
        Second gene module.
        
    Returns
    -------
    tuple of float
        (d(A,B), d(B,A))
    """
    d_ab = shortest_path_distance(G, module_a, module_b)
    d_ba = shortest_path_distance(G, module_b, module_a)
    
    return d_ab, d_ba


def separation_score(
    G: nx.Graph,
    module_a: Set[str],
    module_b: Set[str],
) -> float:
    """
    Compute network separation s_AB between two modules.
    
    Implements: s_AB = <d_AB> - (<d_AA> + <d_BB>)/2
    
    where <d_XY> is the average shortest distance from module X to module Y.
    
    Parameters
    ----------
    G : nx.Graph
        Network graph.
    module_a : set
        First gene module.
    module_b : set
        Second gene module.
        
    Returns
    -------
    float
        Separation score. Positive values indicate separated modules,
        negative values indicate overlapping/proximal modules.
        
    Notes
    -----
    This metric captures whether two modules are complementary (separated)
    or redundant (overlapping) in network topology, which is key for
    the TQAB topological class score.
    """
    # Inter-module distances
    d_ab = shortest_path_distance(G, module_a, module_b)
    d_ba = shortest_path_distance(G, module_b, module_a)
    d_between = (d_ab + d_ba) / 2
    
    # Intra-module distances
    d_aa = shortest_path_distance(G, module_a, module_a)
    d_bb = shortest_path_distance(G, module_b, module_b)
    d_within = (d_aa + d_bb) / 2
    
    s_ab = d_between - d_within
    
    return s_ab


def compute_all_pairwise_distances(
    G: nx.Graph,
    gene_set: Set[str],
) -> Dict[Tuple[str, str], float]:
    """
    Compute all pairwise shortest path distances within a gene set.
    
    Parameters
    ----------
    G : nx.Graph
        Network graph.
    gene_set : set
        Set of genes.
        
    Returns
    -------
    dict
        {(gene_i, gene_j): distance}
    """
    genes = list(gene_set & set(G.nodes()))
    distances = {}
    
    # Use NetworkX all pairs shortest path for efficiency
    lengths = dict(nx.all_pairs_shortest_path_length(G))
    
    for i, gene_i in enumerate(genes):
        for gene_j in genes[i:]:
            if gene_i in lengths and gene_j in lengths[gene_i]:
                dist = lengths[gene_i][gene_j]
                distances[(gene_i, gene_j)] = dist
                distances[(gene_j, gene_i)] = dist
    
    return distances
