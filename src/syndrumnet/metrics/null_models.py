"""
Null model generation and z-score normalization.

Implements degree-preserving network randomization for statistical
significance testing of proximity scores.
"""

import logging
import random
from typing import Dict, List, Set

import networkx as nx
import numpy as np

logger = logging.getLogger(__name__)


def degree_preserving_randomization(
    G: nx.Graph,
    module: Set[str],
    n_random: int = 1000,
    seed: Optional[int] = None,
) -> List[Set[str]]:
    """
    Generate degree-preserving random gene sets.
    
    For each gene in the module, sample a random gene with similar degree
    from the network. This preserves the degree distribution of the module
    under randomization.
    
    Parameters
    ----------
    G : nx.Graph
        Network graph.
    module : set
        Original gene module.
    n_random : int
        Number of random modules to generate.
    seed : int, optional
        Random seed for reproducibility.
        
    Returns
    -------
    list of set
        List of randomized gene sets.
        
    Notes
    -----
    This is the standard approach for network proximity null models
    (Guney et al., 2016; Iida et al., 2024).
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    # Filter module to genes in network
    module = module & set(G.nodes())
    module_size = len(module)
    
    if module_size == 0:
        logger.warning("Empty module after filtering to network")
        return []
    
    # Build degree bins
    degrees = dict(G.degree())
    degree_bins = _build_degree_bins(degrees, n_bins=20)
    
    # Get degree bin for each gene in module
    module_bins = [_get_bin(degrees[gene], degree_bins) for gene in module]
    
    # Generate random modules
    random_modules = []
    all_genes = list(G.nodes())
    
    for _ in range(n_random):
        random_module = set()
        
        for bin_idx in module_bins:
            # Sample from same degree bin
            candidates = degree_bins[bin_idx]
            
            if len(candidates) > 0:
                gene = random.choice(candidates)
                random_module.add(gene)
        
        random_modules.append(random_module)
    
    logger.debug(f"Generated {n_random} random modules of size {module_size}")
    
    return random_modules


def _build_degree_bins(
    degrees: Dict[str, int],
    n_bins: int = 20,
) -> List[List[str]]:
    """
    Build degree bins for stratified sampling.
    
    Parameters
    ----------
    degrees : dict
        {gene: degree}
    n_bins : int
        Number of bins.
        
    Returns
    -------
    list of list
        Bins containing genes with similar degrees.
    """
    # Get degree percentiles
    degree_values = sorted(degrees.values())
    percentiles = np.linspace(0, 100, n_bins + 1)
    bin_edges = np.percentile(degree_values, percentiles)
    
    # Assign genes to bins
    bins = [[] for _ in range(n_bins)]
    
    for gene, degree in degrees.items():
        bin_idx = _get_bin(degree, bins, bin_edges)
        bins[bin_idx].append(gene)
    
    return bins


def _get_bin(degree: int, bin_edges: np.ndarray) -> int:
    """Get bin index for a degree value."""
    for i in range(len(bin_edges) - 1):
        if bin_edges[i] <= degree < bin_edges[i + 1]:
            return i
    return len(bin_edges) - 2  # Last bin


def compute_zscore(
    observed: float,
    null_distribution: List[float],
) -> float:
    """
    Compute z-score of observed value against null distribution.
    
    Z = (observed - mean(null)) / std(null)
    
    Parameters
    ----------
    observed : float
        Observed proximity/distance value.
    null_distribution : list of float
        Null distribution from randomizations.
        
    Returns
    -------
    float
        Z-score. Negative values indicate observed < expected (closer proximity).
    """
    if len(null_distribution) == 0:
        logger.warning("Empty null distribution")
        return 0.0
    
    null_mean = np.mean(null_distribution)
    null_std = np.std(null_distribution)
    
    if null_std == 0:
        logger.warning("Zero standard deviation in null distribution")
        return 0.0
    
    z = (observed - null_mean) / null_std
    
    return z


def compute_normalized_proximity(
    G: nx.Graph,
    disease_module: Set[str],
    drug_module: Set[str],
    n_random: int = 1000,
    seed: Optional[int] = None,
) -> Tuple[float, float, float]:
    """
    Compute z-score normalized proximity between disease and drug modules.
    
    Parameters
    ----------
    G : nx.Graph
        Network graph.
    disease_module : set
        Disease gene module.
    drug_module : set
        Drug gene module.
    n_random : int
        Number of randomizations.
    seed : int, optional
        Random seed.
        
    Returns
    -------
    tuple
        (raw_proximity, z_score, p_value)
    """
    from syndrumnet.metrics.distances import shortest_path_distance
    
    # Observed proximity
    observed = shortest_path_distance(G, disease_module, drug_module)
    
    # Generate random modules
    random_drug_modules = degree_preserving_randomization(
        G, drug_module, n_random, seed
    )
    
    # Compute null distribution
    null_proximities = []
    for random_module in random_drug_modules:
        null_prox = shortest_path_distance(G, disease_module, random_module)
        null_proximities.append(null_prox)
    
    # Z-score
    z = compute_zscore(observed, null_proximities)
    
    # Empirical p-value (one-tailed: observed < null)
    p_value = np.mean([null_p <= observed for null_p in null_proximities])
    
    return observed, z, p_value