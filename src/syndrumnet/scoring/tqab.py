"""
Topological class score (TQAB).

Computes complementarity/redundancy of drug pair localization relative
to disease module based on network topology.
"""

import logging
from typing import Dict, Set, Tuple

import networkx as nx
import numpy as np

from syndrumnet.metrics.distances import separation_score

logger = logging.getLogger(__name__)


class TopologyClass:
    """Enumeration of drug pair topology classes."""
    COMPLEMENTARY = "complementary"
    REDUNDANT = "redundant"
    INTERMEDIATE = "intermediate"


def compute_tqab(
    G: nx.Graph,
    disease_module: Set[str],
    drug_a_module: Set[str],
    drug_b_module: Set[str],
) -> Tuple[float, str]:
    """
    Compute topological class score TQAB.
    
    Classifies drug pair (A,B) relative to disease module Q as:
    - Complementary: Drugs target different parts of disease network
    - Redundant: Drugs target overlapping parts
    - Intermediate: Mixed topology
    
    The score captures whether combining drugs provides complementary
    coverage of disease pathways.
    
    Parameters
    ----------
    G : nx.Graph
        Network graph.
    disease_module : set
        Disease gene module Q.
    drug_a_module : set
        Drug A gene module.
    drug_b_module : set
        Drug B gene module.
        
    Returns
    -------
    tuple
        (tqab_score, topology_class)
        
    Notes
    -----
    Implementation follows Iida et al. (2024) topological classification:
    
    Complementary if:
    - s_AB > 0 (drugs separated)
    - Both drugs proximal to disease (d_AQ, d_BQ small)
    
    Redundant if:
    - s_AB < 0 (drugs overlapping)
    - Or one drug far from disease
    
    Score is computed as combination of:
    - Drug-drug separation s_AB
    - Disease-drug proximities d_AQ, d_BQ
    """
    from syndrumnet.metrics.distances import shortest_path_distance
    
    # Drug-drug separation
    s_ab = separation_score(G, drug_a_module, drug_b_module)
    
    # Disease-drug proximities
    d_aq = shortest_path_distance(G, drug_a_module, disease_module)
    d_bq = shortest_path_distance(G, drug_b_module, disease_module)
    
    # Mean disease-drug proximity
    d_q_mean = (d_aq + d_bq) / 2
    
    # Classify topology
    if s_ab > 0:
        # Drugs are separated
        if d_aq < 3.0 and d_bq < 3.0:
            # Both close to disease
            topology_class = TopologyClass.COMPLEMENTARY
            score = 1.0 - d_q_mean / 10.0  # Normalize proximity
        else:
            # At least one far from disease
            topology_class = TopologyClass.INTERMEDIATE
            score = 0.5 - d_q_mean / 10.0
    else:
        # Drugs are overlapping
        topology_class = TopologyClass.REDUNDANT
        score = -abs(s_ab) / 5.0  # Penalize redundancy
    
    logger.debug(
        f"TQAB: s_AB={s_ab:.3f}, d_AQ={d_aq:.3f}, d_BQ={d_bq:.3f}, "
        f"class={topology_class}, score={score:.3f}"
    )
    
    return score, topology_class


def compute_tqab_batch(
    G: nx.Graph,
    disease_module: Set[str],
    drug_modules: Dict[str, Set[str]],
    drug_pairs: List[Tuple[str, str]],
) -> Dict[Tuple[str, str], Tuple[float, str]]:
    """
    Compute TQAB for multiple drug pairs.
    
    Parameters
    ----------
    G : nx.Graph
        Network graph.
    disease_module : set
        Disease gene module.
    drug_modules : dict
        {drug_name: gene_module}
    drug_pairs : list of tuple
        Drug pair identifiers [(drug_a, drug_b), ...].
        
    Returns
    -------
    dict
        {(drug_a, drug_b): (tqab_score, topology_class)}
    """
    results = {}
    
    for drug_a, drug_b in drug_pairs:
        if drug_a not in drug_modules or drug_b not in drug_modules:
            logger.warning(f"Missing module for pair ({drug_a}, {drug_b})")
            continue
        
        tqab, topo_class = compute_tqab(
            G,
            disease_module,
            drug_modules[drug_a],
            drug_modules[drug_b],
        )
        
        results[(drug_a, drug_b)] = (tqab, topo_class)
    
    return results
