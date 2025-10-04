"""
Proximity score (PQAB).

Aggregates disease-drug network proximities for drug pair.
"""

import logging
from typing import Dict, Set, Tuple

import networkx as nx

logger = logging.getLogger(__name__)


def compute_pqab(
    G: nx.Graph,
    disease_module: Set[str],
    drug_a_module: Set[str],
    drug_b_module: Set[str],
    n_randomizations: int = 1000,
    seed: int = 42,
) -> Tuple[float, float, float]:
    """
    Compute proximity score PQAB.
    
    PQAB = (P_QA + P_QB) / 2
    
    where P_QA and P_QB are z-score normalized proximities.
    
    Parameters
    ----------
    G : nx.Graph
        Network graph.
    disease_module : set
        Disease module Q.
    drug_a_module : set
        Drug A module.
    drug_b_module : set
        Drug B module.
    n_randomizations : int
        Number of random modules for z-score normalization.
    seed : int
        Random seed.
        
    Returns
    -------
    tuple
        (pqab_score, pqa_zscore, pqb_zscore)
    """
    from syndrumnet.metrics.null_models import compute_normalized_proximity
    
    logger.debug("Computing PQAB")
    
    # Compute z-score normalized proximities
    _, z_qa, _ = compute_normalized_proximity(
        G, disease_module, drug_a_module, n_randomizations, seed
    )
    
    _, z_qb, _ = compute_normalized_proximity(
        G, disease_module, drug_b_module, n_randomizations, seed + 1
    )
    
    # Average z-scores (negative z-score = closer than random)
    pqab = (z_qa + z_qb) / 2
    
    # Invert sign so that negative distances (closer) become positive scores
    pqab = -pqab
    
    logger.debug(f"PQAB: P_QA={z_qa:.3f}, P_QB={z_qb:.3f}, PQAB={pqab:.3f}")
    
    return pqab, z_qa, z_qb


def compute_pqab_batch(
    G: nx.Graph,
    disease_module: Set[str],
    drug_modules: Dict[str, Set[str]],
    drug_pairs: List[Tuple[str, str]],
    n_randomizations: int = 1000,
    seed: int = 42,
) -> Dict[Tuple[str, str], Tuple[float, float, float]]:
    """
    Compute PQAB for multiple drug pairs.
    
    Parameters
    ----------
    G : nx.Graph
        Network graph.
    disease_module : set
        Disease module.
    drug_modules : dict
        {drug_name: gene_module}
    drug_pairs : list of tuple
        Drug pair identifiers.
    n_randomizations : int
        Number of randomizations for z-scores.
    seed : int
        Random seed.
        
    Returns
    -------
    dict
        {(drug_a, drug_b): (pqab, pqa, pqb)}
    """
    results = {}
    
    for drug_a, drug_b in drug_pairs:
        if drug_a not in drug_modules or drug_b not in drug_modules:
            continue
        
        pqab, pqa, pqb = compute_pqab(
            G,
            disease_module,
            drug_modules[drug_a],
            drug_modules[drug_b],
            n_randomizations,
            seed,
        )
        
        results[(drug_a, drug_b)] = (pqab, pqa, pqb)
    
    return results
