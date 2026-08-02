"""
Proximity score (PQAB).

Aggregates disease-drug network proximities for drug pair.
"""

import hashlib
import logging
from typing import Dict, List, Set, Tuple

import networkx as nx

logger = logging.getLogger(__name__)


def module_seed(base_seed: int, module: Set[str]) -> int:
    """
    Derive a deterministic null-model seed for one gene module.

    Parameters
    ----------
    base_seed : int
        Run-level seed from the configuration.
    module : set
        Gene module the seed is for.

    Returns
    -------
    int
        Seed depending only on the base seed and the module's contents.

    Notes
    -----
    P_QA is a property of the disease and one drug. It must not change with
    which other drug the pair happens to contain, or with the order the two
    drugs are passed in. Seeding from the module's own contents guarantees
    that: the same module always draws the same null distribution, and two
    different modules draw independent ones.

    SHA-256 is used rather than the builtin hash() because string hashing is
    salted per process, which would make results irreproducible across runs.
    """
    key = f"{base_seed}:" + "|".join(sorted(module))
    digest = hashlib.sha256(key.encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "big")


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

    Notes
    -----
    Each drug is seeded from its own module contents via module_seed(), so a
    drug's z-score is the same whichever side of the pair it is passed as.
    """
    z_qa = proximity_zscore(G, disease_module, drug_a_module, n_randomizations, seed)
    z_qb = proximity_zscore(G, disease_module, drug_b_module, n_randomizations, seed)

    # Average z-scores, then invert the sign: a negative z means the drug is
    # closer to the disease than chance, which should raise the final score.
    pqab = -(z_qa + z_qb) / 2

    logger.debug(f"PQAB: P_QA={z_qa:.3f}, P_QB={z_qb:.3f}, PQAB={pqab:.3f}")

    return pqab, z_qa, z_qb


def proximity_zscore(
    G: nx.Graph,
    disease_module: Set[str],
    drug_module: Set[str],
    n_randomizations: int = 1000,
    seed: int = 42,
) -> float:
    """
    Compute the z-scored disease-drug proximity P_QA for a single drug.

    Parameters
    ----------
    G : nx.Graph
        Network graph.
    disease_module : set
        Disease module Q.
    drug_module : set
        Drug module.
    n_randomizations : int
        Number of random modules for z-score normalization.
    seed : int
        Run-level seed; the actual null-model seed is derived per module.

    Returns
    -------
    float
        Z-score. Negative means closer to the disease than chance.
    """
    from syndrumnet.metrics.null_models import compute_normalized_proximity

    _, z, _ = compute_normalized_proximity(
        G,
        disease_module,
        drug_module,
        n_randomizations,
        module_seed(seed, drug_module),
    )
    return z


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

    Notes
    -----
    P_QA depends only on the disease and one drug, so each drug's z-score is
    computed once and reused across every pair it appears in. Computing it per
    pair instead costs O(n_drugs^2) null-model runs where O(n_drugs) suffice;
    at the 1,488 drugs the default config declares that is 2.2 million
    randomization sweeps rather than 1,488.
    """
    needed = {drug for pair in drug_pairs for drug in pair if drug in drug_modules}

    logger.info(
        f"Computing proximity z-scores for {len(needed)} drugs "
        f"covering {len(drug_pairs)} pairs"
    )

    zscores: Dict[str, float] = {
        drug: proximity_zscore(
            G, disease_module, drug_modules[drug], n_randomizations, seed
        )
        for drug in sorted(needed)
    }

    results = {}

    for drug_a, drug_b in drug_pairs:
        if drug_a not in zscores or drug_b not in zscores:
            continue

        z_qa, z_qb = zscores[drug_a], zscores[drug_b]
        results[(drug_a, drug_b)] = (-(z_qa + z_qb) / 2, z_qa, z_qb)

    return results
