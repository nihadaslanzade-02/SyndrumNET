"""
Transcriptional correlation score (CQAB).

Measures correlation between disease and drug expression signatures.
"""

import logging
from typing import Dict, Set, Tuple

logger = logging.getLogger(__name__)


def compute_cqab(
    disease_signature: Dict[str, float],
    drug_a_signature_up: Set[str],
    drug_a_signature_down: Set[str],
    drug_b_signature_up: Set[str],
    drug_b_signature_down: Set[str],
) -> Tuple[float, float, float]:
    """
    Compute transcriptional correlation score CQAB.
    
    CQAB = (C_QA + C_QB) / 2
    
    where C_QA and C_QB measure inverse correlation between disease
    signature and drug signatures (drugs reverse disease).
    
    Parameters
    ----------
    disease_signature : dict
        {gene: fold_change} for disease.
    drug_a_signature_up : set
        Drug A up-regulated genes.
    drug_a_signature_down : set
        Drug A down-regulated genes.
    drug_b_signature_up : set
        Drug B up-regulated genes.
    drug_b_signature_down : set
        Drug B down-regulated genes.
        
    Returns
    -------
    tuple
        (cqab_score, cqa_score, cqb_score)
    """
    from syndrumnet.metrics.transcription import transcriptional_similarity
    
    logger.debug("Computing CQAB")
    
    # Compute transcriptional similarities (inverse correlation)
    c_qa = transcriptional_similarity(
        disease_signature,
        drug_a_signature_up,
        drug_a_signature_down,
        inverse_correlation=True,
    )
    
    c_qb = transcriptional_similarity(
        disease_signature,
        drug_b_signature_up,
        drug_b_signature_down,
        inverse_correlation=True,
    )
    
    # Average
    cqab = (c_qa + c_qb) / 2
    
    logger.debug(f"CQAB: C_QA={c_qa:.3f}, C_QB={c_qb:.3f}, CQAB={cqab:.3f}")
    
    return cqab, c_qa, c_qb


def compute_cqab_batch(
    disease_signature: Dict[str, float],
    drug_signatures: Dict[str, Dict[str, Set[str]]],
    drug_pairs: List[Tuple[str, str]],
) -> Dict[Tuple[str, str], Tuple[float, float, float]]:
    """
    Compute CQAB for multiple drug pairs.
    
    Parameters
    ----------
    disease_signature : dict
        Disease expression signature.
    drug_signatures : dict
        {drug_name: {'up': set, 'down': set}}
    drug_pairs : list of tuple
        Drug pair identifiers.
        
    Returns
    -------
    dict
        {(drug_a, drug_b): (cqab, cqa, cqb)}
    """
    results = {}
    
    for drug_a, drug_b in drug_pairs:
        if drug_a not in drug_signatures or drug_b not in drug_signatures:
            continue
        
        cqab, cqa, cqb = compute_cqab(
            disease_signature,
            drug_signatures[drug_a]['up'],
            drug_signatures[drug_a]['down'],
            drug_signatures[drug_b]['up'],
            drug_signatures[drug_b]['down'],
        )
        
        results[(drug_a, drug_b)] = (cqab, cqa, cqb)
    
    return results
