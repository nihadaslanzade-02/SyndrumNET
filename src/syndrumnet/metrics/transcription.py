"""
Transcriptional correlation metrics.

Computes correlations between disease and drug expression signatures
for the CQAB component of the prediction score.
"""

import logging
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from scipy.stats import spearmanr, pearsonr

logger = logging.getLogger(__name__)


def compute_correlation(
    signature_a: Dict[str, float],
    signature_b: Dict[str, float],
    method: str = 'spearman',
) -> float:
    """
    Compute correlation between two gene expression signatures.
    
    Parameters
    ----------
    signature_a : dict
        {gene: fold_change} for signature A.
    signature_b : dict
        {gene: fold_change} for signature B.
    method : str
        Correlation method ('spearman' or 'pearson').
        
    Returns
    -------
    float
        Correlation coefficient.
        
    Notes
    -----
    Only genes present in both signatures are used for correlation.
    """
    # Find common genes
    common_genes = set(signature_a.keys()) & set(signature_b.keys())
    
    if len(common_genes) < 3:
        logger.warning(f"Only {len(common_genes)} common genes for correlation")
        return 0.0
    
    # Extract values for common genes
    values_a = np.array([signature_a[g] for g in common_genes])
    values_b = np.array([signature_b[g] for g in common_genes])
    
    # Compute correlation
    if method == 'spearman':
        corr, _ = spearmanr(values_a, values_b)
    elif method == 'pearson':
        corr, _ = pearsonr(values_a, values_b)
    else:
        raise ValueError(f"Unknown correlation method: {method}")
    
    # Handle NaN (e.g., from constant arrays)
    if np.isnan(corr):
        corr = 0.0
    
    return corr


def transcriptional_similarity(
    disease_signature: Dict[str, float],
    drug_signature_up: Set[str],
    drug_signature_down: Set[str],
    inverse_correlation: bool = True,
) -> float:
    """
    Compute transcriptional similarity between disease and drug.
    
    Measures whether drug reverses disease signature (inverse correlation).
    
    Parameters
    ----------
    disease_signature : dict
        Disease expression signature {gene: fold_change}.
    drug_signature_up : set
        Drug up-regulated genes.
    drug_signature_down : set
        Drug down-regulated genes.
    inverse_correlation : bool
        If True, compute -correlation (drug reversal of disease).
        
    Returns
    -------
    float
        Transcriptional similarity score (higher = more similar/reversed).
        
    Notes
    -----
    For synergy prediction, we expect drugs to have INVERSE correlation
    with disease (i.e., drugs reverse the disease signature).
    """
    # Convert drug signature to signed dictionary
    drug_sig = {}
    for gene in drug_signature_up:
        drug_sig[gene] = 1.0
    for gene in drug_signature_down:
        drug_sig[gene] = -1.0
    
    # Compute correlation
    corr = compute_correlation(disease_signature, drug_sig, method='spearman')
    
    # Invert if looking for reversal
    if inverse_correlation:
        corr = -corr
    
    return corr


def aggregate_transcriptional_scores(
    scores: List[float],
    method: str = 'mean',
) -> float:
    """
    Aggregate multiple transcriptional scores.
    
    Parameters
    ----------
    scores : list of float
        Individual scores to aggregate.
    method : str
        Aggregation method ('mean', 'median', 'max').
        
    Returns
    -------
    float
        Aggregated score.
    """
    if len(scores) == 0:
        return 0.0
    
    scores = np.array(scores)
    
    if method == 'mean':
        return np.mean(scores)
    elif method == 'median':
        return np.median(scores)
    elif method == 'max':
        return np.max(scores)
    else:
        raise ValueError(f"Unknown aggregation method: {method}")
