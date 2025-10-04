"""
Similarity layers for multi-omics integration.

Computes disease-disease and drug-drug similarities for use as
propagation priors or regularization terms.
"""

import logging
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr

logger = logging.getLogger(__name__)


def compute_disease_similarity(
    disease_modules: Dict[str, Set[str]],
    method: str = 'jaccard',
) -> np.ndarray:
    """
    Compute pairwise disease similarity matrix.
    
    Based on overlap of disease modules (genes).
    
    Parameters
    ----------
    disease_modules : dict
        {disease_name: gene_set}
    method : str
        Similarity metric ('jaccard', 'overlap').
        
    Returns
    -------
    np.ndarray
        Pairwise similarity matrix (diseases x diseases).
    """
    logger.info(f"Computing disease similarity using {method}")
    
    disease_names = list(disease_modules.keys())
    n_diseases = len(disease_names)
    
    similarity = np.zeros((n_diseases, n_diseases))
    
    for i, disease_i in enumerate(disease_names):
        for j, disease_j in enumerate(disease_names):
            if i == j:
                similarity[i, j] = 1.0
            else:
                genes_i = disease_modules[disease_i]
                genes_j = disease_modules[disease_j]
                
                if method == 'jaccard':
                    sim = jaccard_similarity(genes_i, genes_j)
                elif method == 'overlap':
                    sim = overlap_coefficient(genes_i, genes_j)
                else:
                    raise ValueError(f"Unknown method: {method}")
                
                similarity[i, j] = sim
    
    logger.info(f"Disease similarity matrix: {n_diseases}x{n_diseases}")
    
    return similarity


def compute_drug_similarity(
    drug_fingerprints: Dict[str, np.ndarray],
    method: str = 'tanimoto',
) -> np.ndarray:
    """
    Compute pairwise drug similarity matrix.
    
    Based on chemical structure similarity (fingerprints).
    
    Parameters
    ----------
    drug_fingerprints : dict
        {drug_name: fingerprint_array}
    method : str
        Similarity metric ('tanimoto', 'cosine').
        
    Returns
    -------
    np.ndarray
        Pairwise similarity matrix (drugs x drugs).
    """
    logger.info(f"Computing drug similarity using {method}")
    
    drug_names = list(drug_fingerprints.keys())
    n_drugs = len(drug_names)
    
    # Stack fingerprints
    fps = np.vstack([drug_fingerprints[d] for d in drug_names])
    
    if method == 'tanimoto':
        # Tanimoto similarity for binary fingerprints
        similarity = tanimoto_similarity_matrix(fps)
    elif method == 'cosine':
        # Cosine similarity
        from sklearn.metrics.pairwise import cosine_similarity
        similarity = cosine_similarity(fps)
    else:
        raise ValueError(f"Unknown method: {method}")
    
    logger.info(f"Drug similarity matrix: {n_drugs}x{n_drugs}")
    
    return similarity


def kcf_fingerprint_similarity(
    smiles_a: str,
    smiles_b: str,
) -> float:
    """
    Compute KCF-S fingerprint similarity between two molecules.
    
    KCF-S (KEGG Chemical Function and Substructures) fingerprints
    capture chemical structure for drug similarity.
    
    Parameters
    ----------
    smiles_a : str
        SMILES string for molecule A.
    smiles_b : str
        SMILES string for molecule B.
        
    Returns
    -------
    float
        Tanimoto similarity (0 to 1).
        
    Notes
    -----
    This is a placeholder. Full implementation requires:
    - RDKit for SMILES parsing
    - KCF fingerprint generation
    For now, uses Morgan fingerprints as proxy.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import DataStructs
        
        # Parse SMILES
        mol_a = Chem.MolFromSmiles(smiles_a)
        mol_b = Chem.MolFromSmiles(smiles_b)
        
        if mol_a is None or mol_b is None:
            logger.warning("Invalid SMILES string")
            return 0.0
        
        # Generate Morgan fingerprints (as KCF-S proxy)
        fp_a = AllChem.GetMorganFingerprintAsBitVect(mol_a, radius=2, nBits=2048)
        fp_b = AllChem.GetMorganFingerprintAsBitVect(mol_b, radius=2, nBits=2048)
        
        # Tanimoto similarity
        similarity = DataStructs.TanimotoSimilarity(fp_a, fp_b)
        
        return similarity
        
    except ImportError:
        logger.warning("RDKit not available, returning 0 similarity")
        return 0.0


def jaccard_similarity(set_a: Set, set_b: Set) -> float:
    """
    Jaccard similarity coefficient.
    
    J(A,B) = |A ∩ B| / |A ∪ B|
    """
    if len(set_a) == 0 and len(set_b) == 0:
        return 1.0
    
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    
    if union == 0:
        return 0.0
    
    return intersection / union


def overlap_coefficient(set_a: Set, set_b: Set) -> float:
    """
    Overlap coefficient.
    
    OC(A,B) = |A ∩ B| / min(|A|, |B|)
    """
    if len(set_a) == 0 or len(set_b) == 0:
        return 0.0
    
    intersection = len(set_a & set_b)
    min_size = min(len(set_a), len(set_b))
    
    return intersection / min_size


def tanimoto_similarity_matrix(fingerprints: np.ndarray) -> np.ndarray:
    """
    Compute pairwise Tanimoto similarity for binary fingerprints.
    
    Parameters
    ----------
    fingerprints : np.ndarray
        Binary fingerprint matrix (n_molecules x n_bits).
        
    Returns
    -------
    np.ndarray
        Similarity matrix (n_molecules x n_molecules).
    """
    # Tanimoto = (A & B) / (A | B)
    # For binary: (A · B) / (|A| + |B| - A · B)
    
    dot_product = fingerprints @ fingerprints.T
    sizes = fingerprints.sum(axis=1, keepdims=True)
    union = sizes + sizes.T - dot_product
    
    # Avoid division by zero
    union[union == 0] = 1
    
    similarity = dot_product / union
    
    return similarity


def build_similarity_matrix(
    entities: List[str],
    similarity_func,
    **kwargs,
) -> np.ndarray:
    """
    Build pairwise similarity matrix for a list of entities.
    
    Parameters
    ----------
    entities : list
        Entity identifiers.
    similarity_func : callable
        Function that computes similarity(entity_i, entity_j).
    **kwargs
        Additional arguments for similarity_func.
        
    Returns
    -------
    np.ndarray
        Similarity matrix.
    """
    n = len(entities)
    matrix = np.zeros((n, n))
    
    for i, entity_i in enumerate(entities):
        for j, entity_j in enumerate(entities):
            if i == j:
                matrix[i, j] = 1.0
            elif i < j:
                sim = similarity_func(entity_i, entity_j, **kwargs)
                matrix[i, j] = sim
                matrix[j, i] = sim
    
    return matrix