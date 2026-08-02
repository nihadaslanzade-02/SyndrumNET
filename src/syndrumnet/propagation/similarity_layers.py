"""
Similarity layers for multi-omics integration.

Computes disease-disease and drug-drug similarities for use as
propagation priors or regularization terms.

Every function here returns a bare ``np.ndarray`` with no row labels. Rows and
columns follow ``list(...)`` order over the input mapping, which since Python
3.7 is insertion order, so the caller can label the axes by keeping the same
key order it passed in.
"""

import logging
from typing import Callable, Dict, List, Set

import numpy as np

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
        Pairwise similarity matrix (diseases x diseases), rows and columns in
        ``list(disease_modules)`` order, 1.0 on the diagonal.

    Raises
    ------
    ValueError
        If `method` is not a known metric. The check runs before the loop, so
        an unknown name is rejected even for inputs of zero or one disease,
        where no off-diagonal comparison would otherwise be reached.
    """
    similarity_func = _resolve_set_metric(method)

    logger.info(f"Computing disease similarity using {method}")

    disease_names = list(disease_modules)
    n_diseases = len(disease_names)

    similarity = np.eye(n_diseases)

    for i in range(n_diseases):
        genes_i = disease_modules[disease_names[i]]

        for j in range(i + 1, n_diseases):
            sim = similarity_func(genes_i, disease_modules[disease_names[j]])
            similarity[i, j] = sim
            similarity[j, i] = sim

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
        {drug_name: fingerprint_array}. All fingerprints must have the same
        length; for 'tanimoto' every entry must additionally be 0 or 1.
    method : str
        Similarity metric ('tanimoto', 'cosine').

    Returns
    -------
    np.ndarray
        Pairwise similarity matrix (drugs x drugs), rows and columns in
        ``list(drug_fingerprints)`` order.

    Raises
    ------
    ValueError
        If `method` is unknown, if the mapping is empty, or if the
        fingerprints are of unequal length.
    """
    if method not in ('tanimoto', 'cosine'):
        raise ValueError(
            f"Unknown method: {method!r}. Expected 'tanimoto' or 'cosine'."
        )

    if not drug_fingerprints:
        raise ValueError("drug_fingerprints is empty, there is nothing to compare.")

    logger.info(f"Computing drug similarity using {method}")

    drug_names = list(drug_fingerprints)
    n_drugs = len(drug_names)

    lengths = {np.asarray(drug_fingerprints[d]).shape[-1] for d in drug_names}
    if len(lengths) > 1:
        raise ValueError(
            f"Fingerprints have mixed lengths {sorted(lengths)}. "
            "All molecules must be encoded over the same bit vector."
        )

    # Stack fingerprints
    fps = np.vstack([drug_fingerprints[d] for d in drug_names])

    if method == 'tanimoto':
        # Tanimoto similarity for binary fingerprints
        similarity = tanimoto_similarity_matrix(fps)
    else:
        # Cosine similarity
        from sklearn.metrics.pairwise import cosine_similarity
        similarity = cosine_similarity(fps)

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
    For now, uses Morgan fingerprints as proxy. See
    ``docs/PLACEHOLDERS.md`` for what that substitution costs.
    """
    try:
        from rdkit import Chem, DataStructs
        from rdkit.Chem import AllChem

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

    J(A,B) = |A intersect B| / |A union B|

    Two empty sets score 1.0: they are identical, and the alternative leaves
    the diagonal of a similarity matrix inconsistent for empty modules.
    An empty set against a non-empty one scores 0.0.
    """
    if len(set_a) == 0 and len(set_b) == 0:
        return 1.0

    return len(set_a & set_b) / len(set_a | set_b)


def overlap_coefficient(set_a: Set, set_b: Set) -> float:
    """
    Overlap coefficient.

    OC(A,B) = |A intersect B| / min(|A|, |B|)

    Follows the same empty-set convention as `jaccard_similarity`: two empty
    sets score 1.0, an empty set against a non-empty one scores 0.0.
    """
    if len(set_a) == 0 and len(set_b) == 0:
        return 1.0

    if len(set_a) == 0 or len(set_b) == 0:
        return 0.0

    return len(set_a & set_b) / min(len(set_a), len(set_b))


#: Set-overlap metrics accepted by `compute_disease_similarity`.
_SET_METRICS: Dict[str, Callable[[Set, Set], float]] = {
    'jaccard': jaccard_similarity,
    'overlap': overlap_coefficient,
}


def _resolve_set_metric(method: str) -> Callable[[Set, Set], float]:
    """Look up a set-overlap metric by name, or fail with the valid options."""
    try:
        return _SET_METRICS[method]
    except KeyError:
        raise ValueError(
            f"Unknown method: {method!r}. Expected one of {sorted(_SET_METRICS)}."
        ) from None


def tanimoto_similarity_matrix(fingerprints: np.ndarray) -> np.ndarray:
    """
    Compute pairwise Tanimoto similarity for binary fingerprints.

    Parameters
    ----------
    fingerprints : np.ndarray
        Binary fingerprint matrix (n_molecules x n_bits). Any dtype is
        accepted, including bool, but every entry must be 0 or 1.

    Returns
    -------
    np.ndarray
        Similarity matrix (n_molecules x n_molecules), 1.0 on the diagonal.

    Raises
    ------
    ValueError
        If the input is not two-dimensional, or holds a value other than 0/1.

    Notes
    -----
    The input is cast to an integer dtype before the matrix product, and that
    cast is load-bearing. A boolean array is the most natural spelling of
    "binary fingerprint matrix", but numpy's ``@`` on booleans is logical
    rather than arithmetic: every shared-bit count saturates at 1, so the
    numerator stops counting and the whole matrix, diagonal included, comes
    out wrong without raising anything.
    """
    fps = np.asarray(fingerprints)

    if fps.ndim != 2:
        raise ValueError(
            f"Expected a 2-D (n_molecules x n_bits) matrix, got shape {fps.shape}."
        )

    if not np.isin(fps, (0, 1)).all():
        raise ValueError(
            "Tanimoto similarity as defined here needs binary fingerprints, "
            "but the input holds values other than 0 and 1."
        )

    fps = fps.astype(np.int64)

    # Tanimoto = (A & B) / (A | B)
    # For binary: (A . B) / (|A| + |B| - A . B)
    dot_product = fps @ fps.T
    sizes = fps.sum(axis=1, keepdims=True)
    union = sizes + sizes.T - dot_product

    # Avoid division by zero. union == 0 only when both fingerprints are
    # all-zero, which the diagonal fix below then resolves to 1.0.
    union[union == 0] = 1

    similarity = dot_product / union

    # A molecule is identical to itself, including the all-zero case that the
    # ratio above would otherwise report as 0.
    np.fill_diagonal(similarity, 1.0)

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
        Function that computes similarity(entity_i, entity_j). Assumed
        symmetric: it is called once per unordered pair and the result is
        mirrored across the diagonal.
    **kwargs
        Additional arguments for similarity_func.

    Returns
    -------
    np.ndarray
        Similarity matrix, in `entities` order, 1.0 on the diagonal.
    """
    n = len(entities)
    matrix = np.eye(n)

    for i in range(n):
        for j in range(i + 1, n):
            sim = similarity_func(entities[i], entities[j], **kwargs)
            matrix[i, j] = sim
            matrix[j, i] = sim

    return matrix
