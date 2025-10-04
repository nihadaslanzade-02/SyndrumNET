"""
PRINCE (PRopagation and INference for Complex Entities) algorithm.

Implements network propagation with restart for disease/drug module
prioritization, following the formalism from Vanunu et al. (2010)
and applied in Iida et al. (2024).
"""

import logging
from typing import Dict, Optional, Set, Tuple

import networkx as nx
import numpy as np
from scipy import sparse

logger = logging.getLogger(__name__)


class PRINCE:
    """
    PRINCE network propagation algorithm.
    
    Performs random walk with restart on a network to propagate signal
    from seed nodes (disease/drug modules) across the network.
    
    The algorithm iteratively updates node scores according to:
        F^(t+1) = α * W * F^(t) + (1-α) * F^(0)
    
    where:
    - F^(t) is the score vector at iteration t
    - W is the normalized adjacency matrix
    - α is the restart probability (typically 0.5)
    - F^(0) is the initial seed vector
    
    Parameters
    ----------
    network : nx.Graph
        Molecular interaction network.
    alpha : float
        Restart probability (0 to 1). Higher values propagate further.
    tolerance : float
        Convergence threshold.
    max_iterations : int
        Maximum number of iterations.
    normalize : str
        Normalization method ('column', 'row', 'symmetric').
        
    Examples
    --------
    >>> G = load_network()
    >>> prince = PRINCE(G, alpha=0.5)
    >>> seed_genes = {'TP53', 'BRCA1', 'EGFR'}
    >>> scores = prince.propagate(seed_genes)
    >>> top_genes = sorted(scores.items(), key=lambda x: x[1], reverse=True)[:10]
    
    References
    ----------
    Vanunu et al. (2010). "Associating genes and protein complexes with 
    disease via network propagation". PLoS Comput Biol.
    
    Cowen et al. (2017). "Network propagation: a universal amplifier of
    genetic associations". Nat Rev Genet.
    """
    
    def __init__(
        self,
        network: nx.Graph,
        alpha: float = 0.5,
        tolerance: float = 1e-6,
        max_iterations: int = 1000,
        normalize: str = 'column',
    ) -> None:
        """Initialize PRINCE propagator."""
        self.network = network
        self.alpha = alpha
        self.tolerance = tolerance
        self.max_iterations = max_iterations
        self.normalize = normalize
        
        # Build normalized adjacency matrix
        self._build_propagation_matrix()
        
        logger.info(
            f"PRINCE initialized: α={alpha}, tol={tolerance}, "
            f"max_iter={max_iterations}"
        )
    
    def _build_propagation_matrix(self) -> None:
        """Build normalized adjacency matrix for propagation."""
        logger.debug("Building propagation matrix")
        
        # Convert to adjacency matrix
        self.nodes = list(self.network.nodes())
        self.node_index = {node: i for i, node in enumerate(self.nodes)}
        
        n = len(self.nodes)
        A = nx.adjacency_matrix(self.network, nodelist=self.nodes)
        
        # Normalize based on method
        if self.normalize == 'column':
            # Column normalization: W_ij = A_ij / sum_i(A_ij)
            col_sums = np.array(A.sum(axis=0)).flatten()
            col_sums[col_sums == 0] = 1  # Avoid division by zero
            D_inv = sparse.diags(1.0 / col_sums)
            W = A @ D_inv
            
        elif self.normalize == 'row':
            # Row normalization: W_ij = A_ij / sum_j(A_ij)
            row_sums = np.array(A.sum(axis=1)).flatten()
            row_sums[row_sums == 0] = 1
            D_inv = sparse.diags(1.0 / row_sums)
            W = D_inv @ A
            
        elif self.normalize == 'symmetric':
            # Symmetric normalization: W = D^(-1/2) * A * D^(-1/2)
            degrees = np.array(A.sum(axis=1)).flatten()
            degrees[degrees == 0] = 1
            D_inv_sqrt = sparse.diags(1.0 / np.sqrt(degrees))
            W = D_inv_sqrt @ A @ D_inv_sqrt
            
        else:
            raise ValueError(f"Unknown normalization: {self.normalize}")
        
        self.W = W.tocsr()  # Convert to CSR for efficient operations
        
        logger.debug(f"Propagation matrix built: {n} nodes")
    
    def propagate(
        self,
        seed_nodes: Set[str],
        seed_weights: Optional[Dict[str, float]] = None,
    ) -> Dict[str, float]:
        """
        Propagate from seed nodes across the network.
        
        Parameters
        ----------
        seed_nodes : set
            Initial seed genes/proteins.
        seed_weights : dict, optional
            Custom weights for seed nodes. If None, uniform weights.
            
        Returns
        -------
        dict
            {gene: propagated_score} for all network nodes.
        """
        # Initialize seed vector
        F0 = self._build_seed_vector(seed_nodes, seed_weights)
        
        # Iterative propagation
        F = F0.copy()
        converged = False
        
        for iteration in range(self.max_iterations):
            F_new = self.alpha * self.W @ F + (1 - self.alpha) * F0
            
            # Check convergence
            diff = np.linalg.norm(F_new - F)
            
            if diff < self.tolerance:
                converged = True
                logger.debug(f"Converged at iteration {iteration + 1}")
                break
            
            F = F_new
        
        if not converged:
            logger.warning(
                f"Did not converge after {self.max_iterations} iterations "
                f"(diff={diff:.2e})"
            )
        
        # Convert to dictionary
        scores = {self.nodes[i]: F[i] for i in range(len(self.nodes))}
        
        return scores
    
    def _build_seed_vector(
        self,
        seed_nodes: Set[str],
        seed_weights: Optional[Dict[str, float]] = None,
    ) -> np.ndarray:
        """
        Build initial seed vector F^(0).
        
        Parameters
        ----------
        seed_nodes : set
            Seed genes.
        seed_weights : dict, optional
            Custom weights for seeds.
            
        Returns
        -------
        np.ndarray
            Initial score vector.
        """
        n = len(self.nodes)
        F0 = np.zeros(n)
        
        # Filter seeds to network nodes
        seed_nodes = seed_nodes & set(self.nodes)
        
        if len(seed_nodes) == 0:
            logger.warning("No seed nodes in network")
            return F0
        
        # Set seed weights
        for node in seed_nodes:
            idx = self.node_index[node]
            if seed_weights and node in seed_weights:
                F0[idx] = seed_weights[node]
            else:
                F0[idx] = 1.0
        
        # Normalize to sum to 1
        if F0.sum() > 0:
            F0 = F0 / F0.sum()
        
        return F0
    
    def propagate_multiple(
        self,
        modules: Dict[str, Set[str]],
    ) -> Dict[str, Dict[str, float]]:
        """
        Propagate from multiple modules.
        
        Parameters
        ----------
        modules : dict
            {module_name: seed_genes}
            
        Returns
        -------
        dict
            {module_name: {gene: score}}
        """
        results = {}
        
        for module_name, seed_nodes in modules.items():
            logger.debug(f"Propagating {module_name}")
            scores = self.propagate(seed_nodes)
            results[module_name] = scores
        
        return results
    
    def get_top_genes(
        self,
        scores: Dict[str, float],
        k: int = 100,
    ) -> List[Tuple[str, float]]:
        """
        Get top-k genes by propagation score.
        
        Parameters
        ----------
        scores : dict
            {gene: score}
        k : int
            Number of top genes to return.
            
        Returns
        -------
        list of tuple
            [(gene, score), ...] sorted by score descending.
        """
        sorted_genes = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        return sorted_genes[:k]
