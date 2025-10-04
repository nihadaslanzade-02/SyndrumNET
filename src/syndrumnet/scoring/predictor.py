"""
SyndrumNET synergy predictor.

Orchestrates computation of all scoring components and combines them
into final prediction scores.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import networkx as nx
import pandas as pd

from syndrumnet.scoring.tqab import compute_tqab_batch
from syndrumnet.scoring.pqab import compute_pqab_batch
from syndrumnet.scoring.cqab import compute_cqab_batch

logger = logging.getLogger(__name__)


class SynergyPredictor:
    """
    Complete SyndrumNET synergy prediction pipeline.
    
    Combines topological (TQAB), proximity (PQAB), and transcriptional
    (CQAB) scores to predict synergistic drug combinations.
    
    Final prediction: Score_Q,AB = TQAB + PQAB + CQAB
    
    Parameters
    ----------
    network : nx.Graph
        Molecular interaction network.
    n_randomizations : int
        Number of randomizations for proximity z-scores.
    seed : int
        Random seed for reproducibility.
        
    Examples
    --------
    >>> predictor = SynergyPredictor(G)
    >>> predictor.set_disease_modules(disease_modules)
    >>> predictor.set_drug_modules(drug_modules)
    >>> predictions = predictor.predict_all('asthma')
    >>> top_pairs = predictions.nlargest(10, 'prediction_score')
    """
    
    def __init__(
        self,
        network: nx.Graph,
        n_randomizations: int = 1000,
        seed: int = 42,
    ) -> None:
        """Initialize predictor."""
        self.network = network
        self.n_randomizations = n_randomizations
        self.seed = seed
        
        self.disease_modules: Optional[Dict[str, Set[str]]] = None
        self.drug_modules: Optional[Dict[str, Dict[str, Set[str]]]] = None
        self.disease_signatures: Optional[Dict[str, Dict[str, float]]] = None
        
        logger.info("SynergyPredictor initialized")
    
    def set_disease_modules(self, modules: Dict[str, Set[str]]) -> None:
        """Set disease modules."""
        self.disease_modules = modules
        logger.info(f"Loaded {len(modules)} disease modules")
    
    def set_drug_modules(self, modules: Dict[str, Dict[str, Set[str]]]) -> None:
        """Set drug modules."""
        self.drug_modules = modules
        logger.info(f"Loaded {len(modules)} drug modules")
    
    def set_disease_signatures(self, signatures: Dict[str, Dict[str, float]]) -> None:
        """Set disease expression signatures."""
        self.disease_signatures = signatures
        logger.info(f"Loaded {len(signatures)} disease signatures")
    
    def predict_all(
        self,
        disease: str,
        max_pairs: Optional[int] = None,
    ) -> pd.DataFrame:
        """
        Predict synergy scores for all drug pairs for a disease.
        
        Parameters
        ----------
        disease : str
            Disease name.
        max_pairs : int, optional
            Maximum number of drug pairs to evaluate (for testing).
            
        Returns
        -------
        pd.DataFrame
            Predictions with columns:
            - drug_a, drug_b
            - tqab, pqab, cqab
            - prediction_score (sum of components)
            - topology_class
        """
        if self.disease_modules is None or self.drug_modules is None:
            raise ValueError("Modules not set. Call set_disease_modules() and set_drug_modules().")
        
        if disease not in self.disease_modules:
            raise ValueError(f"Unknown disease: {disease}")
        
        logger.info(f"Predicting synergy for {disease}")
        
        disease_module = self.disease_modules[disease]
        drug_names = list(self.drug_modules.keys())
        
        # Generate all drug pairs
        drug_pairs = []
        for i, drug_a in enumerate(drug_names):
            for drug_b in drug_names[i+1:]:
                drug_pairs.append((drug_a, drug_b))
        
        if max_pairs:
            drug_pairs = drug_pairs[:max_pairs]
        
        logger.info(f"Evaluating {len(drug_pairs)} drug pairs")
        
        # Extract drug modules (union of up/down for topological analysis)
        drug_module_sets = {
            drug: self.drug_modules[drug]['up'] | self.drug_modules[drug]['down']
            for drug in drug_names
        }
        
        # Compute TQAB
        logger.info("Computing TQAB (topological)")
        tqab_results = compute_tqab_batch(
            self.network,
            disease_module,
            drug_module_sets,
            drug_pairs,
        )
        
        # Compute PQAB
        logger.info("Computing PQAB (proximity)")
        pqab_results = compute_pqab_batch(
            self.network,
            disease_module,
            drug_module_sets,
            drug_pairs,
            self.n_randomizations,
            self.seed,
        )
        
        # Compute CQAB (if signatures available)
        cqab_results = {}
        if self.disease_signatures and disease in self.disease_signatures:
            logger.info("Computing CQAB (transcriptional)")
            cqab_results = compute_cqab_batch(
                self.disease_signatures[disease],
                self.drug_modules,
                drug_pairs,
            )
        
        # Combine results
        predictions = []
        
        for drug_a, drug_b in drug_pairs:
            tqab, topo_class = tqab_results.get((drug_a, drug_b), (0.0, 'unknown'))
            pqab, pqa, pqb = pqab_results.get((drug_a, drug_b), (0.0, 0.0, 0.0))
            cqab, cqa, cqb = cqab_results.get((drug_a, drug_b), (0.0, 0.0, 0.0))
            
            # Final prediction
            prediction = tqab + pqab + cqab
            
            predictions.append({
                'disease': disease,
                'drug_a': drug_a,
                'drug_b': drug_b,
                'tqab': tqab,
                'pqab': pqab,
                'cqab': cqab,
                'prediction_score': prediction,
                'topology_class': topo_class,
                'pqa': pqa,
                'pqb': pqb,
                'cqa': cqa,
                'cqb': cqb,
            })
        
        df = pd.DataFrame(predictions)
        
        # Sort by prediction score
        df = df.sort_values('prediction_score', ascending=False)
        
        logger.info(f"Generated {len(df)} predictions for {disease}")
        
        return df
    
    def predict_multiple_diseases(
        self,
        diseases: List[str],
    ) -> Dict[str, pd.DataFrame]:
        """
        Predict synergy for multiple diseases.
        
        Parameters
        ----------
        diseases : list
            Disease names.
            
        Returns
        -------
        dict
            {disease: predictions_df}
        """
        results = {}
        
        for disease in diseases:
            logger.info(f"Processing disease: {disease}")
            df = self.predict_all(disease)
            results[disease] = df
        
        return results
    
    def save_predictions(
        self,
        predictions: pd.DataFrame,
        output_path: Path,
    ) -> None:
        """
        Save predictions to CSV.
        
        Parameters
        ----------
        predictions : pd.DataFrame
            Prediction results.
        output_path : Path
            Output file path.
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)
        predictions.to_csv(output_path, index=False)
        logger.info(f"Saved predictions to {output_path}")