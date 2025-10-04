"""
Disease and drug module construction.

Defines modules as sets of genes that characterize diseases or drug effects:
- Disease modules: susceptibility genes + dysregulated genes
- Drug modules: top 5% up/down regulated genes from L1000
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set

import pandas as pd
import networkx as nx

from syndrumnet.io.parsers import parse_creeds, parse_lincs
from syndrumnet.io.id_mapping import IDMapper

logger = logging.getLogger(__name__)


class ModuleBuilder:
    """
    Construct disease and drug modules.
    
    Disease modules are defined as the union of:
    1. Disease susceptibility genes (OMIM, ClinVar, GWAS, etc.)
    2. Dysregulated genes from expression signatures (CREEDS)
    
    Drug modules are defined as:
    - Top 5% up-regulated genes from L1000
    - Top 5% down-regulated genes from L1000
    
    Parameters
    ----------
    network : nx.Graph
        Molecular interaction network (used to filter modules to network genes).
    id_mapper : IDMapper
        Gene ID harmonization service.
        
    Examples
    --------
    >>> G = load_network()
    >>> mapper = IDMapper()
    >>> builder = ModuleBuilder(G, mapper)
    >>> disease_mods = builder.build_disease_modules('data/raw/creeds.txt')
    >>> drug_mods = builder.build_drug_modules('data/raw/lincs_sigs.txt')
    """
    
    def __init__(
        self,
        network: nx.Graph,
        id_mapper: IDMapper,
    ) -> None:
        """Initialize module builder."""
        self.network = network
        self.id_mapper = id_mapper
        self.network_genes = set(network.nodes())
        
        logger.info(f"Module builder initialized with {len(self.network_genes)} network genes")
    
    def build_disease_modules(
        self,
        creeds_file: Path,
        susceptibility_files: Optional[Dict[str, Path]] = None,
    ) -> Dict[str, Set[str]]:
        """
        Build disease modules from expression and susceptibility data.
        
        Parameters
        ----------
        creeds_file : Path
            CREEDS disease signatures file.
        susceptibility_files : dict, optional
            Paths to disease gene sources (OMIM, ClinVar, etc.).
            
        Returns
        -------
        dict
            {disease_name: set of genes}
        """
        logger.info("Building disease modules")
        
        # Parse CREEDS signatures
        signatures = parse_creeds(creeds_file)
        
        disease_modules = {}
        
        for disease, sig in signatures.items():
            # Combine up and down regulated genes
            sig_genes = set(sig['up']) | set(sig['down'])
            
            # Harmonize IDs
            sig_genes = set(self.id_mapper.harmonize_gene_list(list(sig_genes)))
            
            # Filter to network genes
            module_genes = sig_genes & self.network_genes
            
            # Add susceptibility genes if available
            if susceptibility_files:
                susc_genes = self._load_susceptibility_genes(disease, susceptibility_files)
                susc_genes = set(self.id_mapper.harmonize_gene_list(list(susc_genes)))
                module_genes |= (susc_genes & self.network_genes)
            
            if len(module_genes) > 0:
                disease_modules[disease] = module_genes
                logger.info(f"{disease}: {len(module_genes)} genes")
            else:
                logger.warning(f"No genes in network for disease: {disease}")
        
        logger.info(f"Built {len(disease_modules)} disease modules")
        
        return disease_modules
    
    def build_drug_modules(
        self,
        lincs_sig_file: Path,
        lincs_meta_file: Path,
        top_pct: float = 0.05,
    ) -> Dict[str, Dict[str, Set[str]]]:
        """
        Build drug modules from L1000 expression data.
        
        Parameters
        ----------
        lincs_sig_file : Path
            LINCS L1000 signatures file.
        lincs_meta_file : Path
            LINCS metadata file.
        top_pct : float
            Percentage of genes for up/down sets (default: 0.05 = 5%).
            
        Returns
        -------
        dict
            {drug_name: {'up': set, 'down': set}}
        """
        logger.info("Building drug modules")
        
        # Parse LINCS
        signatures = parse_lincs(lincs_sig_file, lincs_meta_file, top_pct)
        
        drug_modules = {}
        
        for drug, sig in signatures.items():
            # Harmonize IDs
            up_genes = set(self.id_mapper.harmonize_gene_list(sig['up']))
            down_genes = set(self.id_mapper.harmonize_gene_list(sig['down']))
            
            # Filter to network genes
            up_in_net = up_genes & self.network_genes
            down_in_net = down_genes & self.network_genes
            
            if len(up_in_net) > 0 or len(down_in_net) > 0:
                drug_modules[drug] = {
                    'up': up_in_net,
                    'down': down_in_net,
                }
                logger.debug(
                    f"{drug}: {len(up_in_net)} up, {len(down_in_net)} down"
                )
            else:
                logger.warning(f"No genes in network for drug: {drug}")
        
        logger.info(f"Built {len(drug_modules)} drug modules")
        
        return drug_modules
    
    def _load_susceptibility_genes(
        self,
        disease: str,
        source_files: Dict[str, Path],
    ) -> Set[str]:
        """
        Load disease susceptibility genes from multiple sources.
        
        Parameters
        ----------
        disease : str
            Disease name.
        source_files : dict
            {source_name: filepath}
            
        Returns
        -------
        set
            Disease susceptibility genes.
        """
        # This is a placeholder. Actual implementation requires parsing
        # OMIM, ClinVar, GWAS, DisGeNET for specific diseases.
        
        logger.debug(f"Loading susceptibility genes for {disease}")
        
        genes = set()
        
        # TODO: Implement parsers for each source
        # For now, return empty set
        
        return genes
    
    def save_modules(
        self,
        modules: Dict[str, Set[str]],
        output_path: Path,
    ) -> None:
        """
        Save modules to file.
        
        Parameters
        ----------
        modules : dict
            {module_name: set of genes}
        output_path : Path
            Output CSV path.
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        rows = []
        for name, genes in modules.items():
            for gene in genes:
                rows.append({'module': name, 'gene': gene})
        
        df = pd.DataFrame(rows)
        df.to_csv(output_path, index=False)
        
        logger.info(f"Saved {len(modules)} modules to {output_path}")
    
    @staticmethod
    def load_modules(filepath: Path) -> Dict[str, Set[str]]:
        """
        Load modules from file.
        
        Parameters
        ----------
        filepath : Path
            Module CSV file.
            
        Returns
        -------
        dict
            {module_name: set of genes}
        """
        df = pd.read_csv(filepath)
        
        modules = {}
        for module_name in df['module'].unique():
            genes = set(df[df['module'] == module_name]['gene'])
            modules[module_name] = genes
        
        logger.info(f"Loaded {len(modules)} modules from {filepath}")
        
        return modules