"""
Gene/protein ID mapping and harmonization.

Converts between different identifier systems (Entrez, Ensembl, HGNC, UniProt)
to ensure consistent gene names across all data sources.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set

import pandas as pd
import mygene

logger = logging.getLogger(__name__)


class IDMapper:
    """
    Gene/protein ID mapping service.
    
    Harmonizes identifiers across data sources using:
    - HGNC as primary authority for gene symbols
    - mygene.info for batch conversions
    - Local caches to minimize API calls
    
    Parameters
    ----------
    hgnc_file : Path, optional
        Path to HGNC complete set file.
    cache_dir : Path, optional
        Directory for ID mapping cache.
        
    Examples
    --------
    >>> mapper = IDMapper()
    >>> mapper.to_hgnc(['ENSG00000141510', 'TP53'])
    ['TP53', 'TP53']
    >>> mapper.to_entrez(['TP53', 'BRCA1'])
    [7157, 672]
    """
    
    def __init__(
        self,
        hgnc_file: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
    ) -> None:
        """Initialize ID mapper."""
        self.cache_dir = Path(cache_dir) if cache_dir else Path("data/interim/id_cache")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Load HGNC data
        self.hgnc = None
        if hgnc_file and hgnc_file.exists():
            self._load_hgnc(hgnc_file)
        
        # Initialize mygene client
        self.mg = mygene.MyGeneInfo()
        
        # Cache for conversions
        self._cache: Dict[str, Dict] = {}
    
    def _load_hgnc(self, filepath: Path) -> None:
        """Load HGNC gene nomenclature data."""
        logger.info(f"Loading HGNC data from {filepath}")
        
        self.hgnc = pd.read_csv(filepath, sep='\t')
        
        # Build lookup dictionaries
        self.symbol_to_entrez = dict(zip(
            self.hgnc['symbol'],
            self.hgnc['entrez_id']
        ))
        
        self.entrez_to_symbol = {v: k for k, v in self.symbol_to_entrez.items()}
        
        logger.info(f"Loaded {len(self.hgnc)} HGNC genes")
    
    def to_hgnc(
        self,
        ids: List[str],
        from_type: str = 'auto',
    ) -> List[Optional[str]]:
        """
        Convert gene IDs to HGNC symbols.
        
        Parameters
        ----------
        ids : list of str
            Gene identifiers to convert.
        from_type : str
            Source ID type ('ensembl', 'entrez', 'uniprot', 'auto').
            If 'auto', attempts to detect type.
            
        Returns
        -------
        list of str or None
            HGNC symbols. None for unmapped IDs.
        """
        results = []
        
        for gene_id in ids:
            # Check cache
            if gene_id in self._cache:
                results.append(self._cache[gene_id].get('symbol'))
                continue
            
            # Try HGNC first
            if self.hgnc is not None and from_type == 'auto':
                if gene_id in self.symbol_to_entrez:
                    results.append(gene_id)
                    self._cache[gene_id] = {'symbol': gene_id}
                    continue
            
            # Query mygene.info
            try:
                result = self.mg.query(gene_id, species='human', fields='symbol')
                
                if result['total'] > 0:
                    symbol = result['hits'][0].get('symbol')
                    results.append(symbol)
                    self._cache[gene_id] = {'symbol': symbol}
                else:
                    results.append(None)
                    self._cache[gene_id] = {}
                    
            except Exception as e:
                logger.warning(f"Failed to map {gene_id}: {e}")
                results.append(None)
        
        return results
    
    def to_entrez(self, symbols: List[str]) -> List[Optional[int]]:
        """
        Convert HGNC symbols to Entrez IDs.
        
        Parameters
        ----------
        symbols : list of str
            HGNC gene symbols.
            
        Returns
        -------
        list of int or None
            Entrez gene IDs. None for unmapped symbols.
        """
        results = []
        
        for symbol in symbols:
            # Try local HGNC first
            if self.hgnc is not None:
                entrez = self.symbol_to_entrez.get(symbol)
                if entrez is not None:
                    results.append(int(entrez))
                    continue
            
            # Query mygene.info
            try:
                result = self.mg.query(f'symbol:{symbol}', species='human', fields='entrezgene')
                
                if result['total'] > 0:
                    entrez = result['hits'][0].get('entrezgene')
                    results.append(int(entrez) if entrez else None)
                else:
                    results.append(None)
                    
            except Exception as e:
                logger.warning(f"Failed to map {symbol}: {e}")
                results.append(None)
        
        return results
    
    def harmonize_gene_list(
        self,
        genes: List[str],
        remove_unmapped: bool = True,
    ) -> List[str]:
        """
        Harmonize a gene list to HGNC symbols.
        
        Parameters
        ----------
        genes : list of str
            Gene identifiers (mixed types OK).
        remove_unmapped : bool
            Whether to remove genes that couldn't be mapped.
            
        Returns
        -------
        list of str
            Harmonized HGNC symbols.
        """
        mapped = self.to_hgnc(genes)
        
        if remove_unmapped:
            mapped = [g for g in mapped if g is not None]
        
        # Remove duplicates while preserving order
        seen = set()
        result = []
        for gene in mapped:
            if gene not in seen and gene is not None:
                seen.add(gene)
                result.append(gene)
        
        n_unmapped = len(genes) - len(result)
        if n_unmapped > 0:
            logger.warning(f"Could not map {n_unmapped}/{len(genes)} genes")
        
        return result
    
    def batch_convert(
        self,
        ids: List[str],
        from_type: str,
        to_type: str = 'symbol',
    ) -> Dict[str, str]:
        """
        Batch convert between ID types.
        
        Parameters
        ----------
        ids : list of str
            Source IDs.
        from_type : str
            Source type ('ensembl.gene', 'entrezgene', 'uniprot', etc.).
        to_type : str
            Target type ('symbol', 'entrezgene', etc.).
            
        Returns
        -------
        dict
            Mapping from source IDs to target IDs.
        """
        results = self.mg.querymany(
            ids,
            scopes=from_type,
            fields=to_type,
            species='human',
            returnall=True,
        )
        
        mapping = {}
        for hit in results['out']:
            if 'notfound' not in hit:
                query = hit['query']
                target = hit.get(to_type)
                if target:
                    mapping[query] = target
        
        logger.info(f"Mapped {len(mapping)}/{len(ids)} IDs from {from_type} to {to_type}")
        
        return mapping