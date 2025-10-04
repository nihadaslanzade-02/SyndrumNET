"""
Parsers for all data sources.

Each parser converts raw data files to standardized pandas DataFrames.
"""

import gzip
import logging
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd

logger = logging.getLogger(__name__)


def parse_huri(filepath: Path) -> pd.DataFrame:
    """
    Parse HuRI protein-protein interaction data.
    
    Parameters
    ----------
    filepath : Path
        Path to HuRI TSV file.
        
    Returns
    -------
    pd.DataFrame
        Columns: ['protein_a', 'protein_b', 'source']
    """
    logger.info(f"Parsing HuRI from {filepath}")
    
    df = pd.read_csv(filepath, sep='\t', comment='#')
    
    # Standardize column names
    df = df.rename(columns={
        'Symbol A': 'protein_a',
        'Symbol B': 'protein_b',
    })
    
    df['source'] = 'HuRI'
    df['interaction_type'] = 'PPI'
    
    # Keep only human, filter by confidence if available
    df = df[df['protein_a'].notna() & df['protein_b'].notna()]
    
    logger.info(f"Parsed {len(df)} HuRI interactions")
    return df[['protein_a', 'protein_b', 'source', 'interaction_type']]


def parse_corum(filepath: Path) -> pd.DataFrame:
    """
    Parse CORUM protein complex data.
    
    Converts complexes to pairwise interactions (all-vs-all within complex).
    
    Parameters
    ----------
    filepath : Path
        Path to CORUM allComplexes.txt file.
        
    Returns
    -------
    pd.DataFrame
        Columns: ['protein_a', 'protein_b', 'source', 'complex_id']
    """
    logger.info(f"Parsing CORUM from {filepath}")
    
    df = pd.read_csv(filepath, sep='\t')
    
    # Filter for human
    df = df[df['Organism'] == 'Human']
    
    interactions = []
    
    for _, row in df.iterrows():
        # Get complex members (gene symbols)
        members = str(row['subunits(Gene name)']).split(';')
        members = [m.strip() for m in members if m.strip()]
        
        # Create all pairwise interactions
        for i, prot_a in enumerate(members):
            for prot_b in members[i+1:]:
                interactions.append({
                    'protein_a': prot_a,
                    'protein_b': prot_b,
                    'source': 'CORUM',
                    'interaction_type': 'complex',
                    'complex_id': row['ComplexID'],
                })
    
    result = pd.DataFrame(interactions)
    logger.info(f"Parsed {len(result)} CORUM interactions from {len(df)} complexes")
    
    return result


def parse_phosphositeplus(filepath: Path) -> pd.DataFrame:
    """
    Parse PhosphoSitePlus kinase-substrate data.
    
    Parameters
    ----------
    filepath : Path
        Path to PhosphoSitePlus file.
        
    Returns
    -------
    pd.DataFrame
        Columns: ['protein_a', 'protein_b', 'source', 'interaction_type']
    """
    logger.info(f"Parsing PhosphoSitePlus from {filepath}")
    
    df = pd.read_csv(filepath, sep='\t', skiprows=3)
    
    # Filter for human
    df = df[df['KIN_ORGANISM'] == 'human']
    df = df[df['SUB_ORGANISM'] == 'human']
    
    # Extract gene names
    df = df.rename(columns={
        'KINASE': 'protein_a',  # Kinase
        'SUBSTRATE': 'protein_b',  # Substrate
    })
    
    df['source'] = 'PhosphoSitePlus'
    df['interaction_type'] = 'phosphorylation'
    
    result = df[['protein_a', 'protein_b', 'source', 'interaction_type']].drop_duplicates()
    
    logger.info(f"Parsed {len(result)} PhosphoSitePlus interactions")
    return result


def parse_kegg_rpair(filepath: Path) -> pd.DataFrame:
    """
    Parse KEGG RPair reaction data.
    
    Note: This is a placeholder. Actual implementation requires KEGG API calls.
    
    Parameters
    ----------
    filepath : Path
        Path to KEGG RPair data.
        
    Returns
    -------
    pd.DataFrame
        Columns: ['protein_a', 'protein_b', 'source', 'interaction_type']
    """
    logger.warning("KEGG RPair parsing requires API integration (placeholder)")
    
    # Return empty DataFrame with correct schema
    return pd.DataFrame(columns=['protein_a', 'protein_b', 'source', 'interaction_type'])


def parse_creeds(filepath: Path) -> Dict[str, Dict[str, List[str]]]:
    """
    Parse CREEDS disease signatures.
    
    Returns disease signatures as up/down gene sets.
    
    Parameters
    ----------
    filepath : Path
        Path to CREEDS disease signatures file.
        
    Returns
    -------
    dict
        {disease_name: {'up': [genes], 'down': [genes]}}
    """
    logger.info(f"Parsing CREEDS from {filepath}")
    
    df = pd.read_csv(filepath, sep='\t')
    
    signatures = {}
    
    # Group by disease
    for disease in df['disease_name'].unique():
        disease_df = df[df['disease_name'] == disease]
        
        # Separate up and down regulated genes
        up_genes = disease_df[disease_df['direction'] == 'up']['gene_symbol'].tolist()
        down_genes = disease_df[disease_df['direction'] == 'down']['gene_symbol'].tolist()
        
        signatures[disease] = {
            'up': up_genes,
            'down': down_genes,
        }
    
    logger.info(f"Parsed CREEDS signatures for {len(signatures)} diseases")
    return signatures


def parse_lincs(
    sig_filepath: Path,
    meta_filepath: Path,
    top_pct: float = 0.05,
) -> Dict[str, Dict[str, List[str]]]:
    """
    Parse LINCS L1000 drug signatures.
    
    Defines drug modules as top/bottom 5% genes by fold-change.
    
    Parameters
    ----------
    sig_filepath : Path
        Path to LINCS signatures file.
    meta_filepath : Path
        Path to LINCS metadata file.
    top_pct : float
        Percentage of genes to include in up/down sets (default: 0.05 = 5%).
        
    Returns
    -------
    dict
        {drug_name: {'up': [genes], 'down': [genes]}}
    """
    logger.info(f"Parsing LINCS L1000 from {sig_filepath}")
    
    # Read signatures
    df = pd.read_csv(sig_filepath, sep='\t', index_col=0)
    
    # Read metadata
    meta = pd.read_csv(meta_filepath, sep='\t')
    
    signatures = {}
    
    # For each drug (column in signatures)
    for drug in df.columns:
        fold_changes = df[drug].dropna()
        
        # Top/bottom percentiles
        n_top = int(len(fold_changes) * top_pct)
        
        top_up = fold_changes.nlargest(n_top)
        top_down = fold_changes.nsmallest(n_top)
        
        signatures[drug] = {
            'up': top_up.index.tolist(),
            'down': top_down.index.tolist(),
        }
    
    logger.info(f"Parsed LINCS signatures for {len(signatures)} drugs")
    return signatures
