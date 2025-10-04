"""
Benchmark data loading for known synergistic combinations.

Loads gold-standard synergy data from resources like DrugCombDB.
"""

import logging
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd

logger = logging.getLogger(__name__)


def load_known_synergies(
    filepath: Path,
    disease_filter: Optional[str] = None,
) -> Set[Tuple[str, str]]:
    """
    Load known synergistic drug combinations.
    
    Parameters
    ----------
    filepath : Path
        Path to synergy data file (CSV format).
    disease_filter : str, optional
        Filter to specific disease.
        
    Returns
    -------
    set of tuple
        Set of (drug_a, drug_b) pairs with known synergy.
    """
    logger.info(f"Loading known synergies from {filepath}")
    
    if not filepath.exists():
        logger.warning(f"Synergy file not found: {filepath}")
        return set()
    
    df = pd.read_csv(filepath)
    
    # Filter by disease if specified
    if disease_filter and 'disease' in df.columns:
        df = df[df['disease'].str.lower() == disease_filter.lower()]
    
    # Extract drug pairs
    synergies = set()
    for _, row in df.iterrows():
        drug_a = str(row['drug_a']).strip()
        drug_b = str(row['drug_b']).strip()
        
        # Canonical ordering (alphabetical)
        pair = tuple(sorted([drug_a, drug_b]))
        synergies.add(pair)
    
    logger.info(f"Loaded {len(synergies)} known synergistic pairs")
    
    return synergies


def load_drugcombdb(filepath: Path) -> pd.DataFrame:
    """
    Load DrugCombDB database.
    
    Parameters
    ----------
    filepath : Path
        Path to DrugCombDB file.
        
    Returns
    -------
    pd.DataFrame
        DrugCombDB data.
    """
    logger.info(f"Loading DrugCombDB from {filepath}")
    
    df = pd.read_csv(filepath)
    
    # Standardize column names
    df = df.rename(columns={
        'Drug1': 'drug_a',
        'Drug2': 'drug_b',
        'Disease': 'disease',
        'Synergy_Score': 'synergy_score',
    })
    
    logger.info(f"Loaded {len(df)} drug combinations from DrugCombDB")
    
    return df
