"""
Parsers for all data sources.

Each parser converts raw data files to standardized pandas DataFrames.

Three of the seven interaction sources in `configs/default.yaml` have a parser
here: HuRI, CORUM and PhosphoSitePlus. KEGG RPair is a stub that returns an
empty frame, and SignaLink, InnateDB and Instruct have no parser at all, so
the network is built from three sources rather than seven. See
`docs/PLACEHOLDERS.md` section 1.
"""

import logging
from pathlib import Path
from typing import Dict, List, Sequence

import pandas as pd

logger = logging.getLogger(__name__)

#: Column names that have carried the signature identifier across LINCS and
#: L1000CDS2 releases, most specific first. The metadata schema is not stable
#: between releases, so the column is detected rather than assumed.
SIGNATURE_ID_COLUMNS = ('sig_id', 'signature_id', 'distil_id', 'pert_sig_id', 'id')

#: Column names that have carried the compound's human-readable name. This is
#: what has to reach the drug module keys: `eval/benchmarks.py` matches
#: predictions against synergy resources on plain compound names, so a
#: `pert_id` (BRD-Kxxxxxxx) joins to nothing and is the last resort.
COMPOUND_NAME_COLUMNS = (
    'pert_iname', 'pert_desc', 'pert_name', 'drug_name', 'compound', 'pert_id',
)


def _detect_column(
    df: pd.DataFrame,
    candidates: Sequence[str],
    description: str,
    source: Path,
) -> str:
    """
    Find the first candidate column present, or fail naming what was there.

    Returns the column name. Raises `KeyError` listing both the candidates
    tried and the columns actually found, because a wrong guess here would
    otherwise surface much later as an empty join.
    """
    for candidate in candidates:
        if candidate in df.columns:
            return candidate

    raise KeyError(
        f"Could not find the {description} column in {source}. "
        f"Tried {list(candidates)}; the file has {list(df.columns)}. "
        f"Pass the right name explicitly if this release uses another."
    )


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

    Placeholder: returns an empty frame with the right schema, so a caller
    that concatenates it gets fewer edges rather than an error.

    Finishing it needs two things. The REST endpoint is already listed in
    `DataDownloader.URLS` but no `download_*` method calls it, and reactant
    pairs are compound-to-compound, so they only become interactome edges via
    the enzymes that catalyse the reaction. That mapping step is the work.

    Note that `NetworkBuilder.add_source` does not register this parser, so
    it is currently unreachable from the pipeline even with a file in hand.

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
    aggregate: str = 'median',
) -> Dict[str, Dict[str, List[str]]]:
    """
    Parse LINCS L1000 drug signatures.

    Defines drug modules as top/bottom 5% genes by fold-change, from a profile
    aggregated per compound across every cell line it was assayed in.

    Parameters
    ----------
    sig_filepath : Path
        Path to LINCS signatures file. Genes on the rows, one column per
        signature.
    meta_filepath : Path
        Path to LINCS metadata file, carrying the signature-to-compound
        mapping.
    top_pct : float
        Percentage of genes to include in up/down sets (default: 0.05 = 5%).
    aggregate : str
        How to combine a compound's signatures per gene, 'median' or 'mean'.
        Median by default, so one anomalous cell line cannot carry a gene into
        the module on its own.

    Returns
    -------
    dict
        {compound_name: {'up': [genes], 'down': [genes]}}

    Raises
    ------
    KeyError
        If the signature-ID or compound-name column cannot be found in the
        metadata. Both are detected from a list of names used across
        releases, and the error reports what the file actually contains.
    ValueError
        If `aggregate` is unknown, if `top_pct` is outside (0, 1], or if no
        signature column matches any metadata row, which means the detected
        ID column is the wrong one.

    Notes
    -----
    A compound is assayed in several cell lines, and each assay is its own
    column. Without the metadata join the keys here are raw signature
    identifiers, which had two consequences: one compound became several
    distinct "drugs" whose pairs are guaranteed to score as redundant, and
    the output could not be evaluated at all, because
    `eval/benchmarks.py` matches predictions against synergy resources keyed
    by plain compound name.

    The metadata schema is not stable across LINCS and L1000CDS2 releases, so
    the two columns are detected from `SIGNATURE_ID_COLUMNS` and
    `COMPOUND_NAME_COLUMNS` rather than assumed. `pert_id` is the last
    candidate for the name, since a BRD identifier joins to no synergy
    resource; if that is what gets picked, expect evaluation to find nothing.

    Signature columns with no metadata row are dropped with a warning. Keeping
    them under their raw ID would silently reintroduce exactly the problem
    this join exists to remove.
    """
    if aggregate not in ('median', 'mean'):
        raise ValueError(
            f"Unknown aggregate: {aggregate!r}. Expected 'median' or 'mean'."
        )

    if not 0 < top_pct <= 1:
        raise ValueError(f"top_pct must be in (0, 1], got {top_pct}.")

    logger.info(f"Parsing LINCS L1000 from {sig_filepath}")

    # Read signatures: genes on the rows, one column per signature
    df = pd.read_csv(sig_filepath, sep='\t', index_col=0)

    # Read metadata and resolve each signature column to its compound
    meta = pd.read_csv(meta_filepath, sep='\t')

    sig_column = _detect_column(
        meta, SIGNATURE_ID_COLUMNS, "signature identifier", meta_filepath
    )
    name_column = _detect_column(
        meta, COMPOUND_NAME_COLUMNS, "compound name", meta_filepath
    )

    logger.info(
        f"Joining on {sig_column!r}, taking compound names from {name_column!r}"
    )

    # Last row wins for a repeated signature ID; they should be unique, and a
    # duplicated index would make the reindex below raise.
    compound_of = (
        meta.drop_duplicates(subset=sig_column, keep='last')
        .set_index(sig_column)[name_column]
        .reindex(df.columns)
    )

    unmapped = compound_of.isna()

    if unmapped.all():
        raise ValueError(
            f"No signature column in {sig_filepath} matched a {sig_column!r} "
            f"value in {meta_filepath}. The detected join column is probably "
            f"the wrong one."
        )

    if unmapped.any():
        logger.warning(
            f"Dropping {int(unmapped.sum())}/{len(compound_of)} signatures "
            f"with no metadata row"
        )
        df = df.loc[:, ~unmapped.to_numpy()]
        compound_of = compound_of[~unmapped]

    # Aggregate each compound's signatures per gene. Transposed because
    # grouping along the column axis was removed in pandas 3.
    profiles = df.T.groupby(compound_of.to_numpy()).agg(aggregate).T

    logger.info(
        f"Aggregated {df.shape[1]} signatures into {profiles.shape[1]} "
        f"compound profiles by {aggregate}"
    )

    signatures = {}

    for drug in profiles.columns:
        fold_changes = profiles[drug].dropna()

        # Top/bottom percentiles. At least one gene each way, so a small
        # matrix yields a usable module instead of an empty one.
        n_top = max(1, int(len(fold_changes) * top_pct))

        signatures[drug] = {
            'up': fold_changes.nlargest(n_top).index.tolist(),
            'down': fold_changes.nsmallest(n_top).index.tolist(),
        }

    logger.info(f"Parsed LINCS signatures for {len(signatures)} compounds")
    return signatures
