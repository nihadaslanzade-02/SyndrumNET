"""
Integrated human molecular interaction network construction.

Combines multiple data sources (PPI, complexes, pathways, signaling) into
a unified network graph for propagation and distance calculations.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import networkx as nx
import pandas as pd

from syndrumnet.io.parsers import (
    parse_huri,
    parse_corum,
    parse_phosphositeplus,
)
from syndrumnet.io.id_mapping import IDMapper

logger = logging.getLogger(__name__)


class NetworkBuilder:
    """
    Build integrated human molecular interaction network.
    
    Combines PPI, protein complexes, kinase-substrate, and pathway data
    into a unified undirected graph for network-based analysis.
    
    Parameters
    ----------
    id_mapper : IDMapper
        Gene ID mapping service.
    filter_organism : str
        Filter to specific organism (default: 'human').
    min_confidence : float, optional
        Minimum confidence score for interactions (source-dependent).
        
    Examples
    --------
    >>> mapper = IDMapper()
    >>> builder = NetworkBuilder(mapper)
    >>> builder.add_source('huri', 'data/raw/huri.tsv')
    >>> builder.add_source('corum', 'data/raw/corum.txt')
    >>> G = builder.build()
    >>> print(f"Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    """
    
    def __init__(
        self,
        id_mapper: IDMapper,
        filter_organism: str = 'human',
        min_confidence: Optional[float] = None,
    ) -> None:
        """Initialize network builder."""
        self.id_mapper = id_mapper
        self.filter_organism = filter_organism
        self.min_confidence = min_confidence
        
        # Store interactions from all sources
        self.interactions: List[pd.DataFrame] = []
        self.source_counts: Dict[str, int] = {}
        
        # Final network
        self.network: Optional[nx.Graph] = None
    
    def add_source(
        self,
        source_name: str,
        filepath: Path,
        **parser_kwargs,
    ) -> None:
        """
        Add interactions from a data source.
        
        Parameters
        ----------
        source_name : str
            Source identifier ('huri', 'corum', etc.).
        filepath : Path
            Path to source data file.
        **parser_kwargs
            Additional arguments for source-specific parser.
        """
        logger.info(f"Adding interactions from {source_name}")
        
        # Parse based on source
        parsers = {
            'huri': parse_huri,
            'corum': parse_corum,
            'phosphositeplus': parse_phosphositeplus,
        }
        
        if source_name not in parsers:
            logger.warning(f"No parser for {source_name}, skipping")
            return
        
        parser = parsers[source_name]
        df = parser(filepath, **parser_kwargs)
        
        if df is not None and not df.empty:
            self.interactions.append(df)
            self.source_counts[source_name] = len(df)
            logger.info(f"Added {len(df)} interactions from {source_name}")
        else:
            logger.warning(f"No interactions loaded from {source_name}")
    
    def build(self, deduplicate: bool = True) -> nx.Graph:
        """
        Build integrated network from all sources.
        
        Parameters
        ----------
        deduplicate : bool
            Remove duplicate edges across sources.
            
        Returns
        -------
        nx.Graph
            Integrated molecular interaction network.
        """
        logger.info("Building integrated network")
        
        if not self.interactions:
            raise ValueError("No interaction sources added")
        
        # Combine all interactions
        all_interactions = pd.concat(self.interactions, ignore_index=True)
        
        logger.info(f"Total interactions before harmonization: {len(all_interactions)}")
        
        # Harmonize gene IDs
        all_genes = set(all_interactions['protein_a']) | set(all_interactions['protein_b'])
        all_genes = list(all_genes)
        
        logger.info(f"Harmonizing {len(all_genes)} unique genes")
        harmonized = self.id_mapper.harmonize_gene_list(all_genes)
        
        gene_map = dict(zip(all_genes, self.id_mapper.to_hgnc(all_genes)))
        
        # Map to harmonized IDs
        all_interactions['gene_a'] = all_interactions['protein_a'].map(gene_map)
        all_interactions['gene_b'] = all_interactions['protein_b'].map(gene_map)
        
        # Remove unmapped
        all_interactions = all_interactions[
            all_interactions['gene_a'].notna() &
            all_interactions['gene_b'].notna()
        ]
        
        logger.info(f"Interactions after ID mapping: {len(all_interactions)}")
        
        # Build NetworkX graph
        G = nx.Graph()
        
        for _, row in all_interactions.iterrows():
            gene_a = row['gene_a']
            gene_b = row['gene_b']
            
            # Skip self-loops
            if gene_a == gene_b:
                continue
            
            # Add edge with metadata
            if G.has_edge(gene_a, gene_b):
                # Edge exists, add source to metadata
                G[gene_a][gene_b]['sources'].append(row['source'])
            else:
                G.add_edge(
                    gene_a,
                    gene_b,
                    sources=[row['source']],
                    interaction_type=row.get('interaction_type', 'PPI'),
                )
        
        # Get largest connected component
        if not nx.is_connected(G):
            logger.warning("Network is not connected, taking largest component")
            components = list(nx.connected_components(G))
            largest = max(components, key=len)
            G = G.subgraph(largest).copy()
            logger.info(f"Largest component: {len(largest)} nodes")
        
        self.network = G
        
        logger.info(
            f"Built network: {G.number_of_nodes()} nodes, "
            f"{G.number_of_edges()} edges"
        )
        
        return G
    
    def get_network_stats(self) -> Dict[str, any]:
        """
        Get network statistics.
        
        Returns
        -------
        dict
            Network properties (nodes, edges, density, etc.).
        """
        if self.network is None:
            raise ValueError("Network not built yet")
        
        G = self.network
        
        stats = {
            'n_nodes': G.number_of_nodes(),
            'n_edges': G.number_of_edges(),
            'density': nx.density(G),
            'avg_degree': sum(dict(G.degree()).values()) / G.number_of_nodes(),
            'n_components': nx.number_connected_components(G),
            'diameter': nx.diameter(G) if nx.is_connected(G) else None,
            'avg_clustering': nx.average_clustering(G),
        }
        
        return stats
    
    def save(self, output_path: Path) -> None:
        """
        Save network to file.
        
        Parameters
        ----------
        output_path : Path
            Output file path (.graphml, .gml, .edgelist).
        """
        if self.network is None:
            raise ValueError("Network not built yet")
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Save based on extension
        ext = output_path.suffix.lower()
        
        if ext == '.graphml':
            nx.write_graphml(self.network, output_path)
        elif ext == '.gml':
            nx.write_gml(self.network, output_path)
        elif ext == '.edgelist':
            nx.write_edgelist(self.network, output_path)
        else:
            raise ValueError(f"Unsupported format: {ext}")
        
        logger.info(f"Saved network to {output_path}")
    
    @staticmethod
    def load(filepath: Path) -> nx.Graph:
        """
        Load network from file.
        
        Parameters
        ----------
        filepath : Path
            Network file path.
            
        Returns
        -------
        nx.Graph
            Loaded network.
        """
        ext = filepath.suffix.lower()
        
        if ext == '.graphml':
            G = nx.read_graphml(filepath)
        elif ext == '.gml':
            G = nx.read_gml(filepath)
        elif ext == '.edgelist':
            G = nx.read_edgelist(filepath)
        else:
            raise ValueError(f"Unsupported format: {ext}")
        
        logger.info(f"Loaded network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
        
        return G
