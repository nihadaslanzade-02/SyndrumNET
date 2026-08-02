"""
Gene/protein ID mapping and harmonization.

Converts between different identifier systems (Entrez, Ensembl, HGNC, UniProt)
to ensure consistent gene names across all data sources.

Lookups are batched. An interactome carries 15,000 to 20,000 distinct
identifiers, and resolving them one HTTP request at a time is what kept this
pipeline from being run at all; `mygene.info` answers the whole list in a
single `querymany` call. Results are also cached on disk, so the cost is paid
once per machine rather than once per run.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import mygene
import pandas as pd

logger = logging.getLogger(__name__)

#: `mygene` scopes searched for each `from_type`. 'auto' searches every
#: identifier system this pipeline encounters, which is what makes a single
#: batch possible over the mixed identifiers a merged interactome contains.
ID_SCOPES = {
    'auto': 'symbol,entrezgene,ensembl.gene,uniprot,alias,retired',
    'symbol': 'symbol',
    'entrez': 'entrezgene',
    'ensembl': 'ensembl.gene',
    'uniprot': 'uniprot',
}

#: Filename of the on-disk symbol cache inside `cache_dir`.
SYMBOL_CACHE_FILE = "hgnc_symbols.json"


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
        Directory for ID mapping cache. Resolved symbols are written here and
        reloaded on the next run.
    use_disk_cache : bool
        Read and write the on-disk cache. Set False to force fresh lookups.

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
        use_disk_cache: bool = True,
    ) -> None:
        """Initialize ID mapper."""
        self.cache_dir = Path(cache_dir) if cache_dir else Path("data/interim/id_cache")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.use_disk_cache = use_disk_cache

        # Load HGNC data
        self.hgnc = None
        self.symbol_to_entrez: Dict[str, int] = {}
        self.entrez_to_symbol: Dict[int, str] = {}

        if hgnc_file and hgnc_file.exists():
            self._load_hgnc(hgnc_file)

        # Initialize mygene client
        self.mg = mygene.MyGeneInfo()

        # Cache for conversions, keyed by identifier. An empty dict records a
        # resolved-and-not-found identifier, so it is not queried again.
        self._cache: Dict[str, Dict] = {}

        if self.use_disk_cache:
            self.load_cache()

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

    # ----------------------------------------------------------------------
    # Disk cache
    # ----------------------------------------------------------------------

    @property
    def cache_file(self) -> Path:
        """Path of the on-disk symbol cache."""
        return self.cache_dir / SYMBOL_CACHE_FILE

    def load_cache(self) -> None:
        """
        Load previously resolved symbols from disk.

        A corrupt or unreadable cache is a performance problem, not a
        correctness one, so it is logged and ignored rather than raised.
        """
        if not self.cache_file.exists():
            return

        try:
            with open(self.cache_file, 'r', encoding='utf-8') as f:
                stored = json.load(f)
        except (OSError, json.JSONDecodeError) as e:
            logger.warning(f"Ignoring unreadable ID cache {self.cache_file}: {e}")
            return

        self._cache.update({
            key: ({'symbol': value} if value else {})
            for key, value in stored.items()
        })

        logger.info(f"Loaded {len(stored)} cached identifiers from {self.cache_file}")

    def save_cache(self) -> None:
        """Write resolved symbols to disk for the next run."""
        stored = {key: value.get('symbol') for key, value in self._cache.items()}

        try:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            with open(self.cache_file, 'w', encoding='utf-8') as f:
                json.dump(stored, f)
        except OSError as e:
            logger.warning(f"Could not write ID cache {self.cache_file}: {e}")
            return

        logger.debug(f"Cached {len(stored)} identifiers in {self.cache_file}")

    # ----------------------------------------------------------------------
    # Conversions
    # ----------------------------------------------------------------------

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
            Source ID type ('ensembl', 'entrez', 'uniprot', 'symbol',
            'auto'). 'auto' searches every scope at once.

        Returns
        -------
        list of str or None
            HGNC symbols, aligned with `ids`. None for unmapped IDs.

        Raises
        ------
        ValueError
            If `from_type` is not a known scope. This used to be accepted and
            then ignored: every value other than 'auto' behaved identically,
            because the query never passed a scope at all.

        Notes
        -----
        Everything not already known locally goes to `mygene.info` in one
        `querymany` call. Identifiers already carried in a loaded HGNC table
        never leave the machine, and anything resolved once is cached in
        memory and on disk, so a repeat call over an overlapping list costs
        nothing.
        """
        if from_type not in ID_SCOPES:
            raise ValueError(
                f"Unknown from_type: {from_type!r}. "
                f"Expected one of {sorted(ID_SCOPES)}."
            )

        unresolved = []

        for gene_id in ids:
            if gene_id in self._cache:
                continue

            # Already an HGNC symbol, no lookup needed.
            if from_type in ('auto', 'symbol') and gene_id in self.symbol_to_entrez:
                self._cache[gene_id] = {'symbol': gene_id}
                continue

            unresolved.append(gene_id)

        if unresolved:
            self._resolve_symbols(unresolved, from_type)

        return [self._cache.get(gene_id, {}).get('symbol') for gene_id in ids]

    def _resolve_symbols(self, ids: List[str], from_type: str) -> None:
        """
        Resolve identifiers to symbols in a single batched query.

        Deduplicates first, since a merged interactome repeats the same gene
        across sources. `querymany` can return several rows for one query
        when an identifier is ambiguous; they arrive best-scoring first, so
        the first symbol seen wins.
        """
        queries = sorted(set(ids))

        logger.info(
            f"Resolving {len(queries)} identifiers via mygene.info "
            f"in one batch (scopes: {ID_SCOPES[from_type]})"
        )

        try:
            hits = self.mg.querymany(
                queries,
                scopes=ID_SCOPES[from_type],
                fields='symbol',
                species='human',
                returnall=False,
            )
        except Exception as e:
            # A failed batch must not be cached as "not found", or a retry
            # would silently return None for everything.
            logger.warning(f"Batch ID lookup failed for {len(queries)} identifiers: {e}")
            return

        for hit in hits:
            query = hit.get('query')

            if query is None or self._cache.get(query, {}).get('symbol'):
                continue

            symbol = None if hit.get('notfound') else hit.get('symbol')
            self._cache[query] = {'symbol': symbol} if symbol else {}

        # Record the misses too, so they are not re-queried.
        for gene_id in queries:
            self._cache.setdefault(gene_id, {})

        n_mapped = sum(1 for q in queries if self._cache[q].get('symbol'))
        logger.info(f"Resolved {n_mapped}/{len(queries)} identifiers")

        if self.use_disk_cache:
            self.save_cache()

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
            Entrez gene IDs, aligned with `symbols`. None for unmapped ones.

        Notes
        -----
        Batched for the same reason as `to_hgnc`. Symbols present in a loaded
        HGNC table are answered locally and never reach the network.
        """
        resolved: Dict[str, Optional[int]] = {}
        unresolved = []

        for symbol in symbols:
            entrez = self.symbol_to_entrez.get(symbol)

            if entrez is not None and not pd.isna(entrez):
                resolved[symbol] = int(entrez)
            else:
                unresolved.append(symbol)

        if unresolved:
            resolved.update(self._resolve_entrez(sorted(set(unresolved))))

        return [resolved.get(symbol) for symbol in symbols]

    def _resolve_entrez(self, symbols: List[str]) -> Dict[str, Optional[int]]:
        """Resolve symbols to Entrez IDs in a single batched query."""
        logger.info(f"Resolving {len(symbols)} symbols to Entrez in one batch")

        try:
            hits = self.mg.querymany(
                symbols,
                scopes='symbol',
                fields='entrezgene',
                species='human',
                returnall=False,
            )
        except Exception as e:
            logger.warning(f"Batch Entrez lookup failed for {len(symbols)} symbols: {e}")
            return {}

        resolved: Dict[str, Optional[int]] = {}

        for hit in hits:
            query = hit.get('query')

            if query is None or resolved.get(query) is not None:
                continue

            entrez = None if hit.get('notfound') else hit.get('entrezgene')
            resolved[query] = int(entrez) if entrez else None

        return resolved

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
            Harmonized HGNC symbols, deduplicated, in first-seen order.

        Notes
        -----
        Deduplication means the returned list is generally shorter than the
        input even when everything maps, so its length says nothing about how
        many identifiers failed. The unmapped count logged below counts the
        failures themselves.
        """
        mapped = self.to_hgnc(genes)

        n_unmapped = sum(1 for symbol in mapped if symbol is None)

        if not remove_unmapped:
            n_unmapped = 0

        # Remove duplicates while preserving order
        seen = set()
        result = []

        for gene in mapped:
            if gene is not None and gene not in seen:
                seen.add(gene)
                result.append(gene)

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
            A raw `mygene` scope string, unlike `to_hgnc`'s `from_type`.
        to_type : str
            Target type ('symbol', 'entrezgene', etc.).

        Returns
        -------
        dict
            Mapping from source IDs to target IDs. Unmapped IDs are absent.
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
