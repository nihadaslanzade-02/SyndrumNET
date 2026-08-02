"""
Tests for gene ID harmonisation.

Nothing here touches the network. `IDMapper.mg` is replaced with a recorder
that returns canned `querymany` responses and counts how often it was called,
because the property that matters is exactly that: how many requests a list of
identifiers costs. The pipeline was unrunnable because the answer used to be
"one per identifier".
"""

from typing import List

import pandas as pd
import pytest

from syndrumnet.io.id_mapping import IDMapper


class FakeMyGene:
    """
    Stand-in for `mygene.MyGeneInfo`, recording every batch it is handed.

    `querymany` answers from a symbol table, mimicking the real response
    shape: one dict per hit, `notfound` for misses, and the ability to return
    several rows for one ambiguous query.
    """

    def __init__(
        self, symbols=None, entrez=None, extra_hits=None, fail=False, omit=(),
    ):
        self.symbols = symbols or {}
        self.entrez = entrez or {}
        self.extra_hits = extra_hits or []
        self.fail = fail
        self.omit = set(omit)
        self.calls: List[List[str]] = []

    def querymany(self, queries, scopes=None, fields=None, species=None, returnall=False):
        self.calls.append(list(queries))

        if self.fail:
            raise ConnectionError("mygene.info unreachable")

        table = self.entrez if fields == 'entrezgene' else self.symbols
        key = 'entrezgene' if fields == 'entrezgene' else 'symbol'

        hits = []
        for query in queries:
            if query in self.omit:
                # The real service does not always echo back every query.
                continue
            if query in table:
                hits.append({'query': query, key: table[query]})
            else:
                hits.append({'query': query, 'notfound': True})

        return hits + list(self.extra_hits)

    @property
    def n_requests(self) -> int:
        return len(self.calls)


@pytest.fixture
def mapper(tmp_path) -> IDMapper:
    """An IDMapper with no HGNC table and an isolated cache directory."""
    instance = IDMapper(cache_dir=tmp_path / "cache")
    instance.mg = FakeMyGene(symbols={
        'ENSG00000141510': 'TP53',
        'TP53': 'TP53',
        '672': 'BRCA1',
        'P38398': 'BRCA1',
    })
    return instance


@pytest.fixture
def hgnc_file(tmp_path) -> "pd.DataFrame":
    """A minimal HGNC table, the local shortcut's source of truth."""
    path = tmp_path / "hgnc.txt"
    pd.DataFrame({
        'symbol': ['TP53', 'BRCA1', 'EGFR'],
        'entrez_id': [7157, 672, 1956],
    }).to_csv(path, sep='\t', index=False)
    return path


# --------------------------------------------------------------------------
# Batching
# --------------------------------------------------------------------------

def test_a_list_of_ids_costs_one_request(mapper):
    """
    The whole point of the change. to_hgnc used to loop and call
    mg.query once per identifier; on an interactome that is tens of
    thousands of sequential round trips and is why the data build was never
    run.
    """
    ids = ['ENSG00000141510', '672', 'P38398', 'NOT_A_GENE']

    symbols = mapper.to_hgnc(ids)

    assert symbols == ['TP53', 'BRCA1', 'BRCA1', None]
    assert mapper.mg.n_requests == 1


def test_duplicates_are_queried_once(mapper):
    """A merged interactome repeats the same gene across sources."""
    mapper.to_hgnc(['672', '672', 'ENSG00000141510', '672'])

    assert mapper.mg.n_requests == 1
    assert mapper.mg.calls[0] == ['672', 'ENSG00000141510']


def test_results_are_aligned_with_the_input(mapper):
    """
    Order and length must survive deduplication and sorting, since the
    caller zips the result against the original list to build a mapping.
    """
    ids = ['NOT_A_GENE', '672', 'ENSG00000141510', 'NOT_A_GENE', '672']

    symbols = mapper.to_hgnc(ids)

    assert len(symbols) == len(ids)
    assert symbols == [None, 'BRCA1', 'TP53', None, 'BRCA1']


def test_a_second_call_does_not_query_again(mapper):
    mapper.to_hgnc(['672', 'ENSG00000141510'])
    mapper.to_hgnc(['672', 'ENSG00000141510'])

    assert mapper.mg.n_requests == 1


def test_misses_are_cached_too(mapper):
    """A miss is a resolved answer; re-querying it every run is waste."""
    mapper.to_hgnc(['NOT_A_GENE'])
    mapper.to_hgnc(['NOT_A_GENE'])

    assert mapper.mg.n_requests == 1


def test_queries_absent_from_the_response_are_still_resolved(tmp_path):
    """
    querymany does not guarantee a row per query. An identifier that simply
    never comes back must be recorded as a miss, or it is re-queried on every
    call and the mapping loop reads a key that was never written.
    """
    instance = IDMapper(cache_dir=tmp_path / "cache")
    instance.mg = FakeMyGene(symbols={'672': 'BRCA1'}, omit={'GHOST'})

    assert instance.to_hgnc(['672', 'GHOST']) == ['BRCA1', None]

    instance.to_hgnc(['GHOST'])

    assert instance.mg.n_requests == 1


def test_only_the_new_identifiers_are_queried(mapper):
    mapper.to_hgnc(['672'])
    mapper.to_hgnc(['672', 'ENSG00000141510'])

    assert mapper.mg.calls == [['672'], ['ENSG00000141510']]


def test_the_first_hit_wins_for_an_ambiguous_query(tmp_path):
    """
    querymany returns several rows for an ambiguous identifier, best score
    first. Taking the last would silently prefer the worst match.
    """
    instance = IDMapper(cache_dir=tmp_path / "cache")
    instance.mg = FakeMyGene(
        symbols={'AMBIG': 'BEST'},
        extra_hits=[{'query': 'AMBIG', 'symbol': 'WORSE'}],
    )

    assert instance.to_hgnc(['AMBIG']) == ['BEST']


def test_a_failed_batch_is_not_cached_as_not_found(tmp_path):
    """
    Caching a network failure as a miss would make every retry return None
    without ever going back to the server.
    """
    instance = IDMapper(cache_dir=tmp_path / "cache")
    instance.mg = FakeMyGene(fail=True)

    assert instance.to_hgnc(['672']) == [None]

    instance.mg = FakeMyGene(symbols={'672': 'BRCA1'})

    assert instance.to_hgnc(['672']) == ['BRCA1']


# --------------------------------------------------------------------------
# Local shortcut and validation
# --------------------------------------------------------------------------

def test_known_symbols_never_reach_the_network(tmp_path, hgnc_file):
    instance = IDMapper(hgnc_file=hgnc_file, cache_dir=tmp_path / "cache")
    instance.mg = FakeMyGene()

    assert instance.to_hgnc(['TP53', 'BRCA1', 'EGFR']) == ['TP53', 'BRCA1', 'EGFR']
    assert instance.mg.n_requests == 0


def test_the_shortcut_does_not_apply_to_other_id_types(tmp_path, hgnc_file):
    """
    'TP53' is a symbol, but under from_type='entrez' it is not an Entrez ID
    and must not be short-circuited to itself.
    """
    instance = IDMapper(hgnc_file=hgnc_file, cache_dir=tmp_path / "cache")
    instance.mg = FakeMyGene()

    assert instance.to_hgnc(['TP53'], from_type='entrez') == [None]
    assert instance.mg.n_requests == 1


def test_unknown_from_type_is_rejected(mapper):
    """
    Every value other than 'auto' used to be accepted and then ignored: the
    query never passed a scope at all, so a typo silently searched
    everything.
    """
    with pytest.raises(ValueError, match="Unknown from_type"):
        mapper.to_hgnc(['TP53'], from_type='ensemble')  # note the typo


def test_missing_hgnc_table_is_not_an_error(tmp_path):
    """symbol_to_entrez must exist even when no HGNC file was supplied."""
    instance = IDMapper(cache_dir=tmp_path / "cache")

    assert instance.symbol_to_entrez == {}


# --------------------------------------------------------------------------
# Disk cache
# --------------------------------------------------------------------------

def test_the_cache_survives_a_new_mapper(tmp_path):
    """
    cache_dir was created in __init__ and then never written to. Persisting
    means the second run of the data build pays nothing for ID mapping.
    """
    cache_dir = tmp_path / "cache"

    first = IDMapper(cache_dir=cache_dir)
    first.mg = FakeMyGene(symbols={'672': 'BRCA1'})
    first.to_hgnc(['672', 'NOT_A_GENE'])

    assert first.cache_file.exists()

    second = IDMapper(cache_dir=cache_dir)
    second.mg = FakeMyGene(symbols={'672': 'BRCA1'})

    assert second.to_hgnc(['672', 'NOT_A_GENE']) == ['BRCA1', None]
    assert second.mg.n_requests == 0


def test_a_corrupt_cache_is_ignored_not_fatal(tmp_path):
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir(parents=True)
    (cache_dir / "hgnc_symbols.json").write_text("{not json")

    instance = IDMapper(cache_dir=cache_dir)
    instance.mg = FakeMyGene(symbols={'672': 'BRCA1'})

    assert instance.to_hgnc(['672']) == ['BRCA1']


def test_disk_cache_can_be_turned_off(tmp_path):
    instance = IDMapper(cache_dir=tmp_path / "cache", use_disk_cache=False)
    instance.mg = FakeMyGene(symbols={'672': 'BRCA1'})
    instance.to_hgnc(['672'])

    assert not instance.cache_file.exists()


# --------------------------------------------------------------------------
# to_entrez and harmonize_gene_list
# --------------------------------------------------------------------------

def test_to_entrez_answers_locally_when_it_can(tmp_path, hgnc_file):
    instance = IDMapper(hgnc_file=hgnc_file, cache_dir=tmp_path / "cache")
    instance.mg = FakeMyGene()

    assert instance.to_entrez(['TP53', 'BRCA1']) == [7157, 672]
    assert instance.mg.n_requests == 0


def test_to_entrez_batches_the_rest(tmp_path, hgnc_file):
    instance = IDMapper(hgnc_file=hgnc_file, cache_dir=tmp_path / "cache")
    instance.mg = FakeMyGene(entrez={'KRAS': 3845, 'MYC': 4609})

    assert instance.to_entrez(['TP53', 'KRAS', 'MYC', 'NOPE']) == [7157, 3845, 4609, None]
    assert instance.mg.n_requests == 1


def test_harmonize_deduplicates_and_preserves_order(mapper):
    result = mapper.harmonize_gene_list(['672', 'ENSG00000141510', 'P38398'])

    assert result == ['BRCA1', 'TP53']


def test_harmonize_counts_failures_not_duplicates(mapper, caplog):
    """
    The unmapped count was len(input) - len(output), so deduplicating two
    copies of a gene that mapped fine was reported as a mapping failure.
    """
    with caplog.at_level("WARNING"):
        mapper.harmonize_gene_list(['672', '672', '672'])

    assert "Could not map" not in caplog.text


def test_harmonize_does_report_real_failures(mapper, caplog):
    with caplog.at_level("WARNING"):
        mapper.harmonize_gene_list(['672', 'NOT_A_GENE'])

    assert "Could not map 1/2" in caplog.text
