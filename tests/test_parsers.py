"""
Tests for `parse_lincs`, the metadata join.

The signature matrix is built so the right answer is arithmetic rather than
whatever the code prints: two cell lines per compound with known fold changes,
so the aggregated profile and therefore the module membership follow from the
numbers by hand.
"""

import pandas as pd
import pytest

from syndrumnet.io.parsers import parse_lincs

GENES = [f"G{i}" for i in range(10)]


def write_signatures(path, columns) -> None:
    """Signature matrix: genes on the rows, one column per signature."""
    pd.DataFrame(columns, index=GENES).to_csv(path, sep='\t')


def write_metadata(path, rows, sig_column='sig_id', name_column='pert_iname') -> None:
    pd.DataFrame(
        [{sig_column: sig, name_column: name} for sig, name in rows]
    ).to_csv(path, sep='\t', index=False)


@pytest.fixture
def lincs(tmp_path):
    """
    Two compounds, two cell lines each.

    imatinib's two profiles disagree on G0 and G9, so the median across them
    is what decides whether those genes make the module. Chosen so the
    per-line answer and the aggregated answer differ: taking either column on
    its own puts G0 at the top, the median does not.
    """
    sig_path = tmp_path / "sigs.tsv"
    meta_path = tmp_path / "meta.tsv"

    write_signatures(sig_path, {
        # imatinib, cell line A: G0 extreme up, G1 mildly up
        'SIG_A1': [9.0, 2.0, 0.5, 0.1, 0.0, -0.1, -0.2, -0.5, -1.0, -3.0],
        # imatinib, cell line B: G0 extreme down, cancels in the median
        'SIG_B1': [-9.0, 2.5, 0.4, 0.2, 0.0, -0.1, -0.3, -0.6, -1.1, -2.0],
        # dasatinib, two consistent lines
        'SIG_A2': [0.1, 0.2, 3.0, 0.0, -0.1, -0.2, -0.3, -0.4, -0.5, -4.0],
        'SIG_B2': [0.2, 0.1, 3.4, 0.1, -0.1, -0.2, -0.4, -0.5, -0.6, -4.2],
    })
    write_metadata(meta_path, [
        ('SIG_A1', 'imatinib'),
        ('SIG_B1', 'imatinib'),
        ('SIG_A2', 'dasatinib'),
        ('SIG_B2', 'dasatinib'),
    ])

    return sig_path, meta_path


# --------------------------------------------------------------------------
# The join
# --------------------------------------------------------------------------

def test_modules_are_keyed_by_compound(lincs):
    """
    The keys used to be the signature matrix's column names, so one compound
    became several drugs and nothing downstream could match it to a synergy
    resource.
    """
    signatures = parse_lincs(*lincs, top_pct=0.2)

    assert set(signatures) == {'imatinib', 'dasatinib'}


def test_cell_lines_are_aggregated_not_taken_separately(lincs):
    """
    G0 is +9 in one imatinib line and -9 in the other. Either column on its
    own puts it at an extreme; the median across them is 0, so it belongs in
    neither the up nor the down set.
    """
    signatures = parse_lincs(*lincs, top_pct=0.2)

    imatinib = signatures['imatinib']

    assert 'G0' not in imatinib['up']
    assert 'G0' not in imatinib['down']
    assert imatinib['up'] == ['G1', 'G2']
    assert imatinib['down'] == ['G9', 'G8']


def test_median_resists_an_outlier_line_where_mean_does_not(tmp_path):
    """
    Why median is the default. Three lines: G1 is mildly up in all three,
    G0 is flat in two and wildly up in one. The mean lets that single line
    carry G0 into the module; the median does not.
    """
    sig_path = tmp_path / "sigs.tsv"
    meta_path = tmp_path / "meta.tsv"

    write_signatures(sig_path, {
        'S1': [0.0, 1.0, 0.5, 0.1, 0.0, -0.1, -0.2, -0.5, -1.0, -3.0],
        'S2': [0.0, 1.1, 0.4, 0.2, 0.0, -0.1, -0.3, -0.6, -1.1, -2.0],
        'S3': [30.0, 1.2, 0.4, 0.2, 0.0, -0.1, -0.3, -0.6, -1.1, -2.5],
    })
    write_metadata(meta_path, [
        ('S1', 'imatinib'), ('S2', 'imatinib'), ('S3', 'imatinib'),
    ])

    by_median = parse_lincs(sig_path, meta_path, top_pct=0.1, aggregate='median')
    by_mean = parse_lincs(sig_path, meta_path, top_pct=0.1, aggregate='mean')
    by_default = parse_lincs(sig_path, meta_path, top_pct=0.1)

    assert by_median['imatinib']['up'] == ['G1']
    assert by_mean['imatinib']['up'] == ['G0']

    # The default has to be the robust one, not merely available.
    assert by_default['imatinib']['up'] == ['G1']


def test_an_unknown_aggregate_is_rejected(lincs):
    with pytest.raises(ValueError, match="Unknown aggregate"):
        parse_lincs(*lincs, aggregate='mode')


def test_a_single_signature_compound_still_works(tmp_path):
    sig_path = tmp_path / "sigs.tsv"
    meta_path = tmp_path / "meta.tsv"

    write_signatures(sig_path, {'ONLY': list(range(10))})
    write_metadata(meta_path, [('ONLY', 'nilotinib')])

    signatures = parse_lincs(sig_path, meta_path, top_pct=0.2)

    assert signatures['nilotinib']['up'] == ['G9', 'G8']
    assert signatures['nilotinib']['down'] == ['G0', 'G1']


# --------------------------------------------------------------------------
# Column detection
# --------------------------------------------------------------------------

@pytest.mark.parametrize("sig_column", ['sig_id', 'signature_id', 'distil_id', 'id'])
@pytest.mark.parametrize("name_column", ['pert_iname', 'pert_desc', 'drug_name'])
def test_column_names_are_detected_across_release_schemas(
    tmp_path, sig_column, name_column
):
    """
    The metadata schema is not stable across LINCS and L1000CDS2 releases, so
    the columns are detected rather than assumed.
    """
    sig_path = tmp_path / "sigs.tsv"
    meta_path = tmp_path / "meta.tsv"

    write_signatures(sig_path, {'S1': list(range(10))})
    write_metadata(meta_path, [('S1', 'imatinib')], sig_column, name_column)

    assert set(parse_lincs(sig_path, meta_path, top_pct=0.2)) == {'imatinib'}


def test_an_undetectable_id_column_names_what_it_found(tmp_path):
    sig_path = tmp_path / "sigs.tsv"
    meta_path = tmp_path / "meta.tsv"

    write_signatures(sig_path, {'S1': list(range(10))})
    pd.DataFrame({'weird_key': ['S1'], 'pert_iname': ['imatinib']}).to_csv(
        meta_path, sep='\t', index=False
    )

    with pytest.raises(KeyError, match="weird_key"):
        parse_lincs(sig_path, meta_path)


def test_an_undetectable_name_column_names_what_it_found(tmp_path):
    sig_path = tmp_path / "sigs.tsv"
    meta_path = tmp_path / "meta.tsv"

    write_signatures(sig_path, {'S1': list(range(10))})
    pd.DataFrame({'sig_id': ['S1'], 'weird_name': ['imatinib']}).to_csv(
        meta_path, sep='\t', index=False
    )

    with pytest.raises(KeyError, match="weird_name"):
        parse_lincs(sig_path, meta_path)


# --------------------------------------------------------------------------
# Partial and failed joins
# --------------------------------------------------------------------------

def test_signatures_without_metadata_are_dropped_with_a_warning(tmp_path, caplog):
    """
    Keeping an unmapped column under its raw ID would quietly reintroduce the
    signature-as-drug problem this join exists to remove.
    """
    sig_path = tmp_path / "sigs.tsv"
    meta_path = tmp_path / "meta.tsv"

    write_signatures(sig_path, {'S1': list(range(10)), 'ORPHAN': list(range(10))})
    write_metadata(meta_path, [('S1', 'imatinib')])

    with caplog.at_level("WARNING"):
        signatures = parse_lincs(sig_path, meta_path, top_pct=0.2)

    assert set(signatures) == {'imatinib'}
    assert "1/2 signatures" in caplog.text


def test_a_join_that_matches_nothing_is_an_error(tmp_path):
    """
    Detecting the wrong ID column produces an empty result rather than a
    wrong one, which is exactly the kind of silence this repository had too
    much of.
    """
    sig_path = tmp_path / "sigs.tsv"
    meta_path = tmp_path / "meta.tsv"

    write_signatures(sig_path, {'S1': list(range(10))})
    write_metadata(meta_path, [('SOMETHING_ELSE', 'imatinib')])

    with pytest.raises(ValueError, match="matched"):
        parse_lincs(sig_path, meta_path)


def test_duplicate_metadata_rows_do_not_break_the_join(tmp_path):
    """A repeated signature ID would otherwise make the reindex raise."""
    sig_path = tmp_path / "sigs.tsv"
    meta_path = tmp_path / "meta.tsv"

    write_signatures(sig_path, {'S1': list(range(10))})
    write_metadata(meta_path, [('S1', 'imatinib'), ('S1', 'imatinib')])

    assert set(parse_lincs(sig_path, meta_path, top_pct=0.2)) == {'imatinib'}


# --------------------------------------------------------------------------
# Module size
# --------------------------------------------------------------------------

def test_top_pct_selects_the_right_number_of_genes(lincs):
    signatures = parse_lincs(*lincs, top_pct=0.3)

    assert len(signatures['dasatinib']['up']) == 3
    assert len(signatures['dasatinib']['down']) == 3


def test_a_small_matrix_still_yields_a_usable_module(tmp_path):
    """
    int(10 * 0.05) is 0, so the modules came back empty. A module of no genes
    is never the intended answer.
    """
    sig_path = tmp_path / "sigs.tsv"
    meta_path = tmp_path / "meta.tsv"

    write_signatures(sig_path, {'S1': list(range(10))})
    write_metadata(meta_path, [('S1', 'imatinib')])

    signatures = parse_lincs(sig_path, meta_path, top_pct=0.05)

    assert signatures['imatinib']['up'] == ['G9']
    assert signatures['imatinib']['down'] == ['G0']


@pytest.mark.parametrize("top_pct", [0.0, -0.1, 1.5])
def test_an_out_of_range_top_pct_is_rejected(lincs, top_pct):
    with pytest.raises(ValueError, match="top_pct"):
        parse_lincs(*lincs, top_pct=top_pct)
