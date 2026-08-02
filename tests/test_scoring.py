"""Tests for scoring functions."""

import pytest

from syndrumnet.scoring.cqab import compute_cqab
from syndrumnet.scoring.tqab import (
    COMPLEMENTARY_EXPOSURE_SCORE,
    TopologyClass,
    classify_topology,
    compute_tqab,
)

# The six classes of Cheng et al. (2019), Figure 2, as (z_QA, z_QB, s_AB)
# sign patterns. z < 0 means the drug is closer to the disease than chance;
# s_AB >= 0 means the two drug modules are topologically separated.
CLASS_TABLE = [
    (-1.0, -1.0, -1.0, TopologyClass.OVERLAPPING_EXPOSURE),    # P1
    (-1.0, -1.0, +1.0, TopologyClass.COMPLEMENTARY_EXPOSURE),  # P2
    (-1.0, +1.0, -1.0, TopologyClass.INDIRECT_EXPOSURE),       # P3
    (-1.0, +1.0, +1.0, TopologyClass.SINGLE_EXPOSURE),         # P4
    (+1.0, +1.0, -1.0, TopologyClass.NON_EXPOSURE),            # P5
    (+1.0, +1.0, +1.0, TopologyClass.INDEPENDENT_ACTION),      # P6
]


@pytest.mark.parametrize("z_qa,z_qb,s_ab,expected", CLASS_TABLE)
def test_classify_topology_matches_the_published_table(z_qa, z_qb, s_ab, expected):
    """Every cell of the six-class table is reproduced exactly."""
    assert classify_topology(z_qa, z_qb, s_ab) == expected


@pytest.mark.parametrize("z_qa,z_qb,s_ab,expected", CLASS_TABLE)
def test_classification_is_symmetric_in_the_drug_pair(z_qa, z_qb, s_ab, expected):
    """A pair is unordered, so swapping the two drugs cannot change its class."""
    assert classify_topology(z_qb, z_qa, s_ab) == expected


def test_only_complementary_exposure_scores():
    """
    T_QAB is binary: 2 for Complementary Exposure, 0 for everything else.

    The magnitude matters because the final prediction is an unweighted sum of
    the three components, so this constant sets the topological axis' weight.
    """
    for z_qa, z_qb, s_ab, expected_class in CLASS_TABLE:
        score, topology_class = compute_tqab(z_qa, z_qb, s_ab)

        assert topology_class == expected_class
        if expected_class == TopologyClass.COMPLEMENTARY_EXPOSURE:
            assert score == COMPLEMENTARY_EXPOSURE_SCORE
        else:
            assert score == 0.0


def test_separated_boundary_counts_as_separated():
    """
    The published condition is s_AB >= 0, not s_AB > 0.

    Exactly zero separation is the boundary case, and it belongs on the
    separated side.
    """
    _, topology_class = compute_tqab(-1.0, -1.0, 0.0)

    assert topology_class == TopologyClass.COMPLEMENTARY_EXPOSURE


def test_cqab_inverse_correlation():
    """Test CQAB computes inverse correlation."""
    disease_sig = {'G1': 2.0, 'G2': 1.5, 'G3': -1.0}

    # Drug A reverses disease (up genes are down in disease)
    drug_a_up = {'G3'}  # Was down in disease
    drug_a_down = {'G1', 'G2'}  # Were up in disease

    # Drug B same as disease (redundant)
    drug_b_up = {'G1', 'G2'}
    drug_b_down = {'G3'}

    cqab, cqa, cqb = compute_cqab(
        disease_sig,
        drug_a_up, drug_a_down,
        drug_b_up, drug_b_down,
    )

    # Drug A should have positive correlation (reverses disease)
    # Drug B should have negative correlation (same as disease)
    assert cqa > cqb
