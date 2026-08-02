"""
Tests for the similarity layers.

The Tanimoto tests carry most of the weight. A binary fingerprint matrix has
several equally natural dtypes, and one of them used to give a silently wrong
answer, so the ground truth here is computed independently with set algebra
rather than by recording what the implementation returns.
"""

import numpy as np
import pytest

from syndrumnet.propagation.similarity_layers import (
    build_similarity_matrix,
    compute_disease_similarity,
    compute_drug_similarity,
    jaccard_similarity,
    overlap_coefficient,
    tanimoto_similarity_matrix,
)

BITS = np.array([
    [1, 1, 0, 0],
    [1, 0, 1, 0],
    [0, 0, 1, 1],
])


def reference_tanimoto(a: np.ndarray, b: np.ndarray) -> float:
    """Tanimoto from set algebra, independent of the implementation."""
    set_a = {i for i, bit in enumerate(a) if bit}
    set_b = {i for i, bit in enumerate(b) if bit}

    if not set_a and not set_b:
        return 1.0

    return len(set_a & set_b) / len(set_a | set_b)


@pytest.fixture
def reference_matrix() -> np.ndarray:
    n = len(BITS)
    return np.array([
        [reference_tanimoto(BITS[i], BITS[j]) for j in range(n)]
        for i in range(n)
    ])


# --------------------------------------------------------------------------
# Tanimoto
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "dtype", [bool, np.int8, np.uint8, np.int32, np.int64, np.float32, np.float64]
)
def test_tanimoto_is_dtype_independent(dtype, reference_matrix):
    """
    Every reasonable encoding of a binary matrix gives the same answer.

    bool is the case that mattered: numpy's ``@`` on booleans is logical, so
    the shared-bit count saturated at 1 and the matrix came out wrong without
    raising. The diagonal was 0.33 instead of 1.0.
    """
    result = tanimoto_similarity_matrix(BITS.astype(dtype))

    assert np.allclose(result, reference_matrix)


def test_tanimoto_diagonal_is_one():
    """A molecule is identical to itself whatever its bit count."""
    fingerprints = np.array([[1, 1, 0], [0, 0, 0], [1, 1, 1]])

    result = tanimoto_similarity_matrix(fingerprints)

    assert np.allclose(np.diag(result), 1.0)


def test_tanimoto_is_symmetric_and_bounded(reference_matrix):
    result = tanimoto_similarity_matrix(BITS)

    assert np.allclose(result, result.T)
    assert result.min() >= 0.0
    assert result.max() <= 1.0


def test_tanimoto_rejects_non_binary_input():
    """
    Counted fingerprints need a different denominator, so accepting them
    would return a plausible-looking wrong number.
    """
    with pytest.raises(ValueError, match="binary"):
        tanimoto_similarity_matrix(np.array([[1, 2], [0, 1]]))


def test_tanimoto_rejects_a_single_fingerprint():
    with pytest.raises(ValueError, match="2-D"):
        tanimoto_similarity_matrix(np.array([1, 0, 1, 0]))


def test_tanimoto_matches_hand_computed_value():
    """Molecules 0 and 1 share bit 0 and cover bits 0, 1, 2 between them."""
    result = tanimoto_similarity_matrix(BITS)

    assert result[0, 1] == pytest.approx(1 / 3)
    assert result[0, 2] == pytest.approx(0.0)


# --------------------------------------------------------------------------
# Set metrics
# --------------------------------------------------------------------------

def test_jaccard_known_values():
    assert jaccard_similarity({"A", "B"}, {"B", "C"}) == pytest.approx(1 / 3)
    assert jaccard_similarity({"A"}, {"A"}) == 1.0
    assert jaccard_similarity({"A"}, {"B"}) == 0.0


def test_overlap_known_values():
    assert overlap_coefficient({"A"}, {"A", "B"}) == 1.0
    assert overlap_coefficient({"A", "B"}, {"B", "C"}) == pytest.approx(0.5)


@pytest.mark.parametrize("metric", [jaccard_similarity, overlap_coefficient])
def test_empty_set_convention_is_shared(metric):
    """
    Both metrics use the same convention, which they did not always: two
    empty modules are identical, one empty module shares nothing.
    """
    assert metric(set(), set()) == 1.0
    assert metric(set(), {"A"}) == 0.0
    assert metric({"A"}, set()) == 0.0


# --------------------------------------------------------------------------
# Matrix builders
# --------------------------------------------------------------------------

def test_disease_similarity_follows_key_order():
    """
    The matrix has no labels, so its only contract is that rows follow the
    input key order. A caller that re-derives the order must get the same one.
    """
    modules = {
        "z_disease": {"A", "B"},
        "a_disease": {"B", "C"},
        "m_disease": {"X"},
    }

    result = compute_disease_similarity(modules)

    assert result.shape == (3, 3)
    assert np.allclose(np.diag(result), 1.0)
    assert result[0, 1] == pytest.approx(1 / 3)   # z_disease against a_disease
    assert result[0, 2] == 0.0                    # z_disease against m_disease


def test_disease_similarity_is_symmetric():
    modules = {f"d{i}": {f"g{i}", f"g{i + 1}"} for i in range(5)}

    result = compute_disease_similarity(modules)

    assert np.allclose(result, result.T)


def test_disease_similarity_overlap_method_differs_from_jaccard():
    modules = {"a": {"X"}, "b": {"X", "Y", "Z"}}

    jaccard = compute_disease_similarity(modules, method='jaccard')
    overlap = compute_disease_similarity(modules, method='overlap')

    assert jaccard[0, 1] == pytest.approx(1 / 3)
    assert overlap[0, 1] == pytest.approx(1.0)


@pytest.mark.parametrize(
    "modules",
    [
        pytest.param({}, id="no_diseases"),
        pytest.param({"only": {"G1"}}, id="one_disease"),
        pytest.param({"a": {"G1"}, "b": {"G2"}}, id="two_diseases"),
    ],
)
def test_unknown_method_always_raises(modules):
    """
    The check has to happen before the pair loop. It used to sit inside it,
    so a typo passed silently whenever there was no off-diagonal pair to
    evaluate, and a caller could get a matrix back from a metric that does
    not exist.
    """
    with pytest.raises(ValueError, match="Unknown method"):
        compute_disease_similarity(modules, method='not_a_metric')


def test_drug_similarity_uses_tanimoto(reference_matrix):
    fingerprints = {name: BITS[i] for i, name in enumerate("abc")}

    result = compute_drug_similarity(fingerprints)

    assert np.allclose(result, reference_matrix)


def test_drug_similarity_handles_boolean_fingerprints(reference_matrix):
    """RDKit bit vectors convert naturally to bool arrays."""
    fingerprints = {name: BITS[i].astype(bool) for i, name in enumerate("abc")}

    assert np.allclose(compute_drug_similarity(fingerprints), reference_matrix)


def test_drug_similarity_rejects_empty_input():
    """Used to surface as 'need at least one array to concatenate'."""
    with pytest.raises(ValueError, match="empty"):
        compute_drug_similarity({})


def test_drug_similarity_rejects_mixed_fingerprint_lengths():
    fingerprints = {"a": np.array([1, 0, 1]), "b": np.array([1, 0])}

    with pytest.raises(ValueError, match="mixed lengths"):
        compute_drug_similarity(fingerprints)


def test_drug_similarity_rejects_unknown_method():
    with pytest.raises(ValueError, match="Unknown method"):
        compute_drug_similarity({"a": BITS[0]}, method='not_a_metric')


def test_build_similarity_matrix_calls_each_pair_once():
    """Symmetry is assumed, so the callback must not be invoked twice."""
    calls = []

    def counting(a, b):
        calls.append((a, b))
        return 0.5

    result = build_similarity_matrix(["x", "y", "z"], counting)

    assert calls == [("x", "y"), ("x", "z"), ("y", "z")]
    assert np.allclose(np.diag(result), 1.0)
    assert result[1, 0] == 0.5


def test_build_similarity_matrix_passes_kwargs():
    def scaled(a, b, factor=1.0):
        return factor

    result = build_similarity_matrix(["x", "y"], scaled, factor=0.25)

    assert result[0, 1] == 0.25
