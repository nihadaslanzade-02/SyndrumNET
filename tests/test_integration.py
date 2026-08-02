"""
End-to-end tests for the prediction pipeline.

The unit tests cover each scoring component in isolation. These run the whole
SynergyPredictor over the synthetic network from conftest.py, which is where
wiring faults show up: a component can be correct on its own and still be
handed the wrong argument, seeded inconsistently, or recomputed per pair when
it only depends on one drug.
"""

import networkx as nx
import pytest

from syndrumnet.metrics.distances import separation_score, shortest_path_distance
from syndrumnet.scoring.pqab import compute_pqab, compute_pqab_batch, proximity_zscore
from syndrumnet.scoring.predictor import SynergyPredictor
from syndrumnet.scoring.tqab import (
    COMPLEMENTARY_EXPOSURE_SCORE,
    TopologyClass,
    compute_tqab_batch,
)

N_RANDOMIZATIONS = 50
SEED = 42


@pytest.fixture
def predictor(network, disease_module, drug_modules, disease_signature):
    """A predictor loaded with the synthetic disease and drugs."""
    predictor = SynergyPredictor(
        network, n_randomizations=N_RANDOMIZATIONS, seed=SEED
    )
    predictor.set_disease_modules({"synth": disease_module})
    predictor.set_drug_modules(drug_modules)
    predictor.set_disease_signatures({"synth": disease_signature})
    return predictor


@pytest.fixture
def predictions(predictor):
    """Predictions for every drug pair, computed once per test."""
    return predictor.predict_all("synth")


# ---------------------------------------------------------------------------
# Shape and bookkeeping
# ---------------------------------------------------------------------------


def test_predicts_every_drug_pair(predictions, drug_modules):
    """Every unordered pair appears exactly once."""
    n_drugs = len(drug_modules)
    assert len(predictions) == n_drugs * (n_drugs - 1) // 2

    pairs = {frozenset(pair) for pair in zip(predictions.drug_a, predictions.drug_b)}
    assert len(pairs) == len(predictions)
    assert all(len(pair) == 2 for pair in pairs)


def test_no_missing_values(predictions):
    """A silent NaN would propagate into the ranking unnoticed."""
    assert not predictions.isna().to_numpy().any()


def test_score_is_the_sum_of_its_components(predictions):
    """prediction_score must equal TQAB + PQAB + CQAB, not drift from it."""
    total = predictions.tqab + predictions.pqab + predictions.cqab
    assert (total - predictions.prediction_score).abs().max() < 1e-12


def test_results_are_ranked(predictions):
    """Output is sorted best first, so head() is the shortlist."""
    scores = predictions.prediction_score.tolist()
    assert scores == sorted(scores, reverse=True)


def test_unknown_disease_is_rejected(predictor):
    with pytest.raises(ValueError, match="Unknown disease"):
        predictor.predict_all("not_a_disease")


def test_modules_must_be_set_before_predicting(network):
    predictor = SynergyPredictor(network)
    with pytest.raises(ValueError, match="Modules not set"):
        predictor.predict_all("synth")


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------


def test_same_seed_gives_identical_predictions(predictor):
    """
    The null model is stochastic, so this is the guard on the seeding path.

    A run that cannot be reproduced cannot be reported.
    """
    first = predictor.predict_all("synth")
    second = predictor.predict_all("synth")

    assert first.equals(second)


def test_proximity_does_not_depend_on_pair_position(
    network, disease_module, drug_modules
):
    """
    P_QA is a property of the disease and one drug.

    It must not change with which other drug shares the pair, nor with the
    order the two are passed in. Seeding the null model off the run seed alone
    broke this: the first drug drew with `seed` and the second with `seed + 1`,
    so a drug scored differently depending on which side it landed on.
    """
    left = drug_modules["left"]["up"] | drug_modules["left"]["down"]
    right = drug_modules["right"]["up"] | drug_modules["right"]["down"]

    _, z_left_first, z_right_second = compute_pqab(
        network, disease_module, left, right, N_RANDOMIZATIONS, SEED
    )
    _, z_right_first, z_left_second = compute_pqab(
        network, disease_module, right, left, N_RANDOMIZATIONS, SEED
    )

    assert z_left_first == z_left_second
    assert z_right_first == z_right_second


def test_batch_matches_single_pair(network, disease_module, drug_modules):
    """
    The batch path caches one z-score per drug rather than recomputing per
    pair. Cheaper only counts if it agrees with the direct computation.
    """
    module_sets = {
        name: sig["up"] | sig["down"] for name, sig in drug_modules.items()
    }
    names = sorted(module_sets)
    pairs = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]]

    batch = compute_pqab_batch(
        network, disease_module, module_sets, pairs, N_RANDOMIZATIONS, SEED
    )

    for drug_a, drug_b in pairs:
        direct = compute_pqab(
            network,
            disease_module,
            module_sets[drug_a],
            module_sets[drug_b],
            N_RANDOMIZATIONS,
            SEED,
        )
        assert batch[(drug_a, drug_b)] == pytest.approx(direct, abs=1e-12)


def test_each_drug_has_one_proximity_score(predictions):
    """A drug's z-score is identical in every row it appears in."""
    per_drug = {}
    for row in predictions.itertuples():
        per_drug.setdefault(row.drug_a, set()).add(round(row.pqa, 12))
        per_drug.setdefault(row.drug_b, set()).add(round(row.pqb, 12))

    assert all(len(values) == 1 for values in per_drug.values())


# ---------------------------------------------------------------------------
# The scores mean what they claim to mean
# ---------------------------------------------------------------------------


def test_distant_drug_is_further_from_the_disease(
    network, disease_module, drug_modules
):
    """FAR is two communities away, LEFT is adjacent; the metric must agree."""
    left = drug_modules["left"]["up"] | drug_modules["left"]["down"]
    far = drug_modules["far"]["up"] | drug_modules["far"]["down"]

    assert shortest_path_distance(network, far, disease_module) > (
        shortest_path_distance(network, left, disease_module)
    )
    assert proximity_zscore(
        network, disease_module, far, N_RANDOMIZATIONS, SEED
    ) > proximity_zscore(network, disease_module, left, N_RANDOMIZATIONS, SEED)


def test_overlapping_drugs_are_less_separated(network, drug_modules):
    """
    Drugs sharing genes must separate less than drugs in different
    communities. This is the ordering the topological class depends on.
    """
    left = drug_modules["left"]["up"] | drug_modules["left"]["down"]
    overlap = drug_modules["left_overlap"]["up"] | drug_modules["left_overlap"]["down"]
    right = drug_modules["right"]["up"] | drug_modules["right"]["down"]

    assert separation_score(network, left, right) > (
        separation_score(network, left, overlap)
    )


def test_every_pair_gets_a_valid_class(predictions):
    """The classification is exhaustive, so no pair can fall through it."""
    valid = {
        value
        for name, value in vars(TopologyClass).items()
        if not name.startswith("_")
    }

    assert set(predictions.topology_class) <= valid


def test_tqab_is_binary_and_tracks_the_class(predictions):
    """
    T_QAB is 2 for Complementary Exposure and 0 otherwise, with no
    intermediate values. Earlier revisions produced a graded score from
    unpublished constants instead.
    """
    for row in predictions.itertuples():
        expected = (
            COMPLEMENTARY_EXPOSURE_SCORE
            if row.topology_class == TopologyClass.COMPLEMENTARY_EXPOSURE
            else 0.0
        )
        assert row.tqab == expected


def test_topology_uses_the_shared_zscores(network, disease_module, drug_modules):
    """
    TQAB classifies on the same proximity z-scores PQAB scores on.

    Passing them in must give the same answer as letting the batch compute
    its own, otherwise the two components disagree about how close a drug is.
    """
    module_sets = {
        name: sig["up"] | sig["down"] for name, sig in drug_modules.items()
    }
    names = sorted(module_sets)
    pairs = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]]

    zscores = {
        name: proximity_zscore(
            network, disease_module, module_sets[name], N_RANDOMIZATIONS, SEED
        )
        for name in names
    }

    supplied = compute_tqab_batch(
        network, disease_module, module_sets, pairs, proximity_zscores=zscores
    )
    computed = compute_tqab_batch(
        network,
        disease_module,
        module_sets,
        pairs,
        n_randomizations=N_RANDOMIZATIONS,
        seed=SEED,
    )

    assert supplied == computed


def test_complementary_exposure_reaches_the_score(network, drug_modules):
    """
    The one class that scores runs end to end through the batch path.

    The synthetic network has no drug that is genuinely closer to the disease
    than chance, so the z-scores are supplied directly here; the separation
    term still comes from the real network.
    """
    module_sets = {
        name: sig["up"] | sig["down"] for name, sig in drug_modules.items()
    }
    # left and right sit in different communities, so they are separated.
    both_near = {"left": -1.5, "right": -1.2}

    results = compute_tqab_batch(
        network,
        set(),
        {name: module_sets[name] for name in both_near},
        [("left", "right")],
        proximity_zscores=both_near,
    )

    score, topology_class = results[("left", "right")]

    assert topology_class == TopologyClass.COMPLEMENTARY_EXPOSURE
    assert score == COMPLEMENTARY_EXPOSURE_SCORE


def test_transcriptional_score_rewards_reversal(predictions):
    """
    `left` moves the LEFT genes opposite to the disease, `far` moves them
    with it. The reversing drug must score higher on the CQAB axis.
    """
    per_drug = {}
    for row in predictions.itertuples():
        per_drug[row.drug_a] = row.cqa
        per_drug[row.drug_b] = row.cqb

    assert per_drug["left"] > per_drug["far"]


# ---------------------------------------------------------------------------
# Degenerate inputs
# ---------------------------------------------------------------------------


def test_genes_outside_the_network_are_ignored(
    network, disease_module, drug_modules, disease_signature
):
    """
    Real modules always contain genes the interactome does not cover. Those
    must be dropped, not crash and not shift the score.
    """
    padded = {
        name: {
            "up": sig["up"] | {f"NOT_IN_NETWORK_{name}"},
            "down": sig["down"],
        }
        for name, sig in drug_modules.items()
    }

    def run(modules):
        predictor = SynergyPredictor(
            network, n_randomizations=N_RANDOMIZATIONS, seed=SEED
        )
        predictor.set_disease_modules({"synth": disease_module})
        predictor.set_drug_modules(modules)
        predictor.set_disease_signatures({"synth": disease_signature})
        return predictor.predict_all("synth")

    baseline = run(drug_modules)
    padded_result = run(padded)

    assert padded_result.tqab.tolist() == pytest.approx(baseline.tqab.tolist())


def test_disconnected_network_does_not_produce_infinities(
    disease_module, drug_modules
):
    """
    An unreachable module must yield the large finite sentinel, never inf or
    NaN, so the averages downstream stay meaningful.
    """
    G = nx.Graph()
    G.add_edges_from([("DISEASE2", "DISEASE3"), ("DISEASE4", "DISEASE5")])
    G.add_edges_from([("LEFT1", "LEFT2"), ("LEFT2", "LEFT3")])

    distance = shortest_path_distance(G, {"LEFT1", "LEFT2"}, disease_module)

    assert distance == 1000.0
