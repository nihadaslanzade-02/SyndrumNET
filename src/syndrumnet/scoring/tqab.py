"""
Topological class score (TQAB).

Classifies where a drug pair sits relative to a disease module in the
interactome, following the six-class scheme of Cheng et al. (2019) as adopted
by Iida et al. (2024).
"""

import logging
from typing import Dict, List, Optional, Set, Tuple

import networkx as nx

from syndrumnet.metrics.distances import separation_score

logger = logging.getLogger(__name__)

#: Score awarded to Complementary Exposure. Every other class scores zero.
COMPLEMENTARY_EXPOSURE_SCORE = 2.0


class TopologyClass:
    """
    The six drug-drug-disease classes of Cheng et al. (2019), Figure 2.

    They partition cleanly on two axes: how many of the two drugs sit closer
    to the disease than chance (both, one, neither), and whether the two drug
    modules occupy the same network neighbourhood or separate ones.

    Only Complementary Exposure is associated with therapeutic synergy: both
    drugs reach the disease, but through different neighbourhoods, so their
    effects add rather than duplicate.
    """

    OVERLAPPING_EXPOSURE = "overlapping_exposure"        # P1
    COMPLEMENTARY_EXPOSURE = "complementary_exposure"    # P2, the synergistic one
    INDIRECT_EXPOSURE = "indirect_exposure"              # P3
    SINGLE_EXPOSURE = "single_exposure"                  # P4
    NON_EXPOSURE = "non_exposure"                        # P5
    INDEPENDENT_ACTION = "independent_action"            # P6


def classify_topology(z_qa: float, z_qb: float, s_ab: float) -> str:
    """
    Assign a drug pair to one of the six topological classes.

    Parameters
    ----------
    z_qa : float
        Proximity z-score of drug A to the disease. Negative means closer
        than chance, which the source calls "overlapping with the disease".
    z_qb : float
        Proximity z-score of drug B to the disease.
    s_ab : float
        Network separation between the two drug modules. Negative means the
        two drug modules share a neighbourhood; non-negative means they are
        topologically separated.

    Returns
    -------
    str
        One of the TopologyClass constants.

    Notes
    -----
    From Cheng et al. (2019), Figure 2:

        P1 Overlapping Exposure     z_QA < 0, z_QB < 0, s_AB <  0
        P2 Complementary Exposure   z_QA < 0, z_QB < 0, s_AB >= 0
        P3 Indirect Exposure        one drug < 0,       s_AB <  0
        P4 Single Exposure          one drug < 0,       s_AB >= 0
        P5 Non-exposure             neither  < 0,       s_AB <  0
        P6 Independent Action       neither  < 0,       s_AB >= 0

    The sign convention on s_AB is inherited from Menche et al. (2015) via
    Cheng et al.: "For sAB < 0, the targets of the two drugs are located in
    the same network neighborhood, while for sAB >= 0, the two drug targets
    are topologically separated."

    Note that Iida et al. (2024) states the Class II condition inline as
    s_AB < 0, which contradicts both its own description of Class II as "two
    separated drug modules" and the source it cites. The condition used here
    is s_AB >= 0, matching Cheng et al.'s Figure 2 panel P2 and the prose in
    both papers.
    """
    a_near, b_near = z_qa < 0, z_qb < 0
    separated = s_ab >= 0

    if a_near and b_near:
        return (
            TopologyClass.COMPLEMENTARY_EXPOSURE
            if separated
            else TopologyClass.OVERLAPPING_EXPOSURE
        )

    if a_near or b_near:
        return (
            TopologyClass.SINGLE_EXPOSURE
            if separated
            else TopologyClass.INDIRECT_EXPOSURE
        )

    return (
        TopologyClass.INDEPENDENT_ACTION
        if separated
        else TopologyClass.NON_EXPOSURE
    )


def compute_tqab(z_qa: float, z_qb: float, s_ab: float) -> Tuple[float, str]:
    """
    Compute the topological class score T_QAB.

    Parameters
    ----------
    z_qa : float
        Proximity z-score of drug A to the disease module.
    z_qb : float
        Proximity z-score of drug B to the disease module.
    s_ab : float
        Network separation between the two drug modules.

    Returns
    -------
    tuple
        (tqab_score, topology_class)

    Notes
    -----
    The score is binary, not graded: T_QAB = 2 for Complementary Exposure and
    0 for every other class, per Iida et al. (2024). Since the final
    prediction is the unweighted sum T_QAB + P_QAB + C_QAB, that constant of 2
    is what sets the topological axis' weight against the other two.
    """
    topology_class = classify_topology(z_qa, z_qb, s_ab)

    score = (
        COMPLEMENTARY_EXPOSURE_SCORE
        if topology_class == TopologyClass.COMPLEMENTARY_EXPOSURE
        else 0.0
    )

    logger.debug(
        f"TQAB: z_QA={z_qa:.3f}, z_QB={z_qb:.3f}, s_AB={s_ab:.3f}, "
        f"class={topology_class}, score={score}"
    )

    return score, topology_class


def compute_tqab_batch(
    G: nx.Graph,
    disease_module: Set[str],
    drug_modules: Dict[str, Set[str]],
    drug_pairs: List[Tuple[str, str]],
    proximity_zscores: Optional[Dict[str, float]] = None,
    n_randomizations: int = 1000,
    seed: int = 42,
) -> Dict[Tuple[str, str], Tuple[float, str]]:
    """
    Compute TQAB for multiple drug pairs.

    Parameters
    ----------
    G : nx.Graph
        Network graph.
    disease_module : set
        Disease gene module.
    drug_modules : dict
        {drug_name: gene_module}
    drug_pairs : list of tuple
        Drug pair identifiers [(drug_a, drug_b), ...].
    proximity_zscores : dict, optional
        Precomputed {drug_name: z_score}. Supply the values already computed
        for PQAB rather than paying for the null model twice.
    n_randomizations : int
        Randomizations per z-score, when they have to be computed here.
    seed : int
        Run-level seed, when z-scores have to be computed here.

    Returns
    -------
    dict
        {(drug_a, drug_b): (tqab_score, topology_class)}

    Notes
    -----
    z-scores are per drug and are cached; separation is genuinely per pair and
    is not.
    """
    from syndrumnet.scoring.pqab import proximity_zscore

    needed = {drug for pair in drug_pairs for drug in pair if drug in drug_modules}

    if proximity_zscores is None:
        proximity_zscores = {
            drug: proximity_zscore(
                G, disease_module, drug_modules[drug], n_randomizations, seed
            )
            for drug in sorted(needed)
        }

    results = {}

    for drug_a, drug_b in drug_pairs:
        if drug_a not in drug_modules or drug_b not in drug_modules:
            logger.warning(f"Missing module for pair ({drug_a}, {drug_b})")
            continue
        if drug_a not in proximity_zscores or drug_b not in proximity_zscores:
            logger.warning(f"Missing z-score for pair ({drug_a}, {drug_b})")
            continue

        s_ab = separation_score(G, drug_modules[drug_a], drug_modules[drug_b])

        results[(drug_a, drug_b)] = compute_tqab(
            proximity_zscores[drug_a], proximity_zscores[drug_b], s_ab
        )

    return results
