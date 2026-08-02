"""
Shared fixtures: a small synthetic network with known topology.

The pipeline's real inputs are tens of gigabytes of downloaded data, so the
integration tests run against a hand-built graph whose structure is chosen so
the expected answers are obvious by inspection rather than by running the code
and recording whatever it printed.

Layout, a chain of four communities of six genes each:

    FAR --- RIGHT --- DISEASE --- LEFT

Every community is internally complete; adjacent communities are joined by a
single bridge edge. So LEFT and RIGHT are both adjacent to DISEASE but sit in
different neighbourhoods, while FAR is two communities away from it.
"""

import networkx as nx
import pytest

COMMUNITIES = ["FAR", "RIGHT", "DISEASE", "LEFT"]
COMMUNITY_SIZE = 6


def _genes(prefix: str) -> list[str]:
    return [f"{prefix}{i}" for i in range(COMMUNITY_SIZE)]


@pytest.fixture
def network() -> nx.Graph:
    """Four complete communities joined in a chain by single bridge edges."""
    G = nx.Graph()

    for prefix in COMMUNITIES:
        genes = _genes(prefix)
        for i, left in enumerate(genes):
            for right in genes[i + 1:]:
                G.add_edge(left, right)

    for first, second in zip(COMMUNITIES, COMMUNITIES[1:]):
        G.add_edge(f"{first}0", f"{second}0")

    return G


@pytest.fixture
def disease_module() -> set[str]:
    """Disease gene module, the four genes of DISEASE furthest from bridges."""
    return {f"DISEASE{i}" for i in range(2, 6)}


@pytest.fixture
def drug_modules() -> dict[str, dict[str, set[str]]]:
    """
    Four drugs with deliberately different relationships to the disease.

    left and right sit in different communities, both adjacent to the disease,
    so they cover complementary territory. left_overlap shares two of its three
    genes with left, so the pair is redundant. far sits two communities away.
    """
    return {
        "left": {"up": {"LEFT1", "LEFT2"}, "down": {"LEFT3"}},
        "left_overlap": {"up": {"LEFT1", "LEFT2"}, "down": {"LEFT4"}},
        "right": {"up": {"RIGHT1", "RIGHT2"}, "down": {"RIGHT3"}},
        "far": {"up": {"FAR1", "FAR2"}, "down": {"FAR3"}},
    }


@pytest.fixture
def disease_signature() -> dict[str, float]:
    """
    Disease fold changes.

    LEFT genes move opposite to how drug `left` moves them, so `left` reverses
    the disease signature and should score well on the transcriptional axis.
    Values are distinct so the rank correlation is well defined.
    """
    return {
        "DISEASE2": 2.4, "DISEASE3": 1.9, "DISEASE4": -1.7, "DISEASE5": -2.2,
        "LEFT1": -1.5, "LEFT2": -2.1, "LEFT3": 1.8, "LEFT4": 0.9,
        "RIGHT1": 1.1, "RIGHT2": -0.7, "RIGHT3": 0.4,
        "FAR1": 0.3, "FAR2": 0.6, "FAR3": -0.2,
    }
