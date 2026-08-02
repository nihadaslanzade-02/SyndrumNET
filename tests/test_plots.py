"""
Tests for the plotting layer.

A chart is hard to assert on directly, so these tests check the three things
that can be wrong without anyone noticing: the file is written, no figure is
left open, and the artists on the axes encode the numbers that went in. The
last one matters most for `plot_top_predictions`, whose bar geometry used to
disagree with the totals it was drawing.
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from syndrumnet.viz import plots  # noqa: E402
from syndrumnet.viz.plots import (  # noqa: E402
    plot_auc_comparison,
    plot_degree_distribution,
    plot_pr_curve,
    plot_roc_curve,
    plot_score_distributions,
    plot_top_predictions,
)


@pytest.fixture(autouse=True)
def no_leaked_figures():
    """
    Every test starts and ends with no open figures.

    This is the assertion the module needed most: a save that raised used to
    leave its figure open, and nothing in a normal run would ever say so.
    """
    plt.close("all")
    yield
    open_figures = plt.get_fignums()
    plt.close("all")
    assert open_figures == [], f"figures left open: {open_figures}"


@pytest.fixture
def axes(monkeypatch):
    """
    Keep the figure alive past the save so its artists can be inspected.

    Returns a one-element list that the plot call fills in. Teardown closes
    the figure, which runs before the leak check above.
    """
    captured = []

    def capture(fig, output_path, description):
        from pathlib import Path

        path = Path(output_path)
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path)
        captured.append(fig)

    monkeypatch.setattr(plots, "_save", capture)

    yield captured

    for fig in captured:
        plt.close(fig)


@pytest.fixture
def predictions() -> pd.DataFrame:
    """
    Predictions with the sign pattern real output has.

    PQAB is an averaged proximity z-score, so it is negative for any pair
    close to the disease. CQAB takes either sign. Only TQAB is non-negative.
    Totals are kept distinct so rank order is unambiguous.
    """
    frame = pd.DataFrame({
        'drug_a': ['a', 'c', 'e'],
        'drug_b': ['b', 'd', 'f'],
        'tqab': [2.0, 2.0, 0.0],
        'pqab': [-1.4, -0.6, -2.0],
        'cqab': [0.5, -0.5, 0.4],
    })
    frame['prediction_score'] = frame['tqab'] + frame['pqab'] + frame['cqab']
    return frame


def bars_by_row(ax):
    """
    Bar geometry grouped by row, as a list of (x, signed_width) per y centre.

    matplotlib stores a negative-width bar as a rectangle whose left edge is
    the baseline and whose width is negative, so the sign is recoverable and
    the occupied span is the min/max of the two ends.
    """
    rows = {}

    for container in ax.containers:
        for patch in getattr(container, 'patches', []):
            centre = round(patch.get_y() + patch.get_height() / 2, 6)
            rows.setdefault(centre, []).append((patch.get_x(), patch.get_width()))

    return [rows[centre] for centre in sorted(rows)]


# --------------------------------------------------------------------------
# Saving
# --------------------------------------------------------------------------

def test_save_creates_missing_directories(tmp_path):
    output = tmp_path / "nested" / "deeper" / "degree.png"

    plot_degree_distribution(nx.barabasi_albert_graph(50, 2, seed=0), output)

    assert output.exists()
    assert output.stat().st_size > 0


def test_figure_is_released_when_saving_fails(tmp_path):
    """
    A directory sitting where the file should go makes savefig raise. The
    figure has to be closed anyway, which the autouse fixture then checks.
    """
    blocked = tmp_path / "degree.png"
    blocked.mkdir()

    with pytest.raises(OSError):
        plot_degree_distribution(nx.path_graph(10), blocked)


def test_roc_and_pr_curves_write_files(tmp_path):
    fpr = np.linspace(0, 1, 20)
    tpr = np.sqrt(fpr)

    plot_roc_curve(fpr, tpr, 0.81, tmp_path / "roc.png")
    plot_pr_curve(tpr, fpr, 0.64, tmp_path / "pr.png")

    assert (tmp_path / "roc.png").exists()
    assert (tmp_path / "pr.png").exists()


# --------------------------------------------------------------------------
# Degree distribution
# --------------------------------------------------------------------------

def test_degree_distribution_uses_log_spaced_bins(tmp_path, axes):
    """
    Linear bins on a log axis crowd every low degree into the left edge, so
    the bin edges have to grow geometrically rather than by a fixed step.
    """
    plot_degree_distribution(
        nx.barabasi_albert_graph(200, 2, seed=1), tmp_path / "d.png", log_scale=True
    )

    ax = axes[0].axes[0]
    edges = np.array(sorted({patch.get_x() for patch in ax.patches}))
    ratios = edges[1:] / edges[:-1]

    assert ax.get_xscale() == 'log'
    assert edges[0] == pytest.approx(1.0)
    assert np.allclose(ratios, ratios[0]), "bin edges are not geometrically spaced"


def test_degree_distribution_reports_isolated_nodes(tmp_path, axes):
    """
    Isolated nodes have no place on a logarithmic axis. Dropping them is the
    only option, so the count has to appear in the title.
    """
    G = nx.barabasi_albert_graph(100, 2, seed=2)
    G.add_nodes_from(f"iso{i}" for i in range(7))

    plot_degree_distribution(G, tmp_path / "d.png", log_scale=True)

    assert "7 isolated nodes" in axes[0].axes[0].get_title()


def test_degree_distribution_linear_scale_keeps_isolated_nodes(tmp_path):
    G = nx.path_graph(10)
    G.add_nodes_from(["iso1", "iso2"])

    plot_degree_distribution(G, tmp_path / "d.png", log_scale=False)

    assert (tmp_path / "d.png").exists()


def test_degree_distribution_rejects_an_empty_graph(tmp_path):
    """Used to emit a matplotlib warning and write a blank axis."""
    with pytest.raises(ValueError, match="no nodes"):
        plot_degree_distribution(nx.Graph(), tmp_path / "d.png")


def test_degree_distribution_rejects_an_all_isolated_graph(tmp_path):
    G = nx.Graph()
    G.add_nodes_from(range(5))

    with pytest.raises(ValueError, match="isolated"):
        plot_degree_distribution(G, tmp_path / "d.png", log_scale=True)


# --------------------------------------------------------------------------
# Top predictions
# --------------------------------------------------------------------------

def test_stacked_segments_do_not_overlap(tmp_path, axes, predictions):
    """
    The failure this pins: stacking every component from one running offset
    drew the negative PQAB segment back across the positive TQAB segment, so
    two components claimed the same span and the bar stopped reading as a
    decomposition. Segments must partition their side of zero.
    """
    plot_top_predictions(predictions, k=3, output_path=tmp_path / "top.png")

    for row in bars_by_row(axes[0].axes[0]):
        spans = [(min(x, x + w), max(x, x + w)) for x, w in row]

        for i, (left_a, right_a) in enumerate(spans):
            for left_b, right_b in spans[i + 1:]:
                overlap = min(right_a, right_b) - max(left_a, left_b)
                assert overlap <= 1e-9, (
                    f"segments [{left_a}, {right_a}] and [{left_b}, {right_b}] overlap"
                )


def test_segment_widths_sum_to_the_total(tmp_path, axes, predictions):
    """The signed widths of a row's segments must add up to its total score."""
    plot_top_predictions(predictions, k=3, output_path=tmp_path / "top.png")

    totals = [sum(width for _, width in row) for row in bars_by_row(axes[0].axes[0])]
    expected = sorted(predictions['prediction_score'], reverse=True)

    assert sorted(totals, reverse=True) == pytest.approx(expected)


def test_best_prediction_is_at_the_top(tmp_path, axes, predictions):
    """
    barh puts index 0 at the bottom, so an uninverted axis buries the
    highest-scoring pair at the foot of the chart.
    """
    plot_top_predictions(predictions, k=3, output_path=tmp_path / "top.png")

    ax = axes[0].axes[0]
    bottom, top = ax.get_ylim()
    labels = [text.get_text() for text in ax.get_yticklabels()]

    assert top < bottom, "y axis is not inverted, the best pair renders at the bottom"
    assert labels == ["a-b", "c-d", "e-f"]


def test_total_marker_matches_the_prediction_score(tmp_path, axes, predictions):
    """The total is drawn, not inferred from where the bars happen to end."""
    plot_top_predictions(predictions, k=3, output_path=tmp_path / "top.png")

    ax = axes[0].axes[0]
    marked = sorted(ax.collections[0].get_offsets().data[:, 0], reverse=True)

    assert marked == pytest.approx(sorted(predictions['prediction_score'], reverse=True))


@pytest.mark.parametrize(
    "dropped", ['tqab', 'pqab', 'cqab', 'prediction_score', 'drug_a']
)
def test_missing_columns_are_named(tmp_path, predictions, dropped):
    with pytest.raises(KeyError, match=dropped):
        plot_top_predictions(
            predictions.drop(columns=[dropped]), k=2, output_path=tmp_path / "x.png"
        )


def test_k_larger_than_the_frame_is_not_an_error(tmp_path, predictions):
    plot_top_predictions(predictions, k=50, output_path=tmp_path / "top.png")

    assert (tmp_path / "top.png").exists()


# --------------------------------------------------------------------------
# Score distributions and AUC comparison
# --------------------------------------------------------------------------

def test_score_distributions_writes_a_file(tmp_path, predictions):
    plot_score_distributions(predictions, tmp_path / "dist.png")

    assert (tmp_path / "dist.png").exists()


def test_score_distributions_tolerates_an_all_nan_component(tmp_path, predictions):
    predictions['cqab'] = np.nan

    plot_score_distributions(predictions, tmp_path / "dist.png")

    assert (tmp_path / "dist.png").exists()


def test_score_distributions_rejects_a_frame_with_no_scores(tmp_path):
    frame = pd.DataFrame({'drug_a': ['a'], 'drug_b': ['b']})

    with pytest.raises(KeyError, match="score columns"):
        plot_score_distributions(frame, tmp_path / "dist.png")


def test_auc_comparison_writes_a_file(tmp_path):
    results = {
        'asthma': {'auc_roc': 0.71, 'auc_pr': 0.44},
        'diabetes_t2': {'auc_roc': 0.66, 'auc_pr': 0.39},
    }

    plot_auc_comparison(results, tmp_path / "auc.png")

    assert (tmp_path / "auc.png").exists()


def test_auc_comparison_names_the_disease_missing_a_metric(tmp_path):
    """A bare KeyError('auc_pr') does not say which disease was short."""
    results = {
        'asthma': {'auc_roc': 0.71, 'auc_pr': 0.44},
        'diabetes_t2': {'auc_roc': 0.66},
    }

    with pytest.raises(KeyError, match="diabetes_t2"):
        plot_auc_comparison(results, tmp_path / "auc.png")


def test_auc_comparison_rejects_empty_results(tmp_path):
    with pytest.raises(ValueError, match="empty"):
        plot_auc_comparison({}, tmp_path / "auc.png")
