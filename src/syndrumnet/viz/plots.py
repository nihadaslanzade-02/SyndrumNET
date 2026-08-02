"""
Visualization functions using matplotlib.

All plots use matplotlib only (no seaborn) as specified.

Every function that takes an `output_path` creates the parent directory,
writes the figure, and releases it even if the write fails, so a failed save
in a loop cannot leak figures.
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Sequence

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

#: Columns `plot_top_predictions` needs to draw a score decomposition.
SCORE_COMPONENTS = ('tqab', 'pqab', 'cqab')
COMPONENT_LABELS = ('TQAB (Topological)', 'PQAB (Proximity)', 'CQAB (Transcriptional)')
REQUIRED_PREDICTION_COLUMNS = ('drug_a', 'drug_b', 'prediction_score') + SCORE_COMPONENTS


def _save(fig: plt.Figure, output_path: Path, description: str) -> None:
    """
    Write a figure to disk and always release it.

    matplotlib figures are owned by the pyplot state machine, not by the
    variable holding them, so an exception between `subplots` and `close`
    leaves the figure open for the rest of the process. In a loop over
    diseases that accumulates until matplotlib starts warning about it.
    """
    output_path = Path(output_path)

    try:
        fig.tight_layout()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
    finally:
        plt.close(fig)

    logger.info(f"Saved {description} to {output_path}")


def _require_columns(df: pd.DataFrame, columns: Sequence[str], what: str) -> None:
    """Fail with the full list of missing columns rather than the first one."""
    missing = [column for column in columns if column not in df.columns]

    if missing:
        raise KeyError(
            f"{what} is missing required column(s) {missing}. "
            f"Present: {list(df.columns)}."
        )


def plot_degree_distribution(
    G: nx.Graph,
    output_path: Path,
    log_scale: bool = True,
) -> None:
    """
    Plot network degree distribution.

    Parameters
    ----------
    G : nx.Graph
        Network graph.
    output_path : Path
        Output file path.
    log_scale : bool
        Use log-log scale.

    Raises
    ------
    ValueError
        If the graph has no nodes, or if `log_scale` is set and every node is
        isolated, leaving nothing that a logarithmic axis can show.

    Notes
    -----
    On a log-log axis the bins are logarithmically spaced. Linear bins crowd
    every low degree into the left edge, and a bin starting at degree 0 drags
    the lower x-limit down towards zero: on a real network with isolated
    nodes that stretched the view across roughly eighteen orders of
    magnitude. Isolated nodes cannot be placed on a logarithmic axis at all,
    so they are counted in the title instead of silently dropped.
    """
    if G.number_of_nodes() == 0:
        raise ValueError("Cannot plot a degree distribution for a graph with no nodes.")

    degrees = np.array([degree for _, degree in G.degree()])
    n_isolated = int((degrees == 0).sum())
    connected = degrees[degrees > 0]

    # Validate before opening a figure, so a rejected input cannot leak one.
    if log_scale and connected.size == 0:
        raise ValueError(
            "Every node is isolated, so a log-scaled degree distribution "
            "would be empty. Pass log_scale=False to plot the zero bar."
        )

    fig, ax = plt.subplots(figsize=(8, 6))

    title = 'Network Degree Distribution'

    if log_scale:
        ax.hist(
            connected,
            bins=np.logspace(0, np.log10(connected.max()), 50),
            edgecolor='black',
            alpha=0.7,
        )
        ax.set_xscale('log')
        ax.set_yscale('log')

        if n_isolated:
            title += f' ({n_isolated} isolated nodes not shown)'
            logger.warning(f"{n_isolated} isolated nodes omitted from the log-scaled plot")
    else:
        ax.hist(degrees, bins=50, edgecolor='black', alpha=0.7)

    ax.set_xlabel('Degree', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    _save(fig, output_path, "degree distribution")


def plot_roc_curve(
    fpr: np.ndarray,
    tpr: np.ndarray,
    auc: float,
    output_path: Path,
    title: str = "ROC Curve",
) -> None:
    """
    Plot ROC curve.

    Parameters
    ----------
    fpr : np.ndarray
        False positive rate.
    tpr : np.ndarray
        True positive rate.
    auc : float
        AUC score.
    output_path : Path
        Output file path.
    title : str
        Plot title.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(fpr, tpr, linewidth=2, label=f'AUC = {auc:.3f}')
    ax.plot([0, 1], [0, 1], 'k--', linewidth=1, label='Random')

    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='lower right', fontsize=11)
    ax.grid(True, alpha=0.3)

    _save(fig, output_path, "ROC curve")


def plot_pr_curve(
    precision: np.ndarray,
    recall: np.ndarray,
    auc_pr: float,
    output_path: Path,
    title: str = "Precision-Recall Curve",
) -> None:
    """
    Plot precision-recall curve.

    Parameters
    ----------
    precision : np.ndarray
        Precision values.
    recall : np.ndarray
        Recall values.
    auc_pr : float
        AUC-PR score.
    output_path : Path
        Output file path.
    title : str
        Plot title.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(recall, precision, linewidth=2, label=f'AUC-PR = {auc_pr:.3f}')

    ax.set_xlabel('Recall', fontsize=12)
    ax.set_ylabel('Precision', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='lower left', fontsize=11)
    ax.grid(True, alpha=0.3)

    _save(fig, output_path, "PR curve")


def plot_score_distributions(
    predictions: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Plot distributions of TQAB, PQAB, CQAB scores.

    Parameters
    ----------
    predictions : pd.DataFrame
        Predictions with score columns.
    output_path : Path
        Output file path.

    Raises
    ------
    KeyError
        If none of the three score columns is present. Individual missing
        columns leave their panel blank, which is informative; all three
        missing means the caller passed the wrong frame.
    """
    if not any(column in predictions.columns for column in SCORE_COMPONENTS):
        raise KeyError(
            f"predictions has none of the score columns {list(SCORE_COMPONENTS)}. "
            f"Present: {list(predictions.columns)}."
        )

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for ax, score, title in zip(axes, SCORE_COMPONENTS, COMPONENT_LABELS):
        ax.set_xlabel('Score', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)

        values = (
            predictions[score].dropna()
            if score in predictions.columns
            else pd.Series(dtype=float)
        )

        if values.empty:
            ax.text(
                0.5, 0.5, 'no data',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=11, color='grey',
            )
            continue

        ax.hist(values, bins=30, edgecolor='black', alpha=0.7)

    _save(fig, output_path, "score distributions")


def plot_top_predictions(
    predictions: pd.DataFrame,
    k: int = 20,
    output_path: Optional[Path] = None,
) -> None:
    """
    Plot top-k predictions as a signed decomposition of the total score.

    Parameters
    ----------
    predictions : pd.DataFrame
        Predictions dataframe.
    k : int
        Number of top predictions to show.
    output_path : Path, optional
        Output file path. If omitted the figure is shown instead of saved,
        and ownership of it passes to the caller.

    Raises
    ------
    KeyError
        If any of the drug, component or total score columns is missing.

    Notes
    -----
    The three components are not all positive, so they cannot be stacked from
    a single running offset. PQAB is an averaged proximity z-score and is
    negative exactly when the two drugs sit closer to the disease than
    chance, which is the case for every pair worth looking at; CQAB is an
    averaged correlation and takes either sign. Stacking naively draws each
    negative segment back over the positive ones, so the segments overlap and
    the bar stops reading as a decomposition even though its end point
    happens to land on the right total.

    Positive and negative parts are therefore stacked from separate
    baselines, and the total is drawn as a separate marker rather than
    inferred from where the bar ends. The y axis is inverted so the
    highest-scoring pair is at the top.
    """
    _require_columns(predictions, REQUIRED_PREDICTION_COLUMNS, "predictions")

    top_k = predictions.nlargest(k, 'prediction_score')
    labels = [f"{a}-{b}" for a, b in zip(top_k['drug_a'], top_k['drug_b'])]

    fig, ax = plt.subplots(figsize=(12, 8))

    y = np.arange(len(labels))
    height = 0.8

    positive_base = np.zeros(len(top_k))
    negative_base = np.zeros(len(top_k))

    for component in SCORE_COMPONENTS:
        values = top_k[component].to_numpy(dtype=float)
        base = np.where(values >= 0, positive_base, negative_base)

        ax.barh(y, values, height, left=base, label=component.upper(), alpha=0.8)

        positive_base += np.clip(values, 0, None)
        negative_base += np.clip(values, None, 0)

    ax.scatter(
        top_k['prediction_score'], y,
        color='black', s=25, zorder=3, label='Total',
    )
    ax.axvline(0, color='black', linewidth=0.8)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('Score', fontsize=12)
    ax.set_title(
        f'Top {k} Predicted Synergistic Combinations',
        fontsize=14, fontweight='bold',
    )
    ax.legend(loc='lower right', fontsize=11)
    ax.grid(True, axis='x', alpha=0.3)

    if output_path is None:
        fig.tight_layout()
        plt.show()
        return

    _save(fig, output_path, "top predictions plot")


def plot_auc_comparison(
    results: Dict[str, Dict[str, float]],
    output_path: Path,
) -> None:
    """
    Plot AUC comparison across diseases.

    Parameters
    ----------
    results : dict
        {disease: {'auc_roc': ..., 'auc_pr': ...}}
    output_path : Path
        Output file path.

    Raises
    ------
    ValueError
        If `results` is empty.
    KeyError
        If any disease is missing a metric, naming the disease and the metric
        rather than surfacing a bare KeyError from inside a comprehension.
    """
    if not results:
        raise ValueError("results is empty, there is nothing to compare.")

    required = ('auc_roc', 'auc_pr')
    incomplete = {
        disease: [metric for metric in required if metric not in metrics]
        for disease, metrics in results.items()
        if any(metric not in metrics for metric in required)
    }

    if incomplete:
        raise KeyError(f"results is missing metrics: {incomplete}")

    diseases = list(results)
    auc_roc = [results[d]['auc_roc'] for d in diseases]
    auc_pr = [results[d]['auc_pr'] for d in diseases]

    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.arange(len(diseases))
    width = 0.35

    ax.bar(x - width / 2, auc_roc, width, label='AUC-ROC', alpha=0.8)
    ax.bar(x + width / 2, auc_pr, width, label='AUC-PR', alpha=0.8)

    ax.set_ylabel('AUC', fontsize=12)
    ax.set_title('Model Performance Across Diseases', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(diseases, rotation=45, ha='right')
    ax.legend(fontsize=11)
    ax.grid(True, axis='y', alpha=0.3)
    ax.set_ylim([0, 1.0])

    _save(fig, output_path, "AUC comparison")
