"""
Visualization functions using matplotlib.

All plots use matplotlib only (no seaborn) as specified.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


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
    """
    degrees = [d for n, d in G.degree()]
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.hist(degrees, bins=50, edgecolor='black', alpha=0.7)
    ax.set_xlabel('Degree', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Network Degree Distribution', fontsize=14, fontweight='bold')
    
    if log_scale:
        ax.set_xscale('log')
        ax.set_yscale('log')
    
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved degree distribution to {output_path}")


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
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved ROC curve to {output_path}")


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
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved PR curve to {output_path}")


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
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    scores = ['tqab', 'pqab', 'cqab']
    titles = ['TQAB (Topological)', 'PQAB (Proximity)', 'CQAB (Transcriptional)']
    
    for ax, score, title in zip(axes, scores, titles):
        if score in predictions.columns:
            ax.hist(predictions[score].dropna(), bins=30, edgecolor='black', alpha=0.7)
            ax.set_xlabel('Score', fontsize=11)
            ax.set_ylabel('Count', fontsize=11)
            ax.set_title(title, fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved score distributions to {output_path}")


def plot_top_predictions(
    predictions: pd.DataFrame,
    k: int = 20,
    output_path: Optional[Path] = None,
) -> None:
    """
    Plot top-k predictions with stacked score components.
    
    Parameters
    ----------
    predictions : pd.DataFrame
        Predictions dataframe.
    k : int
        Number of top predictions to show.
    output_path : Path, optional
        Output file path.
    """
    top_k = predictions.nlargest(k, 'prediction_score')
    
    # Create labels
    labels = [f"{row['drug_a']}-{row['drug_b']}" for _, row in top_k.iterrows()]
    
    # Extract score components
    tqab = top_k['tqab'].values
    pqab = top_k['pqab'].values
    cqab = top_k['cqab'].values
    
    # Create stacked bar chart
    fig, ax = plt.subplots(figsize=(12, 8))
    
    x = np.arange(len(labels))
    width = 0.8
    
    ax.barh(x, tqab, width, label='TQAB', alpha=0.8)
    ax.barh(x, pqab, width, left=tqab, label='PQAB', alpha=0.8)
    ax.barh(x, cqab, width, left=tqab+pqab, label='CQAB', alpha=0.8)
    
    ax.set_yticks(x)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel('Score', fontsize=12)
    ax.set_title(f'Top {k} Predicted Synergistic Combinations', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right', fontsize=11)
    ax.grid(True, axis='x', alpha=0.3)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved top predictions plot to {output_path}")
    else:
        plt.show()


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
    """
    diseases = list(results.keys())
    auc_roc = [results[d]['auc_roc'] for d in diseases]
    auc_pr = [results[d]['auc_pr'] for d in diseases]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x = np.arange(len(diseases))
    width = 0.35
    
    ax.bar(x - width/2, auc_roc, width, label='AUC-ROC', alpha=0.8)
    ax.bar(x + width/2, auc_pr, width, label='AUC-PR', alpha=0.8)
    
    ax.set_ylabel('AUC', fontsize=12)
    ax.set_title('Model Performance Across Diseases', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(diseases, rotation=45, ha='right')
    ax.legend(fontsize=11)
    ax.grid(True, axis='y', alpha=0.3)
    ax.set_ylim([0, 1.0])
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved AUC comparison to {output_path}")