"""
Evaluation metrics for synergy predictions.

Computes AUC-ROC, AUC-PR, and other classification metrics.
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    precision_recall_curve,
    roc_curve,
)

logger = logging.getLogger(__name__)


def compute_auc(
    y_true: np.ndarray,
    y_score: np.ndarray,
) -> float:
    """
    Compute AUC-ROC.
    
    Parameters
    ----------
    y_true : np.ndarray
        True binary labels (0/1).
    y_score : np.ndarray
        Predicted scores.
        
    Returns
    -------
    float
        AUC-ROC score.
    """
    if len(np.unique(y_true)) < 2:
        logger.warning("Only one class in y_true, cannot compute AUC")
        return np.nan
    
    auc = roc_auc_score(y_true, y_score)
    return auc


def compute_pr(
    y_true: np.ndarray,
    y_score: np.ndarray,
) -> float:
    """
    Compute AUC-PR (average precision).
    
    Parameters
    ----------
    y_true : np.ndarray
        True binary labels.
    y_score : np.ndarray
        Predicted scores.
        
    Returns
    -------
    float
        AUC-PR score.
    """
    if len(np.unique(y_true)) < 2:
        logger.warning("Only one class in y_true, cannot compute AUC-PR")
        return np.nan
    
    auc_pr = average_precision_score(y_true, y_score)
    return auc_pr


def compute_roc_curve(
    y_true: np.ndarray,
    y_score: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute ROC curve.
    
    Parameters
    ----------
    y_true : np.ndarray
        True labels.
    y_score : np.ndarray
        Predicted scores.
        
    Returns
    -------
    tuple
        (fpr, tpr, thresholds)
    """
    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    return fpr, tpr, thresholds


def compute_precision_recall_curve(
    y_true: np.ndarray,
    y_score: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute precision-recall curve.
    
    Parameters
    ----------
    y_true : np.ndarray
        True labels.
    y_score : np.ndarray
        Predicted scores.
        
    Returns
    -------
    tuple
        (precision, recall, thresholds)
    """
    precision, recall, thresholds = precision_recall_curve(y_true, y_score)
    return precision, recall, thresholds


def evaluate_predictions(
    predictions: pd.DataFrame,
    known_synergies: Set[Tuple[str, str]],
) -> Dict[str, float]:
    """
    Evaluate predictions against known synergies.
    
    Parameters
    ----------
    predictions : pd.DataFrame
        Predictions with 'drug_a', 'drug_b', 'prediction_score'.
    known_synergies : set
        Known synergistic pairs.
        
    Returns
    -------
    dict
        Evaluation metrics.
    """
    # Create labels
    y_true = []
    y_score = []
    
    for _, row in predictions.iterrows():
        pair = tuple(sorted([row['drug_a'], row['drug_b']]))
        label = 1 if pair in known_synergies else 0
        
        y_true.append(label)
        y_score.append(row['prediction_score'])
    
    y_true = np.array(y_true)
    y_score = np.array(y_score)
    
    # Compute metrics
    metrics = {
        'auc_roc': compute_auc(y_true, y_score),
        'auc_pr': compute_pr(y_true, y_score),
        'n_predictions': len(predictions),
        'n_known_synergies': len(known_synergies),
        'n_true_positives': int(np.sum(y_true)),
    }
    
    logger.info(f"Evaluation: AUC-ROC={metrics['auc_roc']:.3f}, AUC-PR={metrics['auc_pr']:.3f}")
    
    return metrics
