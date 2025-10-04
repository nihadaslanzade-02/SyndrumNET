"""Evaluation: benchmarks, metrics, and reporting."""

from syndrumnet.eval.benchmarks import load_known_synergies
from syndrumnet.eval.metrics import compute_auc, compute_pr
from syndrumnet.eval.reporting import generate_evaluation_report

__all__ = [
    "load_known_synergies",
    "compute_auc",
    "compute_pr",
    "generate_evaluation_report",
]
