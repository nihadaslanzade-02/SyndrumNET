"""
Evaluate predictions against known synergies.

Computes AUC-ROC, AUC-PR, and other metrics.
"""

import argparse
import logging
from pathlib import Path

import pandas as pd

from syndrumnet.eval.benchmarks import load_known_synergies
from syndrumnet.eval.metrics import evaluate_predictions, compute_roc_curve, compute_precision_recall_curve
from syndrumnet.eval.reporting import generate_evaluation_report
from syndrumnet.viz.plots import plot_roc_curve, plot_pr_curve, plot_auc_comparison
from syndrumnet.utils.config import load_config
from syndrumnet.utils.logging import setup_logger


def main():
    parser = argparse.ArgumentParser(description="Evaluate SyndrumNET predictions")
    parser.add_argument('--config', type=str, required=True, help="Config file path")
    parser.add_argument('--synergy-file', type=str, help="Known synergies file")
    args = parser.parse_args()
    
    # Load config
    config = load_config(args.config)
    
    # Setup
    logger = setup_logger('evaluation', Path('logs'))
    
    logger.info("="*60)
    logger.info("SyndrumNET Evaluation")
    logger.info("="*60)
    
    # Load known synergies
    if args.synergy_file:
        synergy_file = Path(args.synergy_file)
    else:
        synergy_file = Path('data/raw/known_synergies.csv')
    
    if not synergy_file.exists():
        logger.error(f"Synergy file not found: {synergy_file}")
        logger.info("Please provide known synergies file via --synergy-file")
        return
    
    known_synergies = load_known_synergies(synergy_file)
    
    # Evaluate each disease
    results = {}
    figures_dir = Path('reports/figures')
    figures_dir.mkdir(parents=True, exist_ok=True)
    
    predictions_dir = Path('reports/tables')
    
    for disease in config.diseases:
        pred_file = predictions_dir / f"predictions_{disease.lower().replace(' ', '_')}.csv"
        
        if not pred_file.exists():
            logger.warning(f"Predictions not found for {disease}: {pred_file}")
            continue
        
        logger.info(f"\nEvaluating {disease}...")
        
        # Load predictions
        predictions = pd.read_csv(pred_file)
        
        # Evaluate
        metrics = evaluate_predictions(predictions, known_synergies)
        results[disease] = metrics
        
        # Plot ROC curve
        y_true = []
        y_score = []
        for _, row in predictions.iterrows():
            pair = tuple(sorted([row['drug_a'], row['drug_b']]))
            y_true.append(1 if pair in known_synergies else 0)
            y_score.append(row['prediction_score'])
        
        import numpy as np
        y_true = np.array(y_true)
        y_score = np.array(y_score)
        
        if len(np.unique(y_true)) > 1:
            fpr, tpr, _ = compute_roc_curve(y_true, y_score)
            plot_roc_curve(
                fpr, tpr, metrics['auc_roc'],
                figures_dir / f"roc_{disease.lower().replace(' ', '_')}.png",
                title=f"ROC Curve - {disease}"
            )
            
            precision, recall, _ = compute_precision_recall_curve(y_true, y_score)
            plot_pr_curve(
                precision, recall, metrics['auc_pr'],
                figures_dir / f"pr_{disease.lower().replace(' ', '_')}.png",
                title=f"Precision-Recall - {disease}"
            )
    
    # Generate summary report
    if results:
        report_file = Path('reports/tables/evaluation_summary.csv')
        generate_evaluation_report(results, report_file)
        
        # Plot comparison
        plot_auc_comparison(results, figures_dir / 'auc_comparison.png')
    
    logger.info("\n" + "="*60)
    logger.info("Evaluation complete!")
    logger.info("="*60 + "\n")


if __name__ == '__main__':
    main()
