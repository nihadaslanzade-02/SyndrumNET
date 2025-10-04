"""
Reporting utilities for evaluation results.

Generates summary tables and reports.
"""

import logging
from pathlib import Path
from typing import Dict, List

import pandas as pd

logger = logging.getLogger(__name__)


def generate_evaluation_report(
    results: Dict[str, Dict[str, float]],
    output_path: Path,
) -> None:
    """
    Generate evaluation summary report.
    
    Parameters
    ----------
    results : dict
        {disease: metrics_dict}
    output_path : Path
        Output CSV path.
    """
    logger.info("Generating evaluation report")
    
    rows = []
    for disease, metrics in results.items():
        row = {'disease': disease}
        row.update(metrics)
        rows.append(row)
    
    df = pd.DataFrame(rows)
    
    # Sort by AUC-ROC
    df = df.sort_values('auc_roc', ascending=False)
    
    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    
    logger.info(f"Saved evaluation report to {output_path}")
    
    # Print summary
    print("\n" + "="*60)
    print("EVALUATION SUMMARY")
    print("="*60)
    print(df.to_string(index=False))
    print("="*60 + "\n")
