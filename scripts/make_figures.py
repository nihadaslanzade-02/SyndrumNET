"""
Generate all visualization figures.
"""

import argparse
import logging
from pathlib import Path

import pandas as pd

from syndrumnet.data.network_builder import NetworkBuilder
from syndrumnet.viz.plots import (
    plot_degree_distribution,
    plot_score_distributions,
    plot_top_predictions,
)
from syndrumnet.utils.config import load_config
from syndrumnet.utils.logging import setup_logger


def main():
    parser = argparse.ArgumentParser(description="Generate SyndrumNET figures")
    parser.add_argument('--config', type=str, required=True, help="Config file path")
    args = parser.parse_args()
    
    # Load config
    config = load_config(args.config)
    
    # Setup
    logger = setup_logger('figures', Path('logs'))
    
    logger.info("="*60)
    logger.info("SyndrumNET Figure Generation")
    logger.info("="*60)
    
    figures_dir = Path('reports/figures')
    figures_dir.mkdir(parents=True, exist_ok=True)
    
    # Network figures
    logger.info("\n[1/3] Generating network figures...")
    
    network_file = Path('data/processed/network.graphml')
    if network_file.exists():
        G = NetworkBuilder.load(network_file)
        plot_degree_distribution(G, figures_dir / 'degree_distribution.png')
    
    # Score distribution figures
    logger.info("\n[2/3] Generating score distribution figures...")
    
    predictions_dir = Path('reports/tables')
    for disease in config.diseases:
        pred_file = predictions_dir / f"predictions_{disease.lower().replace(' ', '_')}.csv"
        
        if pred_file.exists():
            predictions = pd.read_csv(pred_file)
            
            # Score distributions
            plot_score_distributions(
                predictions,
                figures_dir / f"score_dist_{disease.lower().replace(' ', '_')}.png"
            )
            
            # Top predictions
            plot_top_predictions(
                predictions,
                k=20,
                output_path=figures_dir / f"top_pairs_{disease.lower().replace(' ', '_')}.png"
            )
    
    logger.info("\n[3/3] Figure generation complete!")
    logger.info(f"Figures saved to: {figures_dir}")
    logger.info("="*60 + "\n")


if __name__ == '__main__':
    main()
