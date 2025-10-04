"""
Run complete SyndrumNET prediction pipeline.

Computes synergy scores for all drug pairs across specified diseases.
"""

import argparse
import logging
from pathlib import Path

import pandas as pd

from syndrumnet.data.network_builder import NetworkBuilder
from syndrumnet.data.modules import ModuleBuilder
from syndrumnet.scoring.predictor import SynergyPredictor
from syndrumnet.utils.config import load_config
from syndrumnet.utils.logging import setup_logger
from syndrumnet.utils.seeds import set_random_seed


def load_modules(filepath: Path) -> dict:
    """Load modules from CSV."""
    df = pd.read_csv(filepath)
    
    modules = {}
    for name in df['module'].unique() if 'module' in df.columns else df['drug'].unique():
        col = 'module' if 'module' in df.columns else 'drug'
        genes = set(df[df[col] == name]['gene'])
        modules[name] = genes
    
    return modules


def load_drug_modules(filepath: Path) -> dict:
    """Load drug modules with up/down directions."""
    df = pd.read_csv(filepath)
    
    modules = {}
    for drug in df['drug'].unique():
        drug_df = df[df['drug'] == drug]
        modules[drug] = {
            'up': set(drug_df[drug_df['direction'] == 'up']['gene']),
            'down': set(drug_df[drug_df['direction'] == 'down']['gene']),
        }
    
    return modules


def main():
    parser = argparse.ArgumentParser(description="Run SyndrumNET pipeline")
    parser.add_argument('--config', type=str, required=True, help="Config file path")
    parser.add_argument('--diseases', nargs='+', help="Diseases to process (default: all)")
    parser.add_argument('--max-pairs', type=int, help="Max drug pairs per disease (for testing)")
    args = parser.parse_args()
    
    # Load config
    config = load_config(args.config)
    
    # Setup
    logger = setup_logger('pipeline', Path('logs'))
    set_random_seed(config.random_seed)
    
    logger.info("="*60)
    logger.info("SyndrumNET Pipeline")
    logger.info("="*60)
    
    # Determine diseases to process
    diseases = args.diseases if args.diseases else config.diseases
    
    logger.info(f"Processing diseases: {diseases}")
    
    # Load data
    logger.info("\n[1/3] Loading data...")
    
    processed_dir = Path('data/processed')
    
    # Load network
    G = NetworkBuilder.load(processed_dir / 'network.graphml')
    logger.info(f"Loaded network: {G.number_of_nodes()} nodes")
    
    # Load modules
    disease_modules = load_modules(processed_dir / 'disease_modules.csv')
    logger.info(f"Loaded {len(disease_modules)} disease modules")
    
    drug_modules = load_drug_modules(processed_dir / 'drug_modules.csv')
    logger.info(f"Loaded {len(drug_modules)} drug modules")
    
    # Initialize predictor
    logger.info("\n[2/3] Initializing predictor...")
    predictor = SynergyPredictor(
        G,
        n_randomizations=config.scoring.n_randomizations,
        seed=config.random_seed,
    )
    
    predictor.set_disease_modules(disease_modules)
    predictor.set_drug_modules(drug_modules)
    
    # Run predictions
    logger.info("\n[3/3] Computing predictions...")
    
    output_dir = Path('reports/tables')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for disease in diseases:
        logger.info(f"\nProcessing: {disease}")
        
        try:
            predictions = predictor.predict_all(disease, max_pairs=args.max_pairs)
            
            # Save results
            output_file = output_dir / f"predictions_{disease.lower().replace(' ', '_')}.csv"
            predictor.save_predictions(predictions, output_file)
            
            logger.info(f"Saved {len(predictions)} predictions to {output_file}")
            
        except Exception as e:
            logger.error(f"Failed to process {disease}: {e}", exc_info=True)
    
    logger.info("\n" + "="*60)
    logger.info("Pipeline complete!")
    logger.info(f"Results saved to: {output_dir}")
    logger.info("="*60 + "\n")


if __name__ == '__main__':
    main()
