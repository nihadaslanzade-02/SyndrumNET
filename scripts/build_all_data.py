"""
Download and preprocess all required data.

This script:
1. Downloads molecular interaction data
2. Downloads disease/drug expression data
3. Builds the integrated network
4. Constructs disease and drug modules
5. Saves all processed data for pipeline execution
"""

import argparse
import logging
from pathlib import Path

from syndrumnet.io.downloaders import DataDownloader
from syndrumnet.io.id_mapping import IDMapper
from syndrumnet.data.network_builder import NetworkBuilder
from syndrumnet.data.modules import ModuleBuilder
from syndrumnet.utils.config import load_config
from syndrumnet.utils.logging import setup_logger
from syndrumnet.utils.seeds import set_random_seed


def main():
    parser = argparse.ArgumentParser(description="Build all data for SyndrumNET")
    parser.add_argument('--config', type=str, required=True, help="Config file path")
    parser.add_argument('--retry-failed', action='store_true', help="Retry failed downloads")
    args = parser.parse_args()
    
    # Load config
    config = load_config(args.config)
    
    # Setup
    logger = setup_logger('build_data', Path('logs'))
    set_random_seed(config.random_seed)
    
    logger.info("="*60)
    logger.info("SyndrumNET Data Builder")
    logger.info("="*60)
    
    # Create directories
    data_dir = Path('data')
    raw_dir = data_dir / 'raw'
    interim_dir = data_dir / 'interim'
    processed_dir = data_dir / 'processed'
    
    for d in [raw_dir, interim_dir, processed_dir]:
        d.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Download data
    logger.info("\n[1/5] Downloading data sources...")
    downloader = DataDownloader(raw_dir)
    files = downloader.download_all()
    
    # Step 2: ID mapping
    logger.info("\n[2/5] Setting up ID mapper...")
    mapper = IDMapper(hgnc_file=raw_dir / 'hgnc.txt')
    
    # Step 3: Build network
    logger.info("\n[3/5] Building integrated network...")
    builder = NetworkBuilder(mapper)
    
    # Add sources
    for source in config.data.network_sources:
        if source in files and files[source].exists():
            builder.add_source(source, files[source])
    
    # Build and save
    G = builder.build()
    network_file = processed_dir / 'network.graphml'
    builder.save(network_file)
    
    # Log stats
    stats = builder.get_network_stats()
    logger.info(f"Network statistics: {stats}")
    
    # Step 4: Build disease modules
    logger.info("\n[4/5] Building disease modules...")
    module_builder = ModuleBuilder(G, mapper)
    
    if (raw_dir / 'creeds_disease.txt').exists():
        disease_modules = module_builder.build_disease_modules(
            raw_dir / 'creeds_disease.txt'
        )
        module_builder.save_modules(
            disease_modules,
            processed_dir / 'disease_modules.csv'
        )
    else:
        logger.warning("CREEDS data not found, skipping disease modules")
    
    # Step 5: Build drug modules
    logger.info("\n[5/5] Building drug modules...")
    
    lincs_sigs = raw_dir / 'lincs_signatures.txt'
    lincs_meta = raw_dir / 'lincs_metadata.txt'
    
    if lincs_sigs.exists() and lincs_meta.exists():
        drug_modules = module_builder.build_drug_modules(
            lincs_sigs,
            lincs_meta,
            top_pct=config.scoring.top_pct_genes,
        )
        
        # Save as CSV (flatten nested dict structure)
        import pandas as pd
        rows = []
        for drug, sigs in drug_modules.items():
            for gene in sigs['up']:
                rows.append({'drug': drug, 'gene': gene, 'direction': 'up'})
            for gene in sigs['down']:
                rows.append({'drug': drug, 'gene': gene, 'direction': 'down'})
        
        df = pd.DataFrame(rows)
        df.to_csv(processed_dir / 'drug_modules.csv', index=False)
        logger.info(f"Saved {len(drug_modules)} drug modules")
    else:
        logger.warning("LINCS data not found, skipping drug modules")
    
    logger.info("\n" + "="*60)
    logger.info("Data build complete!")
    logger.info(f"Processed data saved to: {processed_dir}")
    logger.info("="*60 + "\n")


if __name__ == '__main__':
    main()