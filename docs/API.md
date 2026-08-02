# API Reference

Every public class and function in `syndrumnet`, grouped by subpackage in the
order the pipeline uses them: load and parse data, build the network and
modules, compute metrics, propagate, score, evaluate, visualise.

This page is a map. Full parameter and return documentation lives in the
docstrings themselves, in NumPy style, and is the authoritative source:

```python
from syndrumnet.propagation.prince import PRINCE
help(PRINCE)
```

Regenerate this file after changing the public surface:

```bash
python scripts/gen_api_docs.py > docs/API.md
```

---

## `syndrumnet.io`

### `syndrumnet.io.downloaders`

- **`DataDownloader`** - Centralized data downloader for SyndrumNET pipeline.
  - methods: `download_file()`, `download_huri()`, `download_corum()`, `download_phosphositeplus()`, `download_creeds()`, `download_lincs()`, `download_disease_genes()`, `download_id_mapping()`, `download_all()`, `get_file_paths()`

### `syndrumnet.io.id_mapping`

- **`IDMapper`** - Gene/protein ID mapping service.
  - methods: `to_hgnc()`, `to_entrez()`, `harmonize_gene_list()`, `batch_convert()`

### `syndrumnet.io.parsers`

- `parse_huri(filepath)` - Parse HuRI protein-protein interaction data.
- `parse_corum(filepath)` - Parse CORUM protein complex data.
- `parse_phosphositeplus(filepath)` - Parse PhosphoSitePlus kinase-substrate data.
- `parse_kegg_rpair(filepath)` - Parse KEGG RPair reaction data.
- `parse_creeds(filepath)` - Parse CREEDS disease signatures.
- `parse_lincs(sig_filepath, meta_filepath, top_pct)` - Parse LINCS L1000 drug signatures.

## `syndrumnet.data`

### `syndrumnet.data.modules`

- **`ModuleBuilder`** - Construct disease and drug modules.
  - methods: `build_disease_modules()`, `build_drug_modules()`, `save_modules()`, `load_modules()`

### `syndrumnet.data.network_builder`

- **`NetworkBuilder`** - Build integrated human molecular interaction network.
  - methods: `add_source()`, `build()`, `get_network_stats()`, `save()`, `load()`

## `syndrumnet.metrics`

### `syndrumnet.metrics.distances`

- `shortest_path_distance(G, source_set, target_set, infinity_value)` - Compute average shortest path distance from source set to target set.
- `module_proximity(G, module_a, module_b)` - Compute bidirectional proximity between two modules.
- `separation_score(G, module_a, module_b)` - Compute network separation s_AB between two modules.
- `compute_all_pairwise_distances(G, gene_set)` - Compute all pairwise shortest path distances within a gene set.

### `syndrumnet.metrics.null_models`

- `degree_preserving_randomization(G, module, n_random, seed)` - Generate degree-preserving random gene sets.
- `compute_zscore(observed, null_distribution)` - Compute z-score of observed value against null distribution.
- `compute_normalized_proximity(G, disease_module, drug_module, n_random, seed)` - Compute z-score normalized proximity between disease and drug modules.

### `syndrumnet.metrics.transcription`

- `compute_correlation(signature_a, signature_b, method)` - Compute correlation between two gene expression signatures.
- `transcriptional_similarity(disease_signature, drug_signature_up, drug_signature_down, inverse_correlation)` - Compute transcriptional similarity between disease and drug.
- `aggregate_transcriptional_scores(scores, method)` - Aggregate multiple transcriptional scores.

## `syndrumnet.propagation`

### `syndrumnet.propagation.prince`

- **`PRINCE`** - PRINCE network propagation algorithm.
  - methods: `propagate()`, `propagate_multiple()`, `get_top_genes()`

### `syndrumnet.propagation.similarity_layers`

- `compute_disease_similarity(disease_modules, method)` - Compute pairwise disease similarity matrix.
- `compute_drug_similarity(drug_fingerprints, method)` - Compute pairwise drug similarity matrix.
- `kcf_fingerprint_similarity(smiles_a, smiles_b)` - Compute KCF-S fingerprint similarity between two molecules.
- `jaccard_similarity(set_a, set_b)` - Jaccard similarity coefficient.
- `overlap_coefficient(set_a, set_b)` - Overlap coefficient.
- `tanimoto_similarity_matrix(fingerprints)` - Compute pairwise Tanimoto similarity for binary fingerprints.
- `build_similarity_matrix(entities, similarity_func)` - Build pairwise similarity matrix for a list of entities.

## `syndrumnet.scoring`

### `syndrumnet.scoring.cqab`

- `compute_cqab(disease_signature, drug_a_signature_up, drug_a_signature_down, drug_b_signature_up, drug_b_signature_down)` - Compute transcriptional correlation score CQAB.
- `compute_cqab_batch(disease_signature, drug_signatures, drug_pairs)` - Compute CQAB for multiple drug pairs.

### `syndrumnet.scoring.pqab`

- `compute_pqab(G, disease_module, drug_a_module, drug_b_module, n_randomizations, seed)` - Compute proximity score PQAB.
- `compute_pqab_batch(G, disease_module, drug_modules, drug_pairs, n_randomizations, seed)` - Compute PQAB for multiple drug pairs.

### `syndrumnet.scoring.predictor`

- **`SynergyPredictor`** - Complete SyndrumNET synergy prediction pipeline.
  - methods: `set_disease_modules()`, `set_drug_modules()`, `set_disease_signatures()`, `predict_all()`, `predict_multiple_diseases()`, `save_predictions()`

### `syndrumnet.scoring.tqab`

- **`TopologyClass`** - Enumeration of drug pair topology classes.
- `compute_tqab(G, disease_module, drug_a_module, drug_b_module)` - Compute topological class score TQAB.
- `compute_tqab_batch(G, disease_module, drug_modules, drug_pairs)` - Compute TQAB for multiple drug pairs.

## `syndrumnet.eval`

### `syndrumnet.eval.benchmarks`

- `load_known_synergies(filepath, disease_filter)` - Load known synergistic drug combinations.
- `load_drugcombdb(filepath)` - Load DrugCombDB database.

### `syndrumnet.eval.metrics`

- `compute_auc(y_true, y_score)` - Compute AUC-ROC.
- `compute_pr(y_true, y_score)` - Compute AUC-PR (average precision).
- `compute_roc_curve(y_true, y_score)` - Compute ROC curve.
- `compute_precision_recall_curve(y_true, y_score)` - Compute precision-recall curve.
- `evaluate_predictions(predictions, known_synergies)` - Evaluate predictions against known synergies.

### `syndrumnet.eval.reporting`

- `generate_evaluation_report(results, output_path)` - Generate evaluation summary report.

## `syndrumnet.viz`

### `syndrumnet.viz.plots`

- `plot_degree_distribution(G, output_path, log_scale)` - Plot network degree distribution.
- `plot_roc_curve(fpr, tpr, auc, output_path, title)` - Plot ROC curve.
- `plot_pr_curve(precision, recall, auc_pr, output_path, title)` - Plot precision-recall curve.
- `plot_score_distributions(predictions, output_path)` - Plot distributions of TQAB, PQAB, CQAB scores.
- `plot_top_predictions(predictions, k, output_path)` - Plot top-k predictions with stacked score components.
- `plot_auc_comparison(results, output_path)` - Plot AUC comparison across diseases.

## `syndrumnet.utils`

### `syndrumnet.utils.config`

- **`Config`** - Configuration container with nested access support.
  - methods: `get()`, `to_dict()`
- `load_config(config_path)` - Load configuration from YAML file.
- `merge_configs(base, override)` - Merge override dictionary into base config.
- `save_config(config, output_path)` - Save configuration to YAML file.

### `syndrumnet.utils.logging`

- `setup_logger(name, log_dir, level, console)` - Setup structured logger with file and console handlers.
- **`LoggerMixin`** - Mixin to add logger attribute to classes.
  - methods: `logger()`

### `syndrumnet.utils.seeds`

- `set_random_seed(seed)` - Set random seed for all libraries to ensure reproducibility.
- `get_random_state(seed)` - Create a NumPy RandomState object for isolated random operations.

---

## Entry points

The four scripts in `scripts/` wire the above together. Each takes
`--config <path>` and accepts dotted overrides such as
`--propagation.alpha 0.7`.

| Script | Stage |
|---|---|
| `build_all_data.py` | Download sources, map IDs, build the network and modules |
| `run_pipeline.py` | Compute TQAB, PQAB and CQAB for every drug pair |
| `evaluate.py` | Score predictions against known synergies (AUC-ROC, AUC-PR) |
| `make_figures.py` | Render figures into `reports/figures/` |

