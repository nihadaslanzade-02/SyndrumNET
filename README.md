"""

\# SyndrumNET Reproduction



A complete, reproducible Python implementation of the SyndrumNET model for 

predicting synergistic drug combinations using network propagation and 

multi-omics integration.



\## Citations



This implementation reproduces:



1\. \*\*Rahim Rahimov (2025)\*\*. "Reproducing and Extending the SyndrumNET Model 

&nbsp;  for Synergistic Drug Combination Prediction". Master's Thesis.



2\. \*\*Iida et al. (2024)\*\*. "A network-based trans-omics approach for predicting 

&nbsp;  synergistic drug combinations (SyndrumNET)". Communications Medicine.



For network propagation formalism:

3\. \*\*Cowen et al. (2017)\*\*. "Network propagation: a universal amplifier of 

&nbsp;  genetic associations". Nature Reviews Genetics.



\## Overview



SyndrumNET predicts synergistic drug combinations by integrating:

\- \*\*Topological analysis\*\* (TQAB): Network-based complementarity of drug targets

\- \*\*Proximity scoring\*\* (PQAB): Disease-drug module distances  

\- \*\*Transcriptional correlation\*\* (CQAB): Expression signature similarity



The final prediction score combines these components:

```

Score\_Q,AB = TQAB + PQAB + CQAB

```



\## Installation



\### Prerequisites

\- Python 3.10+

\- conda or mamba

\- 16GB+ RAM recommended

\- 50GB disk space for data



\### Setup Environment



```bash

\# Clone repository

git clone https://github.com/yourusername/SyndrumNET-repro.git

cd SyndrumNET-repro



\# Create conda environment

conda env create -f environment.yml

conda activate syndrumnet



\# Install package

pip install -e .



\# Verify installation

python -c "import syndrumnet; print(syndrumnet.\_\_version\_\_)"

```



\## Data Sources \& Licenses



This pipeline downloads and integrates data from:



\- \*\*Molecular interactions\*\*: HuRI, CORUM, PhosphoSitePlus, KEGG RPair, 

&nbsp; SignaLink, InnateDB, Instruct

\- \*\*Disease expression\*\*: CREEDS (79 diseases)

\- \*\*Drug expression\*\*: LINCS L1000 (1,488 drugs)

\- \*\*Disease genes\*\*: OMIM, ClinVar, GWAS, PheWAS, GWASdb, DisGeNET

\- \*\*Drug structures\*\*: PubChem, ChEMBL (for KCF-S fingerprints)



See `docs/DATA\_SOURCES.md` for detailed source information and licenses.



\## Pipeline Diagram



```

Data Acquisition → Preprocessing → Feature Computation → Propagation → Prediction

&nbsp;     ↓                 ↓                  ↓                ↓             ↓

&nbsp; Downloads        ID Mapping        Distances         PRINCE      Final Scores

&nbsp; Parsers          Modules           Proximities    Similarities    Rankings

&nbsp;                  Network           Correlations                   Evaluation

```



\## Quickstart



\### 1. Build Data (first time only, ~2-4 hours)



```bash

python scripts/build\_all\_data.py --config configs/default.yaml

```



This downloads all data sources, performs ID mapping, builds the integrated 

network, and constructs disease/drug modules.



\### 2. Run Full Pipeline



```bash

python scripts/run\_pipeline.py --config configs/default.yaml

```



Computes all scores (TQAB, PQAB, CQAB) for drug pairs across diseases 

(asthma, diabetes T2, hypertension, colorectal cancer, AML, CML).



\### 3. Evaluate Results



```bash

python scripts/evaluate.py --config configs/default.yaml

```



Computes AUC/PR metrics against known synergy resources and generates 

evaluation tables.



\### 4. Generate Figures



```bash

python scripts/make\_figures.py --config configs/default.yaml

```



Creates all matplotlib visualizations in `reports/figures/`.



\### 5. View Results



```bash

\# Top predictions per disease

ls reports/tables/predictions\_\*.csv



\# Evaluation metrics

cat reports/tables/evaluation\_summary.csv



\# Figures

ls reports/figures/

```



\## Configuration Guide



All parameters are controlled via YAML configs in `configs/`:



```yaml

\# configs/default.yaml

random\_seed: 42

n\_cores: 4



data:

&nbsp; network\_sources:

&nbsp;   - huri

&nbsp;   - corum

&nbsp;   - phosphositeplus

&nbsp;   # ...

&nbsp; 

propagation:

&nbsp; alpha: 0.5          # Restart probability

&nbsp; tolerance: 1e-6     # Convergence threshold

&nbsp; max\_iterations: 1000



scoring:

&nbsp; n\_randomizations: 1000  # For z-score null models

&nbsp; top\_pct\_genes: 0.05     # L1000 module definition (5%)

&nbsp; 

diseases:

&nbsp; - asthma

&nbsp; - diabetes\_t2

&nbsp; - hypertension

&nbsp; - colorectal\_cancer

&nbsp; - aml

&nbsp; - cml

```



Override config values via CLI:

```bash

python scripts/run\_pipeline.py --config configs/default.yaml \\

&nbsp;   --propagation.alpha 0.7 --scoring.n\_randomizations 2000

```



\## Extending to New Diseases



Add a new disease from GEO expression data:



```bash

python scripts/add\_new\_disease.py \\

&nbsp;   --name "Breast Cancer" \\

&nbsp;   --geo GSE12345 \\

&nbsp;   --config configs/default.yaml

```



This:

1\. Downloads GEO dataset

2\. Builds CREEDS-style up/down gene signatures

3\. Identifies disease susceptibility genes from OMIM/ClinVar/etc.

4\. Creates disease module (union of signatures + susceptibility genes)

5\. Computes all scores for this disease

6\. Saves results to `reports/tables/predictions\_breast\_cancer.csv`



\## Extending to New Drugs



Add new drugs from LINCS L1000:



```bash

\# drugs.txt contains LINCS compound IDs, one per line

python scripts/add\_new\_drugs.py \\

&nbsp;   --lincs drugs.txt \\

&nbsp;   --config configs/default.yaml

```



This:

1\. Fetches L1000 expression profiles

2\. Harmonizes cell lines per thesis specifications

3\. Defines drug modules (top 5% up/down genes by fold-change)

4\. Integrates into existing pipeline

5\. Recomputes predictions with new drugs included



\## Troubleshooting



\### Out of Memory

\- Reduce `n\_cores` in config

\- Process diseases sequentially: `--diseases asthma --diseases diabetes\_t2`

\- Use sparse matrix operations (default in this implementation)



\### Download Failures

\- Check network connectivity

\- Some sources require registration (e.g., PhosphoSitePlus)

\- Retry: `python scripts/build\_all\_data.py --retry-failed`



\### Reproducibility Issues

\- Ensure same `random\_seed` in config

\- Check conda environment matches `conda-lock.yml`

\- Verify data versions in `data/raw/VERSIONS.txt`



\### Missing Dependencies

```bash

conda env update -f environment.yml

pip install -e .

```



\## Reproducibility Notes



\- All random operations are seeded via `configs/default.yaml:random\_seed`

\- Data versions are logged in `data/raw/VERSIONS.txt`

\- Intermediate artifacts cached in `data/processed/`

\- Run logs saved to `logs/syndrumnet\_YYYYMMDD\_HHMMSS.log`

\- Repeated runs with same seed produce identical outputs



To ensure exact reproduction:

```bash

\# Use locked dependencies

conda create --name syndrumnet --file conda-lock.yml

conda activate syndrumnet

pip install -e .



\# Use versioned data (optional)

python scripts/build\_all\_data.py --use-archived-versions

```



\## Development



\### Running Tests



```bash

\# All tests

pytest tests/ -v



\# Specific test module

pytest tests/test\_distances.py -v



\# With coverage

pytest tests/ --cov=syndrumnet --cov-report=html

```



\### Code Style



```bash

\# Format code

black src/ scripts/ tests/



\# Type checking

mypy src/



\# Linting

flake8 src/ scripts/ tests/

```



\## Documentation



\- `docs/METHOD\_NOTES.md`: Mathematical details and implementation choices

\- `docs/DATA\_SOURCES.md`: Data provenance and licenses  

\- `docs/API.md`: API reference for all modules



\## License



MIT License - see LICENSE file



\## Contact



For questions about this implementation, open an issue on GitHub.

For questions about the original method, refer to the Iida et al. (2024) paper.

