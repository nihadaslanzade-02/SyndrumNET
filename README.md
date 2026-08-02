# SyndrumNET: Network-Based Prediction of Synergistic Drug Combinations

[![CI](https://github.com/nihadaslanzade-02/SyndrumNET/actions/workflows/ci.yml/badge.svg)](https://github.com/nihadaslanzade-02/SyndrumNET/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-3776AB)](https://www.python.org/)
[![Ruff](https://img.shields.io/badge/lint-ruff-261230)](https://docs.astral.sh/ruff/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

A Python implementation of **SyndrumNET**, a network-based trans-omics method that predicts which pairs of drugs will act synergistically against a given disease, by asking where each drug's effect lands on the human interactome relative to where the disease lives.

The pipeline integrates seven molecular-interaction resources into one network, builds gene modules for diseases and drugs from expression signatures, propagates signal across the network with PRINCE, and scores every drug pair on three independent axes: **topology**, **proximity**, and **transcriptional reversal**.

---

## The method

SyndrumNET scores a drug pair **(A, B)** against a disease **Q** by summing three components:

```
Score(Q,A,B)  =  T(Q,A,B)  +  P(Q,A,B)  +  C(Q,A,B)
                 topology     proximity    transcription
```

### Network distance and separation

All three components rest on one primitive, the average shortest-path distance from one gene module to another:

```
d(S,T) = (1/|S|) * sum_{s in S}  min_{t in T} dist(s,t)
```

From it comes the **separation score**, which asks whether two modules occupy distinct network neighbourhoods or overlap:

```
s_AB = <d_AB> - (<d_AA> + <d_BB>)/2
```

Positive `s_AB` means the modules are separated, negative means they overlap. Disconnected node pairs contribute a large finite value (1000.0) rather than infinity, so a fragmented network degrades gracefully instead of poisoning every average with `inf`.

### T: topological class

The intuition the whole method is built on: **two drugs help each other when they hit different parts of the disease neighbourhood, and duplicate each other when they hit the same part.** [`scoring/tqab.py`](src/syndrumnet/scoring/tqab.py) turns that into a three-way classification:

| Condition | Class | Score |
|---|---|---|
| `s_AB > 0` and both drugs within distance 3 of the disease | **Complementary** | `1 - d_mean/10` |
| `s_AB > 0`, at least one drug further away | **Intermediate** | `0.5 - d_mean/10` |
| `s_AB <= 0` (drug modules overlap) | **Redundant** | `-abs(s_AB)/5` |

### P: proximity

`P(Q,A,B) = (P_QA + P_QB)/2`, where each term is a disease-drug distance converted to a z-score against a null distribution:

```
z(d_QA) = (d_QA - mean(d_random)) / std(d_random)
```

The null is built by **degree-preserving randomisation**. Genes are binned into 20 degree strata and random modules are resampled within strata, so a "close" drug module cannot score well merely for containing hub genes. This is the standard construction from Guney et al. (2016).

### C: transcriptional reversal

`C(Q,A,B) = (C_QA + C_QB)/2`, comparing each drug's L1000 expression signature against the disease's CREEDS signature. Drug modules are the top 5% of genes by absolute fold change, median-aggregated across cell lines.

### PRINCE propagation

[`propagation/prince.py`](src/syndrumnet/propagation/prince.py) implements random walk with restart:

```
F^(t+1) = alpha * W * F^(t) + (1 - alpha) * F^(0)
```

with `alpha = 0.5` by default, iterating to a 1e-6 tolerance. The adjacency matrix is normalised column-wise (row-wise and symmetric are also selectable), kept sparse in CSR form, and built once per network so repeated propagations reuse it.

---

## Status

The library layer is tested and green. The **full pipeline has not yet been run end to end**, so no predictions or evaluation results are committed. Being specific, because the difference matters if you are about to clone this:

| Layer | State |
|---|---|
| Package import | Clean across Python 3.10, 3.11 and 3.12 |
| Test suite | **9 passing**: distances, null models, propagation, scoring |
| Lint | `ruff` clean, enforced in CI |
| Distances, separation, PRINCE, TQAB, CQAB, null models | Implemented and covered by tests |
| Data acquisition | Downloaders and parsers written; KEGG RPair and several disease-gene sources are explicit placeholders |
| ID harmonisation | Two competing calls in `network_builder.build()`, marked `TODO` in place |
| LINCS parsing | Metadata table is loaded but never joined, so drug keys are still raw signature IDs; marked `TODO` in place |
| End-to-end run | Not yet performed. No committed predictions, figures or evaluation tables |

Earlier revisions of this repository could not be imported at all: seven `typing` names and `pandas` were used without being imported, and because Python evaluates annotations at definition time, five of the eight subpackages raised `NameError` on import. That took every entry-point script down with them and left `pytest` collecting zero tests. Fixed in `650f52b`, along with the defect the test suite had been written to catch but never got to run: the degree-preserving null model called `_get_bin` with the wrong arity, and separately handed it bin contents where it expected numeric bin edges.

CI now guards against exactly that failure mode. Every module is imported in a dedicated job, so a missing import fails the build even when no test happens to exercise the affected function.

---

## Architecture

```
src/syndrumnet/
├── io/              downloaders, format parsers, HGNC/UniProt ID mapping
├── data/            network assembly from 7 interaction sources; disease & drug modules
├── metrics/         shortest-path distances, separation, degree-preserving null models,
│                    transcriptional correlation
├── propagation/     PRINCE random walk with restart; similarity layers
├── scoring/         tqab.py, pqab.py, cqab.py, and the predictor that combines them
├── eval/            AUC-ROC / AUC-PR against known synergy resources, reporting
├── viz/             matplotlib figures
└── utils/           YAML config loading with CLI override, logging, seeding

scripts/             build_all_data · run_pipeline · evaluate · make_figures · gen_api_docs
configs/             default.yaml + one YAML per disease
tests/               distances · null_models · propagation · scoring
```

Design decisions worth noting: `src/` layout so tests run against the installed package rather than the working directory; every public function carries a full NumPy-style docstring with parameters, returns and references; a single YAML config drives every stage with CLI overrides layered on top; seeding is centralised in `utils/seeds.py` rather than scattered across call sites; and scoring is split into one module per component so the three axes stay independently testable.

The package ships a `py.typed` marker, so type checkers resolve annotations for downstream code.

---

## Getting started

```bash
conda env create -f environment.yml
conda activate syndrumnet
pip install -e ".[dev]"
```

Requires Python 3.10+, roughly 16 GB RAM and ~50 GB of disk for the full data build.

```bash
pytest tests/ -v      # 9 tests, no data or network access needed
ruff check .          # same lint gate CI runs
```

The pipeline itself, in order:

```bash
python scripts/build_all_data.py --config configs/default.yaml   # downloads + network build
python scripts/run_pipeline.py   --config configs/default.yaml   # TQAB + PQAB + CQAB
python scripts/evaluate.py       --config configs/default.yaml   # AUC-ROC / AUC-PR
python scripts/make_figures.py   --config configs/default.yaml   # figures
```

or `make all`. Results land in `reports/tables/predictions_<disease>.csv` and `reports/figures/`. See [RUNBOOK.md](RUNBOOK.md) for expected runtimes and troubleshooting.

Using a single component directly, without any pipeline machinery:

```python
import networkx as nx
from syndrumnet.metrics.distances import separation_score
from syndrumnet.propagation.prince import PRINCE

G = nx.karate_club_graph()
G = nx.relabel_nodes(G, {n: f"GENE{n}" for n in G})

separation_score(G, {"GENE0", "GENE1"}, {"GENE32", "GENE33"})

prince = PRINCE(G, alpha=0.5)
scores = prince.propagate({"GENE0", "GENE1"})
prince.get_top_genes(scores, k=5)
```

## Configuration

Everything is driven by [`configs/default.yaml`](configs/default.yaml):

```yaml
random_seed: 42

propagation:
  alpha: 0.5              # restart probability
  tolerance: 1.0e-6
  normalize: 'column'     # 'column' | 'row' | 'symmetric'

scoring:
  n_randomizations: 1000  # null model draws for proximity z-scores
  top_pct_genes: 0.05     # top 5% |fold-change| defines a drug module
  weight_tqab: 1.0        # component weights, exposed for ablation studies
  weight_pqab: 1.0
  weight_cqab: 1.0

diseases: [asthma, diabetes_t2, hypertension, colorectal_cancer, aml, cml]
```

Any value can be overridden per run:

```bash
python scripts/run_pipeline.py --config configs/default.yaml \
    --propagation.alpha 0.7 --scoring.n_randomizations 2000
```

Six diseases ship with their own config in [`configs/diseases/`](configs/diseases/). `scripts/add_new_disease.py` and `scripts/add_new_drugs.py` are the intended extension points for adding a disease from a GEO series or drugs from LINCS; both are currently stubs.

---

## Data sources

| Layer | Sources |
|---|---|
| Molecular interactions | HuRI, CORUM, PhosphoSitePlus, KEGG RPair, SignaLink, InnateDB, Instruct |
| Disease expression | CREEDS (79 diseases) |
| Drug expression | LINCS L1000 (1,488 compounds) |
| Disease genes | OMIM, ClinVar, GWAS, PheWAS, GWASdb, DisGeNET |
| Drug structures | PubChem, ChEMBL |
| ID mapping | HGNC, UniProt |

No data is committed to this repository; `data/` is populated by `build_all_data.py`. Licences vary and several sources are non-commercial or require registration, so read [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) before redistributing anything derived from them.

**Deviations from the published method**, both documented in [`docs/METHOD_NOTES.md`](docs/METHOD_NOTES.md): Morgan fingerprints (RDKit) stand in for KCF-S fingerprints, and sources requiring registration or API keys are placeholders.

---

## Roadmap

1. **Resolve the two `TODO`s in the data layer.** `network_builder.build()` calls `harmonize_gene_list()` and discards the result, building its mapping from a separate `to_hgnc()` call instead; `parse_lincs()` loads the metadata table but never joins it, leaving raw signature IDs as drug keys. Both are marked in place.
2. **Fill in the placeholder parsers** for KEGG RPair and the disease-gene sources, so the network is the full seven-source integration the method specifies.
3. **Run the pipeline once on one disease**, with a reduced `n_randomizations`, to get a first end-to-end result and a committed prediction table.
4. **Add integration tests** on a small synthetic network that exercise `SynergyPredictor.predict_all()` end to end, not just the components.
5. **Validate against the reference results** in the source thesis, disease by disease.

---

## Documentation

| Document | Contents |
|---|---|
| [`docs/API.md`](docs/API.md) | Every public class and function, generated from source by `scripts/gen_api_docs.py` |
| [`docs/METHOD_NOTES.md`](docs/METHOD_NOTES.md) | The maths as implemented, and the choices made where the paper leaves room |
| [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) | Provenance and licence for every source |
| [`RUNBOOK.md`](RUNBOOK.md) | End-to-end reproduction, expected runtimes, troubleshooting |

## Citations

The method implemented here:

1. **Iida et al. (2024).** *A network-based trans-omics approach for predicting synergistic drug combinations (SyndrumNET).* Communications Medicine.

Methodological foundations:

2. **Vanunu et al. (2010).** *Associating genes and protein complexes with disease via network propagation.* PLoS Computational Biology. (PRINCE)
3. **Guney et al. (2016).** *Network-based in silico drug efficacy screening.* Nature Communications. (network proximity and the degree-preserving null model)
4. **Cowen et al. (2017).** *Network propagation: a universal amplifier of genetic associations.* Nature Reviews Genetics.

If you use this code, please cite the Iida et al. paper for the method and this repository for the implementation. [`CITATION.cff`](CITATION.cff) has the machine-readable form, and GitHub's "Cite this repository" button reads from it.

## Acknowledgements

The starting point for this project was **Rahim Rahimov's** Master's thesis, *Reproducing and Extending the SyndrumNET Model for Synergistic Drug Combination Prediction* (2025), which is what got me interested in the method in the first place. The implementation in this repository is my own.

## License

[MIT](LICENSE). The licence covers this implementation only; the integrated data sources carry their own terms.

## Author

Nihad Aslanzade, [github.com/nihadaslanzade-02](https://github.com/nihadaslanzade-02)
