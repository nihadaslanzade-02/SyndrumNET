# SyndrumNET — Network-Based Prediction of Synergistic Drug Combinations

[![Python](https://img.shields.io/badge/python-3.10%2B-3776AB)](https://www.python.org/)
[![Status](https://img.shields.io/badge/status-work%20in%20progress-yellow)](#status--what-runs-today)
[![Tests](https://img.shields.io/badge/tests-9-informational)](tests/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

A Python implementation of **SyndrumNET**, a network-based trans-omics method that predicts which pairs of drugs will act synergistically against a given disease — by asking where each drug's effect lands on the human interactome relative to where the disease lives.

The pipeline integrates seven molecular-interaction resources into one network, builds gene modules for diseases and drugs from expression signatures, propagates signal across the network with PRINCE, and scores every drug pair on three independent axes: **topology**, **proximity**, and **transcriptional reversal**.

---

## Status — what runs today

This repository is a **structurally complete implementation that has not yet been run end to end.** Being specific, because the difference matters if you are about to clone it:

| Layer | State |
|---|---|
| Architecture, module layout, config system, docstrings | Complete — ~4,800 lines across 8 subpackages |
| Test suite | 9 tests written, covering distances, null models, propagation and scoring |
| **Package import** | **Broken.** 8 missing names across 7 files (7 `typing` imports, 1 `pandas`) mean 5 of the 8 subpackages raise `NameError` on import |
| Test run | `pytest` currently collects **0 tests**. With the missing imports added, **8 of 9 pass** |
| Null model | The 1 remaining failure is real: `_get_bin` is called with 3 arguments against a 2-argument signature (`null_models.py:121`), and separately is handed the bin *contents* where it expects numeric bin *edges* (`null_models.py:69`) |
| Data acquisition | Downloaders and parsers written; KEGG RPair and several disease-gene sources are explicit placeholders |
| Committed results | None — no predictions, figures or evaluation tables have been produced yet |
| `docs/API.md`, `conda-lock.yml` | Present but empty |

The import failures are all one-line fixes. The null-model bug is the interesting one: it sits in the degree-preserving randomisation that turns raw network proximity into a z-score, so it would silently invalidate every PQAB value. `tests/test_null_models.py` catches it immediately — the test was right all along, it just never got the chance to run, because the package it imports could not be loaded.

Fixing the imports first, then the null model, is the shortest path to a pipeline that can actually be executed. See [Roadmap](#roadmap).

---

## The method

SyndrumNET scores a drug pair **(A, B)** against a disease **Q** by summing three components:

```
Score(Q,A,B)  =  T(Q,A,B)  +  P(Q,A,B)  +  C(Q,A,B)
                 topology     proximity    transcription
```

### Network distance and separation

All three components rest on one primitive — the average shortest-path distance from one gene module to another:

```
d(S,T) = (1/|S|) · Σ_{s∈S}  min_{t∈T} dist(s,t)
```

From it comes the **separation score**, which asks whether two modules occupy distinct network neighbourhoods or overlap:

```
s_AB = ⟨d_AB⟩ − (⟨d_AA⟩ + ⟨d_BB⟩)/2
```

Positive `s_AB` means the modules are separated; negative means they overlap. Disconnected node pairs contribute a large finite value (1000.0) rather than infinity, so a fragmented network degrades gracefully instead of poisoning every average with `inf`.

### T — topological class

The intuition the whole method is built on: **two drugs help each other when they hit different parts of the disease neighbourhood, and duplicate each other when they hit the same part.** [`scoring/tqab.py`](src/syndrumnet/scoring/tqab.py) turns that into a three-way classification:

| Condition | Class | Score |
|---|---|---|
| `s_AB > 0` and both drugs within distance 3 of the disease | **Complementary** | `1 − d̄_Q/10` |
| `s_AB > 0`, at least one drug further away | **Intermediate** | `0.5 − d̄_Q/10` |
| `s_AB ≤ 0` (drug modules overlap) | **Redundant** | `−|s_AB|/5` |

### P — proximity

`P(Q,A,B) = (P_QA + P_QB)/2`, where each term is a disease-drug distance converted to a z-score against a null distribution:

```
z(d_QA) = (d_QA − μ_random) / σ_random
```

The null is built by **degree-preserving randomisation** — genes are binned into 20 degree strata and random modules are resampled within strata, so a "close" drug module cannot score well merely for containing hub genes. This is the standard construction from Guney et al. (2016), and it is the piece currently broken (see [Status](#status--what-runs-today)).

### C — transcriptional reversal

`C(Q,A,B) = (C_QA + C_QB)/2`, comparing each drug's L1000 expression signature against the disease's CREEDS signature. Drug modules are the top 5% of genes by absolute fold change, median-aggregated across cell lines.

### PRINCE propagation

[`propagation/prince.py`](src/syndrumnet/propagation/prince.py) implements random walk with restart:

```
F^(t+1) = α · W · F^(t) + (1−α) · F^(0)
```

with `α = 0.5` by default, iterating to a 1e-6 tolerance. The adjacency matrix is normalised column-wise (row-wise and symmetric are also selectable), kept sparse in CSR form, and built once per network so repeated propagations reuse it.

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

scripts/             build_all_data · run_pipeline · evaluate · make_figures
configs/             default.yaml + one YAML per disease
tests/               distances · null_models · propagation · scoring
```

Design decisions worth noting: `src/` layout so tests run against the installed package rather than the working directory; every public function carries a full NumPy-style docstring with parameters, returns and references; a single YAML config drives every stage with CLI overrides layered on top; seeding is centralised in `utils/seeds.py` rather than scattered across call sites; and scoring is split into one module per component so the three axes stay independently testable.

---

## Getting started

```bash
conda env create -f environment.yml
conda activate syndrumnet
pip install -e .
```

Requires Python 3.10+, roughly 16 GB RAM and ~50 GB of disk for the full data build.

**Before the pipeline can run, the import errors listed under [Status](#status--what-runs-today) need fixing.** Once they are, the intended flow is:

```bash
python scripts/build_all_data.py --config configs/default.yaml   # downloads + network build
python scripts/run_pipeline.py   --config configs/default.yaml   # TQAB + PQAB + CQAB
python scripts/evaluate.py       --config configs/default.yaml   # AUC-ROC / AUC-PR
python scripts/make_figures.py   --config configs/default.yaml   # figures
```

or `make all`. Results land in `reports/tables/predictions_<disease>.csv` and `reports/figures/`.

Run the tests with:

```bash
pytest tests/ -v
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

No data is committed to this repository — `data/` is populated by `build_all_data.py`. Licences vary and several sources are non-commercial or require registration; see [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) before redistributing anything derived from them.

**Deviations from the published method**, both documented in [`docs/METHOD_NOTES.md`](docs/METHOD_NOTES.md): Morgan fingerprints (RDKit) stand in for KCF-S fingerprints, and sources requiring registration or API keys are placeholders.

---

## Roadmap

1. **Add the 8 missing imports** — 7 `typing` names across 6 files, plus `pandas` in `eval/metrics.py`. This alone takes the test suite from 0 collected to 8 passing.
2. **Fix `_get_bin` in `metrics/null_models.py`** — correct the call arity at line 121, and give `_build_degree_bins` a way to return its bin edges so line 69 has something numeric to pass. Without this the proximity z-scores cannot be computed at all.
3. **Run the pipeline once, on one disease**, with a reduced `n_randomizations`, to get a first end-to-end result and a committed prediction table.
4. **Fill in the placeholder parsers** — KEGG RPair and the disease-gene sources — so the network is the full seven-source integration the method specifies.
5. **Generate `docs/API.md` and `conda-lock.yml`**, both currently empty while the docs reference them.
6. **Validate against the reference results** in the source thesis, disease by disease.

---

## Citations

This implementation follows:

1. **Rahim Rahimov (2025).** *Reproducing and Extending the SyndrumNET Model for Synergistic Drug Combination Prediction.* Master's Thesis.
2. **Iida et al. (2024).** *A network-based trans-omics approach for predicting synergistic drug combinations (SyndrumNET).* Communications Medicine.

Methodological foundations:

3. **Vanunu et al. (2010).** *Associating genes and protein complexes with disease via network propagation.* PLoS Computational Biology. — PRINCE
4. **Guney et al. (2016).** *Network-based in silico drug efficacy screening.* Nature Communications. — network proximity and the degree-preserving null model
5. **Cowen et al. (2017).** *Network propagation: a universal amplifier of genetic associations.* Nature Reviews Genetics.

## License

[MIT](LICENSE). The licence covers this implementation only — the integrated data sources carry their own terms.

## Contact

Nihad Aslanzade — [github.com/nihadaslanzade-02](https://github.com/nihadaslanzade-02)
