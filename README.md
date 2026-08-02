# SyndrumNET: Network-Based Prediction of Synergistic Drug Combinations

[![CI](https://github.com/nihadaslanzade-02/SyndrumNET/actions/workflows/ci.yml/badge.svg)](https://github.com/nihadaslanzade-02/SyndrumNET/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-3776AB)](https://www.python.org/)
[![Ruff](https://img.shields.io/badge/lint-ruff-261230)](https://docs.astral.sh/ruff/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

A Python implementation of **SyndrumNET**, a network-based trans-omics method that predicts which pairs of drugs will act synergistically against a given disease, by asking where each drug's effect lands on the human interactome relative to where the disease lives.

The pipeline integrates molecular-interaction resources into one network, builds gene modules for diseases and drugs from expression signatures, propagates signal across the network with PRINCE, and scores every drug pair on three independent axes: **topology**, **proximity**, and **transcriptional reversal**.

The scoring layer is complete and tested. The data layer is not: three of the seven configured interaction sources have parsers, and disease modules are built from expression signatures alone. [`docs/PLACEHOLDERS.md`](docs/PLACEHOLDERS.md) is a stage-by-stage account of what is missing and what each omission changes about the results.

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

The intuition the whole method is built on: **two drugs help each other when they reach the disease through different neighbourhoods, and duplicate each other when they reach it through the same one.** [`scoring/tqab.py`](src/syndrumnet/scoring/tqab.py) implements the six-class scheme of Cheng et al. (2019), which partitions on how many of the two drugs sit closer to the disease than chance, and whether the drug modules share a neighbourhood:

| Class | `z_QA` | `z_QB` | `s_AB` | `T_QAB` |
|---|---|---|---|---|
| Overlapping Exposure | `< 0` | `< 0` | `< 0` | 0 |
| **Complementary Exposure** | `< 0` | `< 0` | `>= 0` | **2** |
| Indirect Exposure | one `< 0` | | `< 0` | 0 |
| Single Exposure | one `< 0` | | `>= 0` | 0 |
| Non-exposure | `>= 0` | `>= 0` | `< 0` | 0 |
| Independent Action | `>= 0` | `>= 0` | `>= 0` | 0 |

Only Complementary Exposure is associated with synergy, and the score is binary rather than graded. Because the final prediction is an unweighted sum of the three components, that constant of 2 is what sets the topological axis' weight against the other two.

### P: proximity

`P(Q,A,B) = (P_QA + P_QB)/2`, where each term is a disease-drug distance converted to a z-score against a null distribution:

```
z(d_QA) = (d_QA - mean(d_random)) / std(d_random)
```

The null is built by **degree-preserving randomisation**. Genes are binned into 20 degree strata and random modules are resampled within strata, so a "close" drug module cannot score well merely for containing hub genes. This is the standard construction from Guney et al. (2016).

### C: transcriptional reversal

`C(Q,A,B) = (C_QA + C_QB)/2`, comparing each drug's L1000 expression signature against the disease's CREEDS signature. Drug modules are the top 5% of genes by absolute fold change. The method calls for aggregating a compound's profiles across cell lines first; `parse_lincs` does not, and each cell line currently yields its own module ([`docs/PLACEHOLDERS.md`](docs/PLACEHOLDERS.md) section 4).

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
| Test suite | **103 passing**: unit coverage of every metric and scoring module, plus end-to-end runs over a synthetic network |
| Lint | `ruff` clean, enforced in CI |
| Scoring pipeline | `SynergyPredictor.predict_all()` runs end to end and is covered by integration tests |
| Interaction sources | 3 of the 7 in `configs/default.yaml` have a parser. The other 4 are dropped without a warning |
| Disease modules | Expression signatures only; the susceptibility-gene half is a placeholder that returns an empty set |
| Drug modules | Keyed by LINCS signature ID rather than compound, with no aggregation across cell lines. This blocks evaluation |
| ID harmonisation | One HTTP request per identifier, twice over. `batch_convert()` exists and is unused |
| Run on real data | Not yet performed. No committed predictions, figures or evaluation tables |

Each of those data-layer rows is worked through in
[`docs/PLACEHOLDERS.md`](docs/PLACEHOLDERS.md), including what it changes
about the scores rather than just that it is missing.

### What the tests found

The integration suite was written to exercise the wiring rather than the
components, and it earned its keep immediately. Each fix below is pinned by a
test that fails when the fix is reverted.

**Proximity depended on pair position.** `P_QA` is defined as a property of the
disease and one drug, but `compute_pqab` seeded the null model with `seed` for
the first drug and `seed + 1` for the second. A drug therefore scored
differently depending on which side of the pair it landed on, and since pairs
are enumerated in list order, that noise tracked each drug's position in the
input. Null-model seeds are now derived from the module's own contents, so a
module always draws the same null wherever it appears.

**Redundant pairs could never be detected.** `separation_score` computes
`s_AB = <d_AB> - (<d_AA> + <d_BB>)/2`, but the intra-module terms let every
gene match itself at distance 0, so `<d_AA>` and `<d_BB>` were identically 0.
That collapsed `s_AB` to `<d_AB>`, which is non-negative for any two distinct
modules, so the `s_AB <= 0` branch was unreachable: on a random sample of 200
module pairs, 200 classified as complementary and none as redundant. The
intra-module terms now measure the distance to the nearest *other* gene, and
the same sample splits 153 redundant to 47 complementary.

**Per-drug scores were recomputed once per pair.** Both `P_QA` and `C_QA`
depend on a single drug, yet the batch functions recomputed them for every pair
that drug appeared in. At the 1,488 drugs the default config declares, that is
2.2 million null-model sweeps where 1,488 suffice. Both batch paths now compute
one value per drug and reuse it, and the predictor computes the z-scores once
and hands the same values to both TQAB and PQAB.

**The topological score did not match the published definition.** This one came
out of reading the source papers rather than the tests. `T_QAB` was a graded
score built from three unpublished constants: a cutoff of 3 hops on the raw
disease-drug distance, and divisors of 10 and 5 in the class scores. The
published method classifies on the proximity **z-scores** (`z_QA < 0`, meaning
closer to the disease than a degree-matched random module) rather than on an
absolute hop count, distinguishes six classes rather than three, and awards a
**binary** `T_QAB = 2` for Complementary Exposure and 0 for everything else.
The implementation now follows Cheng et al.'s Figure 2 table, and the six
classes are covered by a parametrised test.

While reconciling it, one discrepancy in the sources is worth recording:
Iida et al. give the Class II condition inline as `s_AB < 0`, which contradicts
both their own description of Class II as "two separated drug modules" and
Cheng et al.'s panel P2, which reads `sAB >= 0`. The convention that
`s_AB < 0` means the modules share a neighbourhood is stated explicitly in
Cheng et al., so this implementation uses `s_AB >= 0` and documents the
conflict in [`docs/METHOD_NOTES.md`](docs/METHOD_NOTES.md).

**Two modules had no coverage at all, and both were wrong.** `similarity_layers.py`
and `viz/plots.py` were the last untested files. Covering them turned up eight
more defects, of which two were silent and would never have announced
themselves:

- `tanimoto_similarity_matrix` gave a wrong answer for boolean fingerprints. A
  bool array is the most natural way to write "binary fingerprint matrix", but
  numpy's `@` on booleans is logical rather than arithmetic, so every
  shared-bit count saturated at 1. On 2048-bit Morgan-style fingerprints every
  entry came back near 0.005 and the diagonal read 0.33 instead of 1.0. Integer
  and float inputs were fine, which is why nothing caught it.
- `plot_top_predictions` stacked all three score components from a single
  running offset. PQAB is an averaged proximity z-score, negative exactly when
  the drugs sit close to the disease, so its segment was drawn back across the
  positive TQAB segment. The bar's end point still landed on the right total,
  which is what makes it easy to miss: only the decomposition was wrong. It now
  stacks positive and negative parts from separate baselines, marks the total
  separately, and inverts the y axis so the best pair is at the top rather than
  the bottom.

Also: the log-log degree distribution used linear bins, and a bin starting at
degree 0 pulled the lower x-limit to 7e-18 on a graph with isolated nodes; a
save that failed left its figure open; and `visualization.dpi`,
`figure_format` and `top_k_predictions` were declared in the config and read by
nothing.

Earlier still, the package could not be imported at all: seven `typing` names
and `pandas` were used without being imported, and because Python evaluates
annotations at definition time, five of the eight subpackages raised
`NameError` on import. That took every entry-point script down with them and
left `pytest` collecting zero tests. Fixed in `650f52b`, along with the defect
the existing test suite had been written to catch but never got to run: the
degree-preserving null model called `_get_bin` with the wrong arity, and
separately handed it bin contents where it expected numeric bin edges.

CI guards that failure mode directly. Every module is imported in a dedicated
job, so a missing import fails the build even when no test happens to exercise
the affected function.

---

## Architecture

```
src/syndrumnet/
├── io/              downloaders, format parsers, HGNC/UniProt ID mapping
├── data/            network assembly from interaction sources; disease & drug modules
├── metrics/         shortest-path distances, separation, degree-preserving null models,
│                    transcriptional correlation
├── propagation/     PRINCE random walk with restart; similarity layers
├── scoring/         tqab.py, pqab.py, cqab.py, and the predictor that combines them
├── eval/            AUC-ROC / AUC-PR against known synergy resources, reporting
├── viz/             matplotlib figures
└── utils/           YAML config loading with CLI override, logging, seeding

scripts/             build_all_data · run_pipeline · evaluate · make_figures · gen_api_docs
configs/             default.yaml + one YAML per disease
tests/               config · distances · null_models · propagation · scoring · similarity_layers · plots · integration
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
pytest tests/ -v      # 103 tests, no data or network access needed
ruff check .          # same lint gate CI runs
```

The integration tests run the whole predictor over a small synthetic network
built in [`tests/conftest.py`](tests/conftest.py): four complete communities
chained together, with drugs placed so that "adjacent to the disease",
"two communities away" and "overlapping with another drug" are true by
construction rather than by whatever the code happened to output.

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

`visualization.dpi`, `visualization.figure_format` and `visualization.top_k_predictions` are read by `make_figures.py` and `evaluate.py`; every other key is consumed where you would expect.

[`configs/diseases/`](configs/diseases/) has one file per disease, but only `cml.yaml` has any content and nothing in the codebase reads the directory yet. `scripts/add_new_disease.py` and `scripts/add_new_drugs.py` are the intended extension points for adding a disease from a GEO series or drugs from LINCS; both are currently stubs.

---

## Data sources

| Layer | Sources |
|---|---|
| Molecular interactions | HuRI, CORUM, PhosphoSitePlus. Configured but not implemented: KEGG RPair, SignaLink, InnateDB, Instruct |
| Disease expression | CREEDS (79 diseases) |
| Drug expression | LINCS L1000 (1,488 compounds) |
| Disease genes | OMIM, ClinVar, GWAS, PheWAS, GWASdb, DisGeNET |
| Drug structures | PubChem, ChEMBL |
| ID mapping | HGNC, UniProt |

No data is committed to this repository; `data/` is populated by `build_all_data.py`. Licences vary and several sources are non-commercial or require registration, so read [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) before redistributing anything derived from them.

**Deviations from the published method** are listed in [`docs/METHOD_NOTES.md`](docs/METHOD_NOTES.md) and traced through to their effect on results in [`docs/PLACEHOLDERS.md`](docs/PLACEHOLDERS.md): four of the seven interaction sources are absent, disease modules carry no susceptibility genes, drug modules are keyed by signature rather than compound, and Morgan fingerprints stand in for KCF-S.

---

## Roadmap

Ordered by what blocks what, rather than by size. Full detail in [`docs/PLACEHOLDERS.md`](docs/PLACEHOLDERS.md).

1. **Batch the ID mapping.** `IDMapper.to_hgnc()` issues one HTTP request per identifier, and `network_builder.build()` runs it twice over the full gene list because it discards the result of the first call. At interactome scale that is tens of thousands of sequential round trips, which is the wall a first real run hits. `batch_convert()` already does the job in one request and nothing calls it.
2. **Join the LINCS metadata.** Until drug modules are keyed by compound instead of signature ID, `eval/` cannot match predictions to any known-synergy resource, so AUC-ROC and AUC-PR cannot be computed at all. The same change adds the missing aggregation across cell lines.
3. **Run the pipeline once on one disease**, with a reduced `n_randomizations`, to get a first result on real data and a committed prediction table.
4. **Load susceptibility genes.** The largest deviation from the method: a disease module is currently its expression signature alone. The parsers are easy, the disease-term mapping between CREEDS, ClinVar, DisGeNET and OMIM is not.
5. **Fill in the four missing interaction parsers** for KEGG RPair, SignaLink, InnateDB and Instruct, so the network is the seven-source integration the config already declares.
6. **Validate against the reference results**, disease by disease.

---

## Documentation

| Document | Contents |
|---|---|
| [`docs/API.md`](docs/API.md) | Every public class and function, generated from source by `scripts/gen_api_docs.py` |
| [`docs/METHOD_NOTES.md`](docs/METHOD_NOTES.md) | The maths as implemented, and the choices made where the paper leaves room |
| [`docs/PLACEHOLDERS.md`](docs/PLACEHOLDERS.md) | Every unimplemented stage, what it silently returns, and what it changes about the results |
| [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) | Provenance and licence for every source |
| [`RUNBOOK.md`](RUNBOOK.md) | End-to-end reproduction, expected runtimes, troubleshooting |

## Citations

The method implemented here:

1. **Iida, M., Kuniki, Y., Yagi, K., Goda, M., Namba, S., Takeshita, J.-I., Sawada, R., Iwata, M., Zamami, Y., Ishizawa, K. & Yamanishi, Y. (2024).** *A network-based trans-omics approach for predicting synergistic drug combinations.* Communications Medicine **4**, 154. [doi:10.1038/s43856-024-00571-2](https://doi.org/10.1038/s43856-024-00571-2)

Methodological foundations:

2. **Cheng, F., Kovács, I. A. & Barabási, A.-L. (2019).** *Network-based prediction of drug combinations.* Nature Communications **10**, 1197. [doi:10.1038/s41467-019-09186-x](https://doi.org/10.1038/s41467-019-09186-x) (the six topological classes and Complementary Exposure)
3. **Menche, J. et al. (2015).** *Uncovering disease-disease relationships through the incomplete interactome.* Science **347**, 1257601. (the separation measure `s_AB`)
4. **Vanunu, O. et al. (2010).** *Associating genes and protein complexes with disease via network propagation.* PLoS Computational Biology **6**, e1000641. (PRINCE)
5. **Guney, E. et al. (2016).** *Network-based in silico drug efficacy screening.* Nature Communications **7**, 10331. (network proximity and the degree-preserving null model)
6. **Cowen, L. et al. (2017).** *Network propagation: a universal amplifier of genetic associations.* Nature Reviews Genetics **18**, 551-562.

If you use this code, please cite the Iida et al. paper for the method and this repository for the implementation. [`CITATION.cff`](CITATION.cff) has the machine-readable form, and GitHub's "Cite this repository" button reads from it.

## Acknowledgements

The starting point for this project was **Rahim Rahimov's** Master's thesis, *Reproducing and Extending the SyndrumNET Model for Synergistic Drug Combination Prediction* (2025), which is what got me interested in the method in the first place. The implementation in this repository is my own.

## License

[MIT](LICENSE). The licence covers this implementation only; the integrated data sources carry their own terms.

## Author

Nihad Aslanzade, [github.com/nihadaslanzade-02](https://github.com/nihadaslanzade-02)
