# Placeholders: what is not implemented, and what it costs

Several stages of this pipeline are placeholders. Each one returns an empty
result rather than raising, which is deliberate so that the pipeline can be
developed stage by stage, but it also means the gap between the method as
published and the method as implemented is invisible at runtime. This file is
that gap, traced from the code to the effect on results.

Every claim below was checked against the code rather than inferred from the
docstrings, which in two places said more than the implementation delivers.

**Nothing here has been run on real data.** These are gaps in the code path,
not measured error.

## Summary

| # | Gap | Where | Status |
|---|---|---|---|
| 1 | Interactome uses 3 of the 7 configured sources | `io/downloaders.py`, `io/parsers.py`, `data/network_builder.py` | Open, silent |
| 2 | Disease modules carry no susceptibility genes | `data/modules.py` | Open, `DEBUG` only |
| 3 | Drug keys were LINCS signature IDs, not compounds | `io/parsers.py` | **Fixed** |
| 4 | No per-gene aggregation across cell lines | `io/parsers.py` | **Fixed** |
| 5 | Morgan fingerprints substitute for KCF-S | `propagation/similarity_layers.py` | Open, by design |
| 6 | `build()` harmonised the gene list twice | `data/network_builder.py` | **Fixed** |
| 7 | ID mapping was one HTTP request per identifier | `io/id_mapping.py` | **Fixed** |
| 8 | Extension scripts are stubs | `scripts/add_new_*.py` | Open, self-declared |
| 9 | Five of six per-disease configs are empty files | `configs/diseases/` | Open, now loud |

Items 3, 4, 6 and 7 were the ones blocking a first run and its evaluation.
They are kept below with what they were and what replaced them, because the
reasoning is the same reasoning anyone finishing items 1 and 2 will need.

---

## 1. The interactome is built from three of the seven configured sources

`configs/default.yaml` declares seven interaction sources. Only three survive
to the network:

| Source | URL in `DataDownloader.URLS` | `download_*` method | `parse_*` function | Registered in `add_source` |
|---|---|---|---|---|
| HuRI | yes | yes | yes | yes |
| CORUM | yes | yes | yes | yes |
| PhosphoSitePlus | yes | yes | yes | yes |
| KEGG RPair | yes | no | returns an empty frame | no |
| SignaLink | yes | no | no | no |
| InnateDB | yes | no | no | no |
| Instruct | yes | no | no | no |

`NetworkBuilder.add_source` does log `No parser for <source>, skipping`, but
in the shipped pipeline that line is unreachable. `build_all_data.py` guards
the call with `if source in files`, and `files` comes from
`DataDownloader.download_all()`, which fetches only the first three. The four
missing sources are dropped before `add_source` is ever called, so a full run
produces no warning, no error and no record that four sevenths of the
declared network is absent.

**Effect.** Every quantity in the method is a function of shortest-path
distance, so a sparser interactome moves all of them. Missing edges lengthen
paths, which inflates `d(S,T)` and `<d_AB>`, and disconnected pairs fall back
to the 1000.0 sentinel. The proximity z-scores are not simply shifted by
this, because the degree-preserving null is drawn from the same reduced graph:
both the observed distance and its null distribution change, in the same
direction but not by the same amount. Absolute scores from this build are
therefore not comparable to the paper's, and the sign of `s_AB`, which decides
the topological class outright, can flip for module pairs that a denser
network would have connected.

**To finish.** KEGG RPair needs the REST endpoint already listed in `URLS`
(`https://rest.kegg.jp/get/rpair`) wired to a downloader method, and
`parse_kegg_rpair` needs to convert reactant pairs into gene-level edges via
the enzymes that catalyse them, which is the part that is not mechanical.
SignaLink, InnateDB and Instruct each need a parser, a downloader method, and
an entry in the `parsers` dict inside `add_source`. InnateDB ships PSI-MITAB,
so that one is close to standard.

## 2. Disease modules carry no susceptibility genes

`ModuleBuilder._load_susceptibility_genes` returns an empty set:

```python
genes = set()
# TODO: Implement parsers for each source
# For now, return empty set
return genes
```

Its caller does `module_genes |= (susc_genes & self.network_genes)`, which
against an empty set is a no-op. The only trace is a `logger.debug` line, so
at the default log level a caller who passes `susceptibility_files` gets no
indication that the argument did nothing.

In practice the argument is never passed anyway: `build_all_data.py` calls
`build_disease_modules(raw_dir / 'creeds_disease.txt')` with no
`susceptibility_files`, so the branch is not even entered. Two independent
reasons for the same outcome.

**Effect.** This is the largest deviation in the pipeline. The method defines
a disease module as the union of two different kinds of evidence:
susceptibility genes (OMIM, ClinVar, GWAS, PheWAS, GWASdb, DisGeNET) and
genes dysregulated in the disease (CREEDS). As implemented, a disease module
is the second half only, intersected with the network.

Those two halves do not sit in the same place on the interactome. Differential
expression picks up downstream response, much of it shared across diseases and
concentrated in high-degree, well-annotated genes; susceptibility genes are
causal and typically sit further upstream. Dropping the causal half moves
where "the disease" is, and since every proximity z-score measures distance to
that module, it changes what the drugs are being scored against. Predictions
would remain internally consistent, but they would answer a different question
from the one the paper asks.

Worth noting: `DataDownloader.download_disease_genes()` does fetch ClinVar and
DisGeNET, so those files land in `data/raw/` on every build and are never
read. OMIM needs an API key and is only warned about.

**To finish.** The parsers themselves are easy. The real work is the disease
vocabulary: CREEDS names diseases in free text, ClinVar uses MedGen concepts,
DisGeNET uses UMLS CUIs, and OMIM uses its own numbering. Joining them needs a
disease-term mapping, and how that mapping is drawn will change module
membership more than the choice of source will.

## 3 and 4. Drug keys, and aggregation across cell lines: FIXED

One gap and one fix, so they are written together.

**What it was.** `parse_lincs` loaded the metadata table and never joined it:

```python
meta = pd.read_csv(meta_filepath, sep='\t')  # noqa: F841

for drug in df.columns:
    ...
```

The metadata carries the signature-to-compound mapping, so `df.columns` were
raw signature identifiers, and everything downstream inherited those as drug
names. Separately, `docs/METHOD_NOTES.md` described a per-gene median across
cell lines as a deliberate implementation choice; no such aggregation existed
anywhere in the data path.

**What it cost.** A compound assayed in several cell lines became several
distinct "drugs". Three consequences, in increasing order of seriousness:

1. Pair enumeration was quadratic in the wrong count.
2. Same-compound pairs were generated and scored. Two signatures of one drug
   have near-identical modules, so `s_AB` comes out negative and they classify
   as Overlapping Exposure. They score 0 and do no harm to the ranking, but
   they consume the compute.
3. The output could not be evaluated at all. `eval/benchmarks.py` matches
   predictions against synergy resources keyed by compound name. Signature IDs
   join to nothing, so AUC-ROC and AUC-PR were unreachable.

**What it is now.** `parse_lincs` joins the metadata, groups the signature
matrix by compound, and aggregates per gene before taking the percentiles.
Keys are compound names. Aggregation is `median` by default and `mean` on
request; median because a single anomalous cell line should not be able to
carry a gene into the module on its own, which
`test_median_resists_an_outlier_line_where_mean_does_not` demonstrates
directly.

Three details worth knowing:

- **The metadata schema is detected, not assumed.** It is not stable across
  LINCS and L1000CDS2 releases, so the signature-ID and compound-name columns
  are found from the candidate lists `SIGNATURE_ID_COLUMNS` and
  `COMPOUND_NAME_COLUMNS`. A file matching neither fails with a `KeyError`
  naming both what was tried and what the file contains. **This has not been
  run against a real download**; the candidate lists cover the column names
  these releases have used, and the failure mode if they are wrong is a loud
  error rather than a wrong answer.
- **`pert_id` is deliberately the last name candidate.** A BRD identifier
  joins to no synergy resource, so if detection falls through to it, expect
  evaluation to find nothing.
- **Signatures with no metadata row are dropped with a warning.** Keeping them
  under their raw ID would silently reintroduce exactly the problem the join
  removes. A join that matches nothing at all raises instead of returning an
  empty result.

## 5. Morgan fingerprints substitute for KCF-S

`kcf_fingerprint_similarity` computes an RDKit Morgan fingerprint (radius 2,
2048 bits) and calls it a KCF-S proxy. KCF-S encodes KEGG-defined chemical
functional substructures; Morgan encodes circular atom environments. Both
produce Tanimoto-comparable bit vectors, which is why the substitution is
mechanically possible, but they are not measuring the same notion of chemical
similarity.

**Effect on results today: none.** `propagation/similarity_layers.py` is
exported from the package but no scoring path calls it. The similarity layers
are infrastructure for using chemical and disease similarity as propagation
priors, which is not wired up. It would matter the moment they are.

**One sharp edge to know about.** If RDKit is not installed the function logs
a warning and returns `0.0`, which is indistinguishable from "these molecules
are genuinely dissimilar". A missing dependency currently looks like data.

## 5. Morgan fingerprints substitute for KCF-S

`kcf_fingerprint_similarity` computes an RDKit Morgan fingerprint (radius 2,
2048 bits) and calls it a KCF-S proxy. KCF-S encodes KEGG-defined chemical
functional substructures; Morgan encodes circular atom environments. Both
produce Tanimoto-comparable bit vectors, which is why the substitution is
mechanically possible, but they are not measuring the same notion of chemical
similarity.

**Effect on results today: none.** `propagation/similarity_layers.py` is
exported from the package but no scoring path calls it. The similarity layers
are infrastructure for using chemical and disease similarity as propagation
priors, which is not wired up. It would matter the moment they are.

**One sharp edge to know about.** If RDKit is not installed the function logs
a warning and returns `0.0`, which is indistinguishable from "these molecules
are genuinely dissimilar". A missing dependency currently looks like data.

## 6. `build()` harmonised the gene list twice: FIXED

**What it was**, marked `TODO` in place:

```python
harmonized = self.id_mapper.harmonize_gene_list(all_genes)  # noqa: F841
gene_map = dict(zip(all_genes, self.id_mapper.to_hgnc(all_genes)))
```

The first result was discarded and the mapping rebuilt from a second call.
Both go through `to_hgnc`, so they could not disagree about the mapping
itself. The question the `TODO` posed was which one was authoritative, and the
answer is `to_hgnc`: `harmonize_gene_list` deduplicates and drops unmapped
identifiers, which loses exactly the positional correspondence to the original
names that the `zip` on the next line depends on. It was never the right tool
there.

**What it is now.** The discarded call is gone, with a comment recording why
rather than a `TODO` asking the question again.

## 7. ID mapping issued one HTTP request per identifier: FIXED

Not a placeholder, but the reason the data build had never been run.

**What it was.** `IDMapper.to_hgnc` looped over the identifiers and called
`self.mg.query(...)` once per gene. `IDMapper.batch_convert` did exactly this
job through `mygene`'s `querymany`, in one request, and nothing called it. On
an interactome of the expected size, the non-symbol fraction of 15,000 to
20,000 identifiers became that many sequential HTTP round trips, twice over
because of point 6.

**What it is now.** One `querymany` per call, over the deduplicated set of
identifiers not already known locally. Identifiers already in a loaded HGNC
table still never leave the machine. Beyond the batching:

- **Results are cached to disk** under `cache_dir`, which was created in
  `__init__` and then only ever used in memory. The second run of the data
  build pays nothing for ID mapping. A corrupt cache is logged and ignored,
  since it is a performance problem and not a correctness one.
- **Misses are cached too**, so a genuinely unknown identifier is resolved
  once rather than re-queried every run.
- **A failed batch is not cached.** Recording a network error as "not found"
  would make every retry return `None` without going back to the server, which
  is the kind of silent-wrong this repository has had enough of.
- **`from_type` is validated.** Every value other than `'auto'` used to be
  accepted and then ignored, because the old query never passed a scope at
  all; a typo silently searched everything.
- **The local symbol shortcut respects `from_type`.** `'TP53'` is an HGNC
  symbol, but under `from_type='entrez'` it is not an Entrez ID and must not
  short-circuit to itself.
- `harmonize_gene_list`'s unmapped count was `len(input) - len(output)`, so
  deduplicating two copies of a gene that mapped perfectly well was reported
  as a mapping failure. It now counts the failures themselves.

`to_entrez` is batched the same way. `batch_convert` is kept as the raw
`mygene`-scope escape hatch.

## 8. Extension scripts are stubs

`scripts/add_new_disease.py` and `scripts/add_new_drugs.py` parse their
arguments, print what they would do, and exit. They are honest about it on
stdout, so they are listed here only for completeness.

## 9. Five of the six per-disease configs are empty files

`configs/diseases/` holds one file per disease in `config.diseases`. Only
`cml.yaml` has content, 262 bytes of disease name, key genes and known drugs.
`aml.yaml`, `asthma.yaml`, `colorectal_cancer.yaml`, `diabetes_t2.yaml` and
`hypertension.yaml` are zero bytes.

Nothing in the codebase reads the directory, so this breaks no run today. It
did make `load_config` fail badly, though: `yaml.safe_load` returns `None`
for an empty file, which reached `Config.__init__` and raised
`'NoneType' object has no attribute 'items'` without naming the file. That is
now an explicit `ValueError`, covered in `tests/test_config.py`.

The files are left empty on purpose. Filling them means naming each disease's
key genes and known drugs, which are domain facts, and a plausible-looking
guess in a config file is worse than an obviously empty one.

## Not a placeholder, but adjacent: config keys nothing read

Three keys under `visualization` in `configs/default.yaml` were declared and
consumed by nothing: `dpi`, `figure_format` and `top_k_predictions`. The plot
functions hardcoded 300 dpi and `.png`, and `make_figures.py` hardcoded
`k=20`. All three are now threaded through, with a test asserting the block
is still read so they cannot quietly drift back out of use.

---

## What this adds up to

The two that blocked anything from happening are closed. Point 7 was the wall
a first run hit, and point 3 was why its output could not have been evaluated
even if it had finished.

What remains changes the answer rather than preventing one. Points 1 and 2 mean
a smaller network and a differently constituted disease module, so the scores
this pipeline would produce are not the scores the paper describes, even though
every formula between them is implemented and tested. Point 1 is mechanical
work; point 2 is a genuine research decision about disease vocabularies.

None of this affects the scoring layer, which is tested end to end over a
synthetic network in `tests/`. The gap is entirely in the data path between a
real download and a module.

**Still true: nothing here has been run on real data.** The fixes above are
covered by tests over synthetic inputs, and every one of them was verified by
reverting it alone and confirming the matching test fails. That is not the same
as a successful download, and the LINCS column detection in particular has
never met a real metadata file.
