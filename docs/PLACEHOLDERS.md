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

| # | Gap | Where | Visible at runtime? |
|---|---|---|---|
| 1 | Interactome uses 3 of the 7 configured sources | `io/downloaders.py`, `io/parsers.py`, `data/network_builder.py` | No |
| 2 | Disease modules carry no susceptibility genes | `data/modules.py` | Only at `DEBUG` |
| 3 | Drug keys are LINCS signature IDs, not compounds | `io/parsers.py` | No |
| 4 | No per-gene aggregation across cell lines | `io/parsers.py` | No |
| 5 | Morgan fingerprints substitute for KCF-S | `propagation/similarity_layers.py` | Only if RDKit is absent |
| 6 | `build()` harmonises the gene list twice | `data/network_builder.py` | No, marked `TODO` |
| 7 | ID mapping is one HTTP request per identifier | `io/id_mapping.py` | As runtime |
| 8 | Extension scripts are stubs | `scripts/add_new_*.py` | Yes, they print it |

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

## 3. Drug keys are LINCS signature IDs, not compounds

`parse_lincs` loads the metadata table and never joins it:

```python
meta = pd.read_csv(meta_filepath, sep='\t')  # noqa: F841

for drug in df.columns:
    ...
```

The metadata carries the signature-to-compound mapping, so `df.columns` are
raw signature identifiers. Everything downstream inherits those as drug names.

**Effect.** A compound profiled in several cell lines becomes several distinct
"drugs". Three consequences, in increasing order of seriousness:

1. Pair enumeration is quadratic in the wrong count. Where the config declares
   1,488 compounds, the signature matrix has more columns than that, and the
   pair count grows with the square of the difference.
2. Same-compound pairs are generated and scored. Two signatures of one drug
   have near-identical modules, so `s_AB` comes out negative and they are
   classified Overlapping Exposure. They score 0 and do no harm to the
   ranking, but they consume the compute.
3. The output cannot be evaluated. `eval/` matches predictions against known
   synergy resources, which are keyed by compound name. Signature IDs join to
   nothing, so AUC-ROC and AUC-PR cannot be computed at all until this is
   fixed.

**To finish.** Join `meta` on the signature ID column, group the signature
matrix by compound, and use the compound name as the module key. Point 4 below
is the same fix.

## 4. No aggregation across cell lines

Related to the above and worth stating separately, because the documentation
claimed otherwise until this commit. `docs/METHOD_NOTES.md` described a
per-gene median across cell lines as a deliberate implementation choice. There
is no such aggregation anywhere in the data path: `parse_lincs` takes each
column of the signature matrix on its own, and no `groupby` or `median` over
cell lines exists in `io/`, `data/` or `metrics/`.

The `median` in `metrics/transcription.py` is an unrelated score aggregator
and does not touch cell lines.

**Effect.** Whatever noise a single cell line carries goes straight into the
drug module. The top 5% of genes by fold change in one cell line is a noisier
set than the top 5% of a median profile, and the drug module is the input to
both the proximity and the transcriptional axis.

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

## 6. `build()` harmonises the gene list twice

Marked `TODO` in place:

```python
harmonized = self.id_mapper.harmonize_gene_list(all_genes)  # noqa: F841
gene_map = dict(zip(all_genes, self.id_mapper.to_hgnc(all_genes)))
```

The first result is discarded and the mapping is rebuilt from a second call.
Both go through `to_hgnc`, so they cannot disagree about the mapping itself;
`harmonize_gene_list` additionally deduplicates and drops unmapped
identifiers, and that filtered list is what is thrown away.

**Effect.** Correctness is unaffected, cost is not: this doubles the most
expensive step in the data build, for the reason in point 7. The `noqa` and
the `TODO` are there so it fails review rather than lint.

## 7. ID mapping issues one HTTP request per identifier

Not a placeholder, but the reason the data build has not been run.

`IDMapper.to_hgnc` loops over the identifiers and calls `self.mg.query(...)`
once per gene. Identifiers already present in a loaded HGNC table take a local
shortcut, so the cost falls on everything that is not already an HGNC symbol.
`IDMapper.batch_convert` does exactly this job through `mygene`'s `querymany`,
in one request, and nothing in the codebase calls it.

**Effect.** On an interactome of the expected size, the unmapped fraction of
15,000 to 20,000 identifiers becomes that many sequential HTTP round trips,
twice over because of point 6. This is the wall a first real run hits.

**To finish.** Route `to_hgnc` through `batch_convert`, and persist the cache
to `cache_dir`, which is created in `__init__` and then only used in memory.

## 8. Extension scripts are stubs

`scripts/add_new_disease.py` and `scripts/add_new_drugs.py` parse their
arguments, print what they would do, and exit. They are honest about it on
stdout, so they are listed here only for completeness.

---

## What this adds up to

Points 1 and 2 change the answer: a smaller network and a differently
constituted disease module mean the scores this pipeline would produce are not
the scores the paper describes, even though every formula between them is
implemented and tested. Point 3 blocks evaluation outright. Point 7 blocks the
first run.

None of this affects the scoring layer, which is tested end to end over a
synthetic network in `tests/`. The gap is entirely in the data path between a
real download and a module.
