# SyndrumNET Implementation Notes

Mathematical details and the implementation choices made where the published
method leaves room for interpretation.

## Mathematical details

### Network proximity

Disease-drug proximity is the average, over the source module, of the distance
to the nearest gene in the target module:

```
d(S,T) = (1/|S|) * sum_{s in S} min_{t in T} dist(s,t)
```

Raw distances are not comparable across modules of different size and degree
composition, so each is converted to a z-score against a null distribution:

```
z(d_QA) = (d_QA - mean(d_random)) / std(d_random)
```

### Topological class score (TQAB)

A drug pair is classified relative to the disease module using two quantities:

1. **Drug-drug separation**: `s_AB = d_AB - (d_AA + d_BB)/2`
2. **Disease-drug proximities**: `d_AQ`, `d_BQ`

| Condition | Class |
|---|---|
| `s_AB > 0` and both drugs close to the disease | Complementary |
| `s_AB > 0`, at least one drug far from the disease | Intermediate |
| `s_AB <= 0` | Redundant |

The threshold for "close" is a distance of 3 hops, set in `scoring/tqab.py`.

### PRINCE propagation

Random walk with restart, iterated to convergence:

```
F^(t+1) = alpha * W * F^(t) + (1 - alpha) * F^(0)
```

where `F^(t)` is the score vector at iteration `t`, `W` the normalised adjacency
matrix, `alpha` the restart probability (default 0.5), and `F^(0)` the seed
vector. Iteration stops when the L2 norm of the update falls below 1e-6, or at
1000 iterations.

### Final prediction

```
Score_Q,AB = TQAB + PQAB + CQAB
```

where `PQAB = (P_QA + P_QB)/2` is the average proximity z-score and
`CQAB = (C_QA + C_QB)/2` the average transcriptional correlation.

## Implementation choices

### Cell line handling (L1000)

A compound is profiled in several cell lines, and the method does not specify
how to reconcile them. Here:

- Fold changes are aggregated per gene by **median across all cell lines**,
  which is robust to a single anomalous line.
- The drug module is the **top 5% of genes by absolute fold change**.

### Null model

Degree-preserving randomisation, so that a module cannot score as "proximal"
merely because it contains hub genes:

- 1000 randomisations by default
- Genes binned into 20 degree strata
- Random modules resampled within strata

### Distance handling

Disconnected node pairs are assigned a large finite value (1000.0) rather than
infinity. An `inf` anywhere in a module would make every average containing it
`inf`, discarding the information carried by the connected pairs; a large finite
sentinel keeps the average dominated by real distances while still penalising
disconnection heavily.

## Deviations from the paper

1. **KCF-S fingerprints** - Morgan fingerprints (RDKit) are used as a proxy.
2. **Some data sources** - placeholders where registration or an API key is
   required (KEGG RPair, several disease-gene resources).

All deviations are also marked in the code comments at the point where they
apply.
