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

1. **Drug-drug separation**: `s_AB = <d_AB> - (<d_AA> + <d_BB>)/2`
2. **Disease-drug proximity z-scores**: `z_QA`, `z_QB`

The six classes come from Cheng et al. (2019), Figure 2. They partition on how
many of the two drugs sit closer to the disease than chance, and whether the
two drug modules share a network neighbourhood:

| Class | Name | `z_QA` | `z_QB` | `s_AB` |
|---|---|---|---|---|
| P1 | Overlapping Exposure | `< 0` | `< 0` | `< 0` |
| **P2** | **Complementary Exposure** | `< 0` | `< 0` | `>= 0` |
| P3 | Indirect Exposure | one `< 0` | | `< 0` |
| P4 | Single Exposure | one `< 0` | | `>= 0` |
| P5 | Non-exposure | `>= 0` | `>= 0` | `< 0` |
| P6 | Independent Action | `>= 0` | `>= 0` | `>= 0` |

Only Complementary Exposure is associated with synergy: both drugs reach the
disease, but through different neighbourhoods, so their effects add rather
than duplicate. The score is binary.

```
T_QAB = 2  for Complementary Exposure
        0  for every other class
```

**Sign convention on `s_AB`**, inherited from Menche et al. (2015) via Cheng
et al. (2019), quoting the latter directly:

> For sAB < 0, the targets of the two drugs are located in the same network
> neighborhood, while for sAB >= 0, the two drug targets are topologically
> separated.

Note that Iida et al. (2024) writes the Class II condition inline as
`s_AB < 0`. That contradicts its own description of Class II as "two separated
drug modules", and contradicts Cheng et al.'s Figure 2 panel P2, which reads
`ZDA < 0, ZDB < 0, sAB >= 0`. This implementation follows the latter.

**Intra-module distances exclude self-comparisons.** `<d_AA>` is the mean, over
each gene in module A, of the distance to the nearest *other* gene in A. Letting
a gene match itself makes `<d_AA>` identically 0, which collapses `s_AB` to
`<d_AB>` and leaves it non-negative for any two distinct modules, so the
overlapping classes become unreachable. The cross-module term `<d_AB>` does
count shared genes at distance 0; Cheng et al. state this explicitly: "In
<dAB>, targets associated with both drugs A and B have a zero distance by
definition."

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
