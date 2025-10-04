\# SyndrumNET Implementation Notes



\## Mathematical Details



\### Network Proximity



Disease-drug proximity is computed as:



```

d(S,T) = (1/|S|) \* sum\_{s in S} min\_{t in T} dist(s,t)

```



Z-score normalization uses degree-preserving randomization:



```

z(d\_QA) = (d\_QA - mean(d\_random)) / std(d\_random)

```



\### Topological Class Score (TQAB)



Drug pair topology relative to disease is classified based on:



1\. \*\*Drug-drug separation\*\*: s\_AB = d\_AB - (d\_AA + d\_BB)/2

2\. \*\*Disease-drug proximities\*\*: d\_AQ, d\_BQ



Complementary: s\_AB > 0 and both drugs close to disease

Redundant: s\_AB < 0 or drugs far from disease



\### PRINCE Propagation



Iterative update:

```

F^(t+1) = α \* W \* F^(t) + (1-α) \* F^(0)

```



where:

\- F^(t): node scores at iteration t

\- W: normalized adjacency matrix

\- α: restart probability (default: 0.5)

\- F^(0): initial seed vector



\### Final Prediction



```

Score\_Q,AB = TQAB + PQAB + CQAB

```



where:

\- TQAB: topological class score

\- PQAB = (P\_QA + P\_QB)/2: average proximity

\- CQAB = (C\_QA + C\_QB)/2: average transcriptional correlation



\## Implementation Choices



\### Cell Line Handling (L1000)



Drug signatures are aggregated across cell lines using:

\- Median fold-change per gene across all cell lines

\- Top 5% genes by |fold-change| define drug modules



\### Null Model



Degree-preserving randomization:

\- 1000 randomizations by default

\- Genes binned by degree (20 bins)

\- Random sampling within bins



\### Distance Handling



Disconnected nodes:

\- Set to large finite value (1000.0) rather than infinity

\- Allows graceful handling of network components



\## Deviations from Paper



1\. \*\*KCF-S fingerprints\*\*: Using Morgan fingerprints as proxy (RDKit)

2\. \*\*Some data sources\*\*: Placeholders where registration/API required



All deviations documented in code comments.

