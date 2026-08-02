# Data Sources and Licenses

No data is committed to this repository. Everything below is fetched by
`scripts/build_all_data.py` into `data/raw/`.

Licences vary, and several sources are non-commercial or require registration.
Check the individual terms before redistributing anything derived from them.

## Molecular interactions

| Source | Description | License | Citation |
|---|---|---|---|
| [HuRI](http://www.interactome-atlas.org/) | Human Reference Interactome | Academic use permitted | Luck et al. (2020) *Nature* |
| [CORUM](https://mips.helmholtz-muenchen.de/corum/) | Protein complexes | Free for academic use | Giurgiu et al. (2019) *Nucleic Acids Res* |
| [PhosphoSitePlus](https://www.phosphosite.org/) | Kinase-substrate interactions | **Registration required** | Hornbeck et al. (2015) *Nucleic Acids Res* |

## Expression data

| Source | Description | License | Citation |
|---|---|---|---|
| [CREEDS](https://maayanlab.cloud/Harmonizome/) | Disease expression signatures | CC BY-NC-SA 4.0 | Wang et al. (2016) *Database* |
| [LINCS L1000](https://lincsproject.org/) | Drug perturbation profiles | Public domain | Subramanian et al. (2017) *Cell* |

## Disease genes

| Source | Description | License | Citation |
|---|---|---|---|
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Clinical variant-disease associations | Public domain | - |
| [DisGeNET](https://www.disgenet.org/) | Gene-disease associations | CC BY-NC-SA 4.0 | Piñero et al. (2020) *Nucleic Acids Res* |

## ID mapping

| Source | Description | License |
|---|---|---|
| [HGNC](https://www.genenames.org/) | Gene nomenclature | Free for all use |
| [UniProt](https://www.uniprot.org/) | Protein ID mapping | CC BY 4.0 |

---

**Note:** always cite the original data sources in any publication that uses
output from this pipeline. The CC BY-NC-SA sources above restrict commercial
use, and that restriction propagates to derived results.
