# Data access and availability

This document distinguishes data archived by the authors from reused public datasets and from inputs that are not currently distributed through this repository. The GitHub repository contains code and documentation, not the underlying study datasets.

## Author-archived study output

| Resource | Persistent identifier | Current contents | Important limitation |
|---|---|---|---|
| Processed snATAC-seq object | [Zenodo DOI 10.5281/zenodo.18620722](https://doi.org/10.5281/zenodo.18620722) | `Combined_snATAC_Annotated_AddMotifsed_RunChromVARed.rds` | This record does not contain every dataset, intermediate object, figure source-data table, or script used in the manuscript. |

The manuscript draft currently gives accession `18620721`. The public Zenodo record is `18620722`; the manuscript, repository, and submission metadata should use the same verified identifier.

No data license is asserted here because one was not displayed for the current Zenodo record at the time this documentation was prepared. The authors should select and display an appropriate data license on Zenodo.

## Reused public datasets and resources

| Dataset or resource | Identifier or access route used in the study | Repository coverage |
|---|---|---|
| Human GBM single-cell RNA-seq | GEO [GSE235676](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235676) | Selected downstream scripts |
| Gliomagenesis single-cell RNA-seq | GEO [GSE278511](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278511) | Selected downstream and inferCNV scripts |
| Brain-metastasis single-cell RNA-seq | GEO [GSE186344](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186344) | Selected downstream script |
| Published GBM Seurat object | [Figshare DOI 10.6084/m9.figshare.22434341](https://doi.org/10.6084/m9.figshare.22434341) | Selected downstream script; verify the exact file and version before release |
| Human GBM spatial transcriptomics | Ravi *et al.*, *Cancer Cell* (2022) | Selected downstream scripts; add the exact stable data link and version |
| Glioma GEMM single-cell and spatial data | Kloosterman *et al.*, *Cell* (2024) | Original exploratory material retained only in `legacy/`; add the exact stable data link and version |
| GBmap reference atlas | Source publication and project distribution | Selected reference-mapping script; add the exact stable release URL and version |
| Ivy Glioblastoma Atlas Project | [Ivy GAP](https://glioblastoma.alleninstitute.org/) | No complete download/preprocessing workflow in this release |
| TCGA glioblastoma cohort | [NCI Genomic Data Commons](https://portal.gdc.cancer.gov/) | Clinical/genomic workflow not included in this release |
| Rembrandt cohort | Source repository used by the authors | Clinical workflow not included; add the exact accession and download route |
| Single-nucleus ATAC-seq | GEO [GSE240822](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240822) | A processed object is archived on Zenodo; preprocessing code is not included |

Users must follow the access conditions and licenses set by each original data provider. References to external resources do not redistribute or relicense those datasets.

## Inputs that are not currently archived here

- Complete source-data tables underlying all main and supplementary figures.
- Raw and intermediate experimental measurements and associated sample metadata.
- Exact FeSig gene coefficients and the complete clinical-analysis inputs needed to reproduce model construction and validation.
- All intermediate Seurat, inferCNV, NMF, CellChat, and spatial-analysis objects referenced by hard-coded paths in exploratory scripts.
- A machine-readable manifest connecting every figure panel to its input, script, output, and software environment.

## Recommended manuscript wording

> The processed single-nucleus ATAC-seq object used in this study is available on Zenodo at https://doi.org/10.5281/zenodo.18620722. Public datasets were obtained from their original repositories under the identifiers and access routes listed in the accompanying GitHub repository (https://github.com/Du-GBM-Lab/Fer-in-GBM). Selected analysis code is available in that repository under the MIT License. Data and code not covered by these deposits are identified in the repository's reproducibility statement and should not be described as publicly available unless they are deposited before publication.

## Author actions before publication

1. Correct `18620721` to `18620722` everywhere, after confirming that record 18620722 is the intended deposit.
2. Add a README, file manifest, provenance, object-creation description, checksums, R/package versions, and an explicit license to the Zenodo record.
3. Deposit figure-level source-data tables and a de-identified sample metadata dictionary.
4. Supply stable accession numbers or URLs and exact versions for the Ravi, Kloosterman, GBmap, and Rembrandt resources.
5. Archive a tagged GitHub release in a long-term repository and cite its software DOI separately from the data DOI.
