# Fer-in-GBM

## Spatially resolved ferroptosis-associated states link mesenchymal-like, macrophage-rich niches to adverse clinical features in glioblastoma

**Status:** manuscript submitted to *Briefings in Bioinformatics*.

This project examines ferroptosis-associated heterogeneity in glioblastoma by integrating single-cell RNA sequencing, spatial transcriptomics, single-nucleus chromatin accessibility, clinicogenomic analyses, and experimental validation.

The repository provides selected R scripts for single-cell, spatial, regulatory, and cell-cell communication analyses. It also documents the current boundaries of data and code availability. Spatial co-localization and expression-score associations are not, by themselves, direct measurements of ferroptosis or evidence of causation.

## Resources

- [GitHub repository](https://github.com/Du-GBM-Lab/Fer-in-GBM)
- [Data access and availability](https://github.com/Du-GBM-Lab/Fer-in-GBM/blob/main/DATA.md)
- [Reproducibility status](https://github.com/Du-GBM-Lab/Fer-in-GBM/blob/main/REPRODUCIBILITY_STATUS.md)
- [Analysis scripts](https://github.com/Du-GBM-Lab/Fer-in-GBM/tree/main/scripts)
- [Citation metadata](https://github.com/Du-GBM-Lab/Fer-in-GBM/blob/main/CITATION.cff)
- [Zenodo data record](https://doi.org/10.5281/zenodo.18620722)

## Reproducibility note

The current code release is not a complete one-command workflow. Several scripts retain hard-coded local paths or depend on objects created in earlier exploratory sessions. Complete workflows for snATAC-seq preprocessing, FeSig construction and validation, clinicogenomic and drug-response analyses, experimental statistics, and all figure source data are not yet included. See the repository's reproducibility statement before reusing the code.

## Authors

Wenbo Zhao, Qingbei Lian, Yihang Yang, and Jianyang Du.

Code is available under the [MIT License](https://github.com/Du-GBM-Lab/Fer-in-GBM/blob/main/LICENSE). Third-party datasets remain subject to their original terms.
