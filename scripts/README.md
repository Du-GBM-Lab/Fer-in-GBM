# Analysis scripts

The curated scripts are grouped by their main analytical role:

| Directory | Scope |
|---|---|
| `01_single_cell/` | Brain-metastasis, published GBM, GBmap, gliomagenesis, and inferCNV analyses |
| `02_spatial/` | Human spatial transcriptomics, spatial segmentation, and semla/NMF analyses |
| `03_regulatory/` | hdWGCNA and decoupleR pathway/transcription-factor activity |
| `04_communication/` | Spatial CellChat analysis |

These files were renamed for navigation. The present curation splits a malformed, concatenated first line in `02_spatial/01_human_spatial_transcriptomics.R`, normalizes four case-sensitive package names (`Seurat`, `UCell`, `HGNChelper`, and `CellChat`), and removes trailing whitespace. It does not otherwise claim to refactor or validate the scientific logic.

Before execution, inspect each script for:

- hard-coded paths and filenames;
- objects expected from an earlier interactive session;
- package and version requirements;
- sample exclusions and factor levels;
- random seeds and parallel-worker settings;
- output paths that may overwrite existing files.

The automated workflow checks syntax only. It does not download data, install the complete analysis environment, run the analyses, or compare outputs with manuscript figures.
