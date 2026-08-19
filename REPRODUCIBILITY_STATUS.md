# Reproducibility status

This repository is a curated release of selected exploratory R scripts. It improves discoverability and provenance but is not yet a complete, one-command reproduction package for the manuscript.

## Coverage matrix

| Analysis module | Code in this release | Inputs available from documented public routes | End-to-end rerun status |
|---|---:|---:|---|
| Single-cell downstream analyses | Yes, selected scripts | Partial | Not validated end to end |
| Human spatial transcriptomics | Yes, selected scripts | Partial | Not validated end to end |
| Spatial segmentation and NMF | Yes, selected scripts | Partial | Not validated end to end |
| hdWGCNA and pathway/TF activity | Yes, selected scripts | Partial | Not validated end to end |
| CellChat analysis | Yes, selected script | Partial | Not validated end to end |
| snATAC-seq processing | No complete workflow | Processed R object only | Not reproducible from raw data here |
| FeSig construction and validation | No complete workflow | No complete analysis-ready inputs | Not reproducible here |
| Mutation/CNV and clinical analyses | No complete workflow | Public cohort routes only | Not reproducible here |
| Drug-response and sorafenib analyses | No complete workflow | Incomplete | Not reproducible here |
| Experimental statistics and source data | No complete workflow | Not archived here | Not reproducible here |

## What has been checked

- Curated `.R` files under `scripts/` are checked for R parse errors by GitHub Actions.
- Original files are retained in `legacy/` for provenance, including one non-executable file that contains R Markdown-style material in an `.R` file.
- The malformed first line of the curated human spatial-transcriptomics script was split into valid R expressions, four case-sensitive package names were normalized, and trailing whitespace was removed, without changing analysis logic.
- Repository documentation distinguishes syntax validation from successful execution and scientific validation.

## Known execution barriers

- Many scripts still contain machine-specific absolute paths.
- Several scripts expect objects created interactively or by earlier, unpublished steps.
- Gene-set files, object schemas, sample manifests, random seeds, and package versions are incomplete.
- Some packages are available only from Bioconductor or GitHub and may require version-specific installation.
- Third-party data downloads and preprocessing are not fully scripted.

## Interpretation boundary

Code availability does not make every manuscript result independently reproducible. The current release is best treated as a transparent record of selected analyses. A complete reproduction package would require immutable input manifests, environment locks, figure-panel mappings, missing workflows, and reload-validated outputs.

## Priorities for a reproducible release

1. Add a figure-to-code-to-data manifest covering every main and supplementary panel.
2. Replace hard-coded paths with `config/paths.local.R` and validate all expected files before analysis.
3. Add the missing snATAC-seq, FeSig, clinical/genomic, drug-response, and experimental-analysis workflows.
4. Record exact package versions and seeds, preferably in a locked R environment.
5. Deposit source data and tagged software releases with separate persistent identifiers.
