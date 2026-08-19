# R environment

`check_packages.R` reports whether packages referenced across the curated scripts are installed in the active R library. It does not install packages and does not reconstruct the original analysis environment.

Run from the repository root:

```sh
Rscript environment/check_packages.R
Rscript environment/check_syntax.R
```

The original package versions were not recorded in the repository. A future validated release should create an environment lock only after each workflow has been rerun against its documented inputs; an automatically generated lockfile from an unrelated machine would provide false precision.
