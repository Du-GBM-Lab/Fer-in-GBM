# Contributing

Contributions that improve reproducibility, documentation, or portability are welcome.

Before opening a pull request:

1. Do not commit controlled, identifiable, licensed, or large research data.
2. Put machine-specific paths in an untracked `config/paths.local.R` file.
3. Preserve the biological unit of analysis and document sample-level exclusions.
4. Record the input object schema, package versions, random seeds, and generated outputs.
5. Run `Rscript environment/check_packages.R` and parse all edited R scripts.
6. Explain whether the change affects numerical results or only organization/documentation.

Bug reports should name the script, R version, package versions, input object type, and complete error message without exposing credentials or sensitive paths.
