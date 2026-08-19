# Legacy scripts

This directory preserves byte-for-byte copies of the scripts that were present at the repository root before the 0.1.0 reorganization. Their staged Git blob identifiers were checked against `origin/main` during curation.

The files are retained for provenance and should not be treated as the recommended entry points. They include machine-specific paths, implicit interactive objects, and incomplete dependencies. In particular, `SC adn ST GEMM .R` contains R Markdown-style front matter, code fences, author metadata, and workflow blocks associated with an external study in a file with an `.R` extension. It is therefore intentionally excluded from automated R parsing and should not be reused or relicensed until its provenance and original license have been confirmed.

Use the organized files under [`../scripts/`](../scripts/) for future maintenance. Do not silently edit a legacy file; make a curated copy and document the change instead.
