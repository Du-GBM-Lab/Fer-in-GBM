required_packages <- c(
  "bigmemory", "BiocParallel", "CellChat", "clusterProfiler", "clustree",
  "corrplot", "COSG", "cowplot", "data.table", "decoupleR", "dittoSeq",
  "dplyr", "export", "fgsea", "forcats", "future", "genenmf", "ggplot2",
  "ggpubr", "ggrastr", "ggsci", "GSVA", "harmony", "hdf5r", "hdWGCNA",
  "HGNChelper", "Hmisc", "igraph", "infercnv", "irGSEA", "irlba", "Matrix",
  "msigdbr", "nichenetr", "OmnipathR", "OneGene", "openxlsx", "patchwork",
  "PerformanceAnalytics", "pheatmap", "plot1cell", "purrr", "R.utils",
  "RColorBrewer", "RcppML", "readr", "readxl", "reticulate", "rlog",
  "scater", "scRNAtoolVis", "semla", "Seurat", "SeuratData", "SeuratObject",
  "singlet", "spacexr", "SPATA2", "SPATAData", "stringr", "tibble", "tidyr",
  "tidyverse", "UCell", "viridis", "WGCNA"
)

available <- vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
status <- data.frame(
  package = required_packages,
  installed = unname(available),
  version = vapply(
    required_packages,
    function(pkg) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(utils::packageVersion(pkg))
      } else {
        NA_character_
      }
    },
    character(1)
  ),
  row.names = NULL
)

print(status, row.names = FALSE)

if (any(!status$installed)) {
  message(
    "\nMissing packages: ",
    paste(status$package[!status$installed], collapse = ", ")
  )
  quit(status = 1L)
}

message("\nAll listed packages are available in the active R library.")
