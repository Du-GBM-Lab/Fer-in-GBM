# Part 1: Environment Setup and Data Preparation ------------------------------------------------
library(COSG) 
library(harmony) 
library(ggsci) 
library(dplyr) 
library(future) 
library(seurat) 
library(clustree) 
library(cowplot) 
library(data.table) 
library(ggplot2) 
library(patchwork) 
library(stringr) 
library(Matrix) 
library(R.utils) 
library(BiocParallel) 

# Configure parallel processing
plan("multisession", workers = 1)
options(future.globals.maxSize = 1000 * 1024^5)
options(stringsAsFactors = FALSE)
set.seed(123)

# Set working directory
setwd("GSE186344_raw/") # Enter the decompressed directory

# Part 2: Data Loading and Preprocessing ----------------------------------------
# Create sample directory mapping table (handling special characters)
file_df <- data.frame(
  gz_file = list.files(pattern = "_(barcodes|features|matrix)\\.(tsv|mtx)\\.gz$"),
  stringsAsFactors = FALSE
) %>% 
  mutate(
    sample = str_replace(gz_file, "^GSM\\d+_(.*?)_(barcodes|features|matrix).*", "\\1"),
    file_type = case_when(
      str_detect(gz_file, "barcodes") ~ "barcodes.tsv",
      str_detect(gz_file, "features") ~ "features.tsv",
      str_detect(gz_file, "matrix") ~ "matrix.mtx"
    )
  )

# Define target files and compression function
target_files <- c("barcodes.tsv", "features.tsv", "matrix.mtx")

compress_files <- function(dir) {
  for (f in target_files) {
    original_path <- file.path(dir, f)
    gz_path <- paste0(original_path, ".gz")
    if (file.exists(original_path)) {
      gzip(filename = original_path, 
           destname = gz_path, 
           overwrite = TRUE, 
           remove = FALSE)
      cat("Compressed:", original_path, "\n")
    }
  }
}

# Process sample directories
folder_list <- list.dirs(full.names = TRUE, recursive = FALSE)
lapply(folder_list, compress_files)

# Part 3: Seurat Object Creation and QC -----------------------------------------
# Initialize Seurat object list
seurat_list <- list()

# Read and process each sample
for (folder in folder_list) {
  data <- Read10X(data.dir = folder)
  seurat_obj <- CreateSeuratObject(counts = data, project = basename(folder))
  seurat_list[[basename(folder)]] <- seurat_obj
}

# Merge Seurat objects
seurat_merge <- merge(seurat_list[], y = seurat_list[c(2:16)], 
                      add.cell.ids = names(seurat_list))

# Calculate QC metrics
mito_genes <- rownames(seurat_merge)[grep("^mt-", rownames(seurat_merge), ignore.case = T)]
seurat_merge <- PercentageFeatureSet(seurat_merge, features = mito_genes, col.name = "percent_mito")

ribo_genes <- rownames(seurat_merge)[grep("^RP[SL]", rownames(seurat_merge), ignore.case = T)]
seurat_merge <- PercentageFeatureSet(seurat_merge, features = ribo_genes, col.name = "percent_ribo")

hb_genes <- rownames(seurat_merge)[grep("^HB[^(P)]", rownames(seurat_merge), ignore.case = T)]
seurat_merge <- PercentageFeatureSet(seurat_merge, features = hb_genes, col.name = "percent_hb")

# Visualization of QC metrics
feats <- c("nFeature_RNA", "nCount_RNA")
p1 <- VlnPlot(seurat_merge, group.by = "orig.ident", features = feats, 
              pt.size = 0, ncol = 2) + NoLegend()

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2 <- VlnPlot(seurat_merge, group.by = "orig.ident", features = feats,
              pt.size = 0, ncol = 3, same.y.lims = T) + 
  scale_y_continuous(breaks = seq(0, 100, 5)) + 
  NoLegend()

# Filter low-quality cells
selected_mito <- WhichCells(seurat_merge, expression = percent_mito < 25)
selected_ribo <- WhichCells(seurat_merge, expression = percent_ribo > 3)
selected_hb <- WhichCells(seurat_merge, expression = percent_hb < 1)

seurat_merge.filt <- subset(seurat_merge, cells = selected_mito)
seurat_merge.filt <- subset(seurat_merge.filt, cells = selected_ribo)
seurat_merge.filt <- subset(seurat_merge.filt, cells = selected_hb)

# Post-filter visualization
p1_filtered <- VlnPlot(seurat_merge.filt, group.by = "orig.ident", 
                       features = c("nFeature_RNA", "nCount_RNA"), 
                       pt.size = 0, ncol = 2) + NoLegend()

# Part 4: Data Normalization and Integration ------------------------------------
# Standard workflow
seurat_merge.filt <- NormalizeData(seurat_merge.filt)
seurat_merge.filt <- FindVariableFeatures(seurat_merge.filt)
seurat_merge.filt <- ScaleData(seurat_merge.filt)
seurat_merge.filt <- RunPCA(seurat_merge.filt)

# Harmony integration
seuratobj <- RunHarmony(seurat_merge.filt, "orig.ident")

# Dimensionality reduction
seuratobj <- RunUMAP(seuratobj, dims = 1:15, reduction = "harmony")
seuratobj <- RunTSNE(seuratobj, dims = 1:15, reduction = "harmony")

# Add metadata
seuratobj$cancer <- str_split_fixed(seuratobj$orig.ident, "_", n = 2)[, 1]

# Save processed data
saveRDS(seuratobj, file = 'seuratobj_harmony.rds')
library(irGSEA)
Local-Ferro <- read.gmt(Local-Ferro)
seuratObj <- irGSEA.score(object = seuratObj,assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 1,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, geneset = Local-Ferro, geneid = "symbol",
                             method = c("AUCell","UCell","singscore","ssgsea", "JASMINE", "viper"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')

result.dge <- irGSEA.integrate(object = seuratObj,
                               group.by = "cancer",
                               method = c("AUCell","UCell","singscore","ssgsea", "JASMINE"))
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge,
                                      method = "RRA",
                                      top = 50,
                                      show.geneset = NULL)

irGSEA.heatmap.plot
# 查看RRA识别的在多种打分方法中都普遍认可的差异基因集
geneset.show <- result.dge$RRA %>% 
  dplyr::filter(pvalue <= 0.05) %>% 
  dplyr::pull(Name) %>% unique(.)

scatterplot <- irGSEA.density.scatterplot(object = seuratObj,
                                          method = "ssgsea",
                                          show.geneset = "HALLMARK-E2F-TARGETS",
                                          reduction = "tsne")

scatterplot

















seurat_merge <- merge(seurat_list[[1]], y = seurat_list[c(2:17)], add.cell.ids = names(seurat_list))
colnames(seurat_merge)
# 第四部分：整合细胞注释 -------------------------------------------------------
