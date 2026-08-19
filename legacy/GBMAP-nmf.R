.libPaths(c('~/r/x86_64-pc-linux-gnu-library/4.4', '/refdir/rlib', '/usr/local/lib/r/library')) 
.libPaths("~/r/x86_64-pc-linux-gnu-library/4.4/") 
getwd() 
setwd("/home/data/t060202/spatial-gsva-duijian/programs origins of myeloid/") 
.libPaths() 
library(Seurat) 
library(patchwork) 
library(ggplot2) 
library(bigmemory) 
library(patchwork) 
library(tidyr) 
library(dplyr) 
library(RColorBrewer) 
library(ggsci) 
library(ggrastr) 
library(clusterProfiler) 
library(genenmf) 
library(BiocParallel) 
library(viridis) 

custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))))(10), rev(magma(323, begin = 0.18))) 

bpparam() # Check current parallelization parameters
plan("multisession", workers = 16) 
options(future.globals.maxSize = 1000 * 1024^5) 
options(stringsAsFactors = FALSE) 
set.seed(123) 

mycol <- c(pal_d3()(7),pal_aaas()(7),pal_uchicago()(7),pal_jama()(7)) 

gbmap <- readRDS("/home/data/t060202/spatial-gsva-duijian/gbmap-data/azimuth_core_gbmap.rds") 
colnames(gbmap@meta.data) 
table(gbmap$annotation_level_1) 
table(gbmap$annotation_level_2) 
table(gbmap$annotation_level_3) 
table(gbmap$annotation_level_4) 

# Inspect object structure
str(gbmap, max.level = 2) 
# Check if assay type is assay5

DimPlot(gbmap,group.by = 'annotation_level_2') 

# table(gbmap$patient) 
# DimPlot(gbmap,group.by = 'patient') 

# select_cell <- rownames(gbmap@meta.data)[as.vector(gbmap@meta.data[,"annotation_level_1"]) %in% c("neoplastic")] 
# seu <- subset(gbmap, cells = select_cell) 

seu <- gbmap 
rm(gbmap) 
colnames(seu@meta.data) 
DimPlot(seu, group.by="annotation_level_3") 
DimPlot(seu, group.by="annotation_level_4") 
DimPlot(seu) 

DefaultAssay(seu) <- "RNA" 
colnames(seu@meta.data) 

seu.list <- SplitObject(seu, split.by = "patient") 

genenmf.programs <- MultiNMF(seu.list, assay="RNA", k=4:9, min.exp = 0.05) 
genenmf.metaprograms <- getMetaprograms(genenmf.programs, metric = "cosine", weight.explained = 0.5, nmp=6) 
# genenmf.metaprograms <- getMetaprograms(genenmf.programs, nmp=4) 

ph <- plotMetaprograms(genenmf.metaprograms,palette = custom_magma, similarity.cutoff = c(0.1,1)) 
ph 

# States in high-impact articles
plotMetaprograms(genenmf.metaprograms,similarity.cutoff = c(0.1,1), palette = c('white','#bfd3e6','#9ebcda','#8c96c6', '#8c6bb1','#88419d','#810f7c','#4d004b')) 

genenmf.metaprograms_jaccard <- getMetaprograms(genenmf.programs, metric = "jaccard", weight.explained = 0.5,nmp=5) 
ph_j <- plotMetaprograms(genenmf.metaprograms_jaccard,palette = custom_magma,similarity.cutoff = c(0, 0.3)) 
ph_j 

genenmf.metaprograms$metaprograms.metrics 
rm(seu,seu.list) 
genenmf.metaprograms$metaprograms.metrics 
lapply(genenmf.metaprograms$metaprograms.genes, head) 
genenmf.metaprograms$metaprograms.genes.weights$mp1 

# Annotation step
library(msigdbr) 
library(fgsea) 

top_p <- lapply(genenmf.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(seu), category = "C5", subcategory = "GO:BP")
}) 
head(top_p$mp1) 
head(top_p$mp2) 
head(top_p$mp3) 
head(top_p$mp4) 

library(ucell) 
mp.genes <- genenmf.metaprograms$metaprograms.genes 
seu <- AddModuleScore_UCell(seu, features = mp.genes, slot = "data",name = "", ncores=20) 

VlnPlot(seu, features=names(mp.genes), group.by = "annotation_level_2", pt.size = 0, ncol=3,raster=FALSE) 

matrix <- seu@meta.data[,names(mp.genes)] 
# dimred <- scale(matrix) 
dimred <- as.matrix(matrix) 
colnames(dimred) <- paste0("mp_",seq(1, ncol(dimred))) 

# Create new dim reduction
seu@reductions[["mpsignatures"]] <- new("DimReduc", cell.embeddings = dimred, assay.used = "RNA", key = "mp_", global = FALSE) 

set.seed(123) 
seu <- RunUMAP(seu, reduction="mpsignatures", dims=1:length(seu@reductions[["mpsignatures"]]), metric = "euclidean", reduction.name = "umap_mp") 

colnames(seu@meta.data) 
DimPlot(seu, reduction = "umap_mp", group.by = "patient") + theme(aspect.ratio = 1) 
DimPlot(seu, reduction = "umap_mp", group.by = "annotation_level_2") + theme(aspect.ratio = 1) 
DimPlot(seu, reduction = "umap_mp", group.by = "annotation_level_3") + theme(aspect.ratio = 1) 
DimPlot(seu, reduction = "umap_mp", group.by = "celltype_original") + theme(aspect.ratio = 1) 
DimPlot(seu, reduction = "umap_mp", group.by = "annotation_level_1") + theme(aspect.ratio = 1) 

library(viridis) 
FeaturePlot(seu, features = names(mp.genes), reduction = "umap_mp", ncol=3) & scale_color_viridis(option="A") & theme(aspect.ratio = 1, axis.text=element_blank(), axis.ticks=element_blank()) 

# Extract metadata for mp1-mp5
# Load packages
library(Hmisc) 
library(pheatmap) 
library(RColorBrewer) 

# Extract data (example)
mp_data <- seu@meta.data[, c("mp1", "mp2", "mp3", "mp4", "mp5")] 

# Convert to numeric matrix (prevent factor issues)
mp_matrix <- as.matrix(apply(mp_data, 2, as.numeric)) 

# Calculate correlation matrix and p-values (handle missing values)
cor_result <- rcorr(mp_matrix, type = "pearson") 

# Extract results
cor_matrix <- cor_result$r 
p_matrix <- cor_result$p 

# Create significance stars matrix
signif_stars <- ifelse(p_matrix < 0.001, "***", 
                       ifelse(p_matrix < 0.01, "**", 
                              ifelse(p_matrix < 0.05, "*", ""))) 

# Create combined label matrix
label_matrix <- matrix(
  paste0(
    sprintf("%.2f", cor_matrix), # Keep 2 decimals
    "\n", # New line
    signif_stars # Add stars
  ),
  nrow = nrow(cor_matrix),
  dimnames = dimnames(cor_matrix)
) 

# Set color palette
color_palette <- colorRampPalette(c("#3794bf", "white", "#df8640"))(100) 

# Plot heatmap
pheatmap(
  mat = cor_matrix,
  color = color_palette,
  display_numbers = label_matrix, # Show combined labels
  number_color = "black", # Text color
  fontsize_number = 8, # Font size
  cluster_rows = FALSE, # Disable clustering
  cluster_cols = FALSE,
  border_color = "gray60", # Border color
  main = "MP1-5 Correlation Matrix\n(Coefficient with Significance Stars)",
  breaks = seq(-1, 1, length.out = 100), # Color breaks
  legend_breaks = c(-1, -0.5, 0, 0.5, 1), # Legend ticks
  legend_labels = c("-1.0", "-0.5", "0", "0.5", "1.0"),
  angle_col = 45, # Column label angle
  cellwidth = 40, # Cell width
  cellheight = 40 # Cell height
) 
