library(Seurat)
library(hdWGCNA)
library(WGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)

# Load results from before
load("UKF269TS2.rdata")

SpatialPlot(UKF269T)

## Enable parallel processing for network analysis (optional)
enableWGCNAThreads(nThreads = 16)

# Re-run normalization pipeline (not using SCT method)
UKF269T <- UKF269T %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

# Clustering and UMAP
UKF269T <- FindNeighbors(UKF269T, dims = 1:30)
UKF269T <- FindClusters(UKF269T, verbose = TRUE)
UKF269T <- RunUMAP(UKF269T, dims = 1:30)

DimPlot(UKF269T, label=TRUE, reduction = "umap", group.by = "seurat_clusters") 
SpatialDimPlot(UKF269T, label = TRUE, label.size = 3, group.by = "Region")

# Extract spatial coordinates and add to metadata
coord <- UKF269T@images[["UKF269T"]]@coordinates
UKF269T <- AddMetaData(UKF269T, metadata = coord)
UKF269T@meta.data[1:5,1:12]

# Construct metaspots
UKF269T <- SetupForWGCNA(
  UKF269T,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "vis"
)

UKF269T <- MetaspotsByGroups(
  UKF269T,
  group.by = c("Region"),
  ident.group = "Region",
  assay = 'Spatial'
)

UKF269T <- NormalizeMetacells(UKF269T)

# Set up expression matrix (set group.by and group_name to NULL to include all spots)
UKF269T <- SetDatExpr(
  UKF269T,
  assay='Spatial',
  group.by=NULL,
  group_name = NULL
)

# Test different soft thresholds
UKF269T <- TestSoftPowers(UKF269T)
plot_list <- PlotSoftPowers(UKF269T)
wrap_plots(plot_list, ncol=2)

# Build co-expression network
UKF269T <- ConstructNetwork(
  UKF269T,
  tom_name='test',
  overwrite_tom=TRUE
)

# Plot dendrogram
PlotDendrogram(UKF269T, main='Spatial hdWGCNA dendrogram')

UKF269T <- ModuleEigengenes(UKF269T)
UKF269T <- ModuleConnectivity(UKF269T)

# Reset module names with "M" prefix (optional)
UKF269T <- ResetModuleNames(
  UKF269T,
  new_name = "M"
)

# Get module eigengenes and gene-module assignments
MEs <- GetMEs(UKF269T)
modules <- GetModules(UKF269T)

# Remove grey module
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# Add MEs to Seurat metadata for visualization
UKF269T@meta.data <- cbind(UKF269T@meta.data, MEs)

# Plot using Seurat's DotPlot function
p <- DotPlot(UKF269T, features=mods, group.by = 'Region', dot.min=0.1)

# Customize plot: flip axes, rotate labels, adjust colors
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  xlab('') + ylab('')

p

# Visualize modules on spatial coordinates
SpatialFeaturePlot(
  UKF269T,
  features = mods,
  alpha = c(0.1, 1),
  ncol = 3
)

# Perform UMAP embedding on co-expression network
UKF269T <- RunModuleUMAP(
  UKF269T,
  n_hubs = 5,
  n_neighbors=15,
  min_dist=0.3,
  spread=1
)

# Plot network graph
ModuleUMAPPlot(
  UKF269T,
  edge.alpha=0.5,
  sample_edges=TRUE,
  keep_grey_edges=FALSE,
  edge_prop=0.075, 
  label_hubs=5 
)

library(igraph)

# Plot module networks (default: top 25 hub genes per module)
ModuleNetworkPlot(
  UKF269T,
  outdir = './ModuleNetworks/'
)

# Save results
save(UKF269T, file = 'UKF269T_hdWGCNA.rdata')

# Check the number of genes in each module
table(modules$module)
