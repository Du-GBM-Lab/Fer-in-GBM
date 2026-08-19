library(decoupleR)
library(OmnipathR)
library(Seurat)
library(dplyr)
library(pheatmap)

# Load results saved from before
load("UKF269TS2.rdata")
p1 <- SpatialPlot(UKF269T, label = TRUE, label.size = 5, group.by = "seurat_clusters")
p2 <- SpatialPlot(UKF269T, pt.size.factor = 0.6) + NoLegend()
p1 + p2
# Region definition
UKF269T@meta.data$Region <- NA
UKF269T@meta.data$Region[UKF269T@meta.data$seurat_clusters %in% c('0','3','4','5')] <- "High"
UKF269T@meta.data$Region[UKF269T@meta.data$seurat_clusters %in% c('1','6','2','7')] <- "Low"
SpatialPlot(UKF269T, label = TRUE, label.size = 5, group.by = 'Region', cols = c('Low'='#4b5cc4','High'='#AA0000'))

Idents(UKF269T) <- UKF269T$Region  

net <- decoupleR::get_progeny(organism = 'human', top = 500)
net

# Extract normalized transformed counts
mat <- as.matrix(UKF269T@assays$SCT@data)

# Run mlm
acts <- decoupleR::run_mlm(mat = mat, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor = 'weight', 
                           minsize = 5)

# Extract mlm and store in pathwaysmlm assay
UKF269T[['pathwaysmlm']] <- acts %>%
  tidyr::pivot_wider(id_cols = 'source', 
                     names_from = 'condition',
                     values_from = 'score') %>%
  tibble::column_to_rownames(var = 'source') %>%
  Seurat::CreateAssayObject(.)

# Switch assay
Seurat::DefaultAssay(object = UKF269T) <- "pathwaysmlm"

# Scale the data
UKF269T <- Seurat::ScaleData(UKF269T)
UKF269T@assays$pathwaysmlm@data <- UKF269T@assays$pathwaysmlm@scale.data
# Visualization
# Spatial visualization of pathway activity
p1 <- Seurat::SpatialPlot(UKF269T, 
                          group.by = "Region") + ggplot2::ggtitle('Region')

# Visualize MAPK pathway (related to apoptosis)
p2 <- Seurat::SpatialPlot(UKF269T, features = c("Trail")) + ggplot2::ggtitle('Trail activity')
p2 <- Seurat::SpatialPlot(UKF269T, features = c("JAK-STAT")) + ggplot2::ggtitle('JAK-STAT activity')
p0 <- Seurat::SpatialPlot(UKF269T, features = c("TNFa")) + ggplot2::ggtitle('TNFa activity')

p1 + p2 + p0

# UMAP visualization of pathway activity
p3 <- Seurat::DimPlot(UKF269T, 
                      reduction = "umap", 
                      label = TRUE, 
                      pt.size = 0.5) + 
  Seurat::NoLegend() + 
  ggplot2::ggtitle('Region')

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p4 <- Seurat::FeaturePlot(UKF269T, features = c("TNFa")) + 
  ggplot2::scale_colour_gradient2(low = colors, mid = 'white', high = colors) +
  ggplot2::ggtitle('TNFa activity')
p3 + p4

# Extract pathway activities from object
df <- t(as.matrix(UKF269T@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  dplyr::mutate(cluster = Seurat::Idents(UKF269T)) %>%
  tidyr::pivot_longer(cols = -cluster, 
                      names_to = "source", 
                      values_to = "score") %>%
  dplyr::group_by(cluster, source) %>%
  dplyr::summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  tidyr::pivot_wider(id_cols = 'cluster', 
                     names_from = 'source',
                     values_from = 'mean') %>%
  tibble::column_to_rownames(var = 'cluster') %>%
  as.matrix()

# Color scale
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-1.25, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 1.25, length.out = floor(100 / 2)))

# Plot
pheatmap::pheatmap(mat = top_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 20,
                   cellheight = 20,
                   treeheight_row = 20,
                   treeheight_col = 20) 

# Save results
save(UKF269T, file = 'UKF269T_Pathway_activity.rdata')

# BiocManager::install("OmnipathR", force = TRUE)

# Update OmnipathR if errors occur
net1 <- decoupleR::get_collectri(organism = 'human', 
                                 split_complexes = FALSE)

net1

# Remove NA values in net to prevent errors
Net2 <- subset(net1, !is.na(source) & !is.na(target))

# Extract normalized transformed counts
mat <- as.matrix(UKF269T@assays$SCT@data)

# Run ulm
acts <- decoupleR::run_ulm(mat = mat, 
                           net = Net2, 
                           .source = 'source', 
                           .target = 'target',
                           .mor='mor', 
                           minsize = 5)

acts

# Extract ulm and store in tfsulm assay
UKF269T[['tfsulm']] <- acts %>%
  tidyr::pivot_wider(id_cols = 'source', 
                     names_from = 'condition',
                     values_from = 'score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Switch assay to 'tfsulm'
DefaultAssay(object = UKF269T) <- "tfsulm"

# Scale the data
UKF269T <- Seurat::ScaleData(UKF269T)
UKF269T@assays$tfsulm@data <- UKF269T@assays$tfsulm@scale.data

# Spatial visualization of TF activity
p1 <- Seurat::SpatialPlot(UKF269T, 
                          group.by = "Region") + ggplot2::ggtitle('Region')

# Visualize STAT3 transcription factor activity
p2 <- Seurat::SpatialPlot(UKF269T, features = c("STAT3")) + ggplot2::ggtitle('STAT3 activity')

p1 + p2

# UMAP visualization of TF activity
p3 <- Seurat::DimPlot(UKF269T, 
                      reduction = "umap", 
                      label = TRUE, 
                      pt.size = 0.5) + 
  Seurat::NoLegend() + 
  ggplot2::ggtitle('Region')

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p4 <- Seurat::FeaturePlot(UKF269T, features = c("STAT3")) + 
  ggplot2::scale_colour_gradient2(low = colors, mid = 'white', high = colors) +
  ggplot2::ggtitle('STAT3 activity')

## Check corresponding expression levels
DefaultAssay(object = UKF269T) <- "SCT"
p5 <- Seurat::FeaturePlot(UKF269T, 
                          features = c("STAT3")) + 
  ggplot2::ggtitle('STAT3 expression')

p3 + p4 + p5
n_tfs <- 50

# Extract TF activities
df <- t(as.matrix(UKF269T@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  dplyr::mutate(cluster = UKF269T@meta.data$Region) %>%
  tidyr::pivot_longer(cols = -cluster, 
                      names_to = "source", 
                      values_to = "score") %>%
  dplyr::group_by(cluster, source) %>%
  dplyr::summarise(mean = mean(score))

# Retrieve highly variable TFs across spatial regions
tfs <- df %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(mean)) %>%
  dplyr::arrange(-abs(std)) %>%
  head(n_tfs) %>%
  dplyr::pull(source)

# Extract top TFs long dataframe and convert to wide matrix
top_acts_mat <- df %>%
  dplyr::filter(source %in% tfs) %>%
  tidyr::pivot_wider(id_cols = 'cluster', 
                     names_from = 'source',
                     values_from = 'mean') %>%
  tibble::column_to_rownames('cluster') %>%
  as.matrix()

# Color palette
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))

# Heatmap visualization
pheatmap::pheatmap(mat = top_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 20,
                   treeheight_col = 20) 

# UMAP visualization of TF activity
p3 <- Seurat::DimPlot(UKF269T, 
                      reduction = "umap", 
                      label = TRUE, 
                      pt.size = 0.5) + 
  Seurat::NoLegend() + 
  ggplot2::ggtitle('Region')

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p4 <- Seurat::FeaturePlot(UKF269T, features = c("MYCN")) + 
  ggplot2::scale_colour_gradient2(low = colors, mid = 'white', high = colors) +
  ggplot2::ggtitle('MYCN activity')

## Check corresponding expression levels
DefaultAssay(object = UKF269T) <- "SCT"
p5 <- Seurat::FeaturePlot(UKF269T, 
                          features = c("MYCN")) + 
  ggplot2::ggtitle('MYCN expression')

p3 + p4 + p5
