#remotes::install_github("ludvigla/semla") 
library(tibble)
library(semla)
library(tibble)

# Load data
#load("ukf269ts2.rdata") 

# Create Seurat object
ukf269 <- Load10X_Spatial(
  data.dir = "./ukf269/outs/",  # Data directory
  filename = "filtered_feature_bc_matrix.h5",  # h5 filename
  slice = "ukf269")  # h&e image name
ukf269$orig.ident <- "ukf269"
SpatialFeaturePlot(ukf269, features = "nCount_Spatial")

# Convert Seurat to semla object
ukf269 <- UpdateSeuratForSemla(ukf269)

# Load h&e images
ukf269 <- LoadImages(ukf269)

# Normalize and find variable features
ukf269 <- ukf269 |> 
  NormalizeData() |> 
  FindVariableFeatures()

cols <- viridis::rocket(11, direction = -1)

# Plot with semla
MapFeatures(ukf269, features = "nFeature_Spatial", image_use = "raw", override_plot_dims = TRUE, pt_size = 2) & theme_legend_right()

# Function MapFeatures for feature visualization
cols <- viridis::rocket(11, direction = -1)
MapFeatures(ukf269, features = c("KRT2", "SLPI"), colors = cols, pt_size = 1.5)

# Function MapFeaturesSummary() adds subplots (box/violin/density) on the right
MapFeaturesSummary(ukf269, features = "SLPI", subplot_type = "violin", pt_size = 1.4, colors = cols)

# Plot multiple features
MapMultipleFeatures(ukf269, pt_size = 2, max_cutoff = 0.99, override_plot_dims = TRUE, 
                    features = c('SLPI', 'COL1A2','IGHGP','DCD'))

# Crop image
MapFeatures(ukf269, section_number = 1, features = "SLPI", image_use = "raw", 
            color = cols, pt_size = 1.4, pt_alpha = 0.5) & 
  labs(x="x-axis", y="y-axis") & 
  theme(panel.grid.major = element_line(linetype = "dashed"),
        axis.text = element_text(),
        axis.title = element_text())

MapFeatures(ukf269, features = c("SLPI"), image_use = "raw", pt_size = 3, 
            section_number = 1, color = cols, crop_area = c(0.5, 0.3, 0.72, 0.6))

library(singlet)

# Run NMF
ukf269 <- RunNMF(ukf269)
k <- ncol(ukf269@reductions$nmf@feature.loadings)
k

# Plot cross-validation results to find optimal value
RankPlot(ukf269, reduction = "nmf")

# Spatial visualization
MapFeatures(ukf269, features = paste0("nmf_", 1:6), override_plot_dims = TRUE, 
            pt_size = 1.5, colors = viridis::magma(n = 11, direction = -1)) & 
  theme(plot.title = element_blank())

factor_colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
                   '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
                   '#9a6324', '#fffac8', '#800000', '#aaffc3')

# Select factors 1-6
selected_factors <- c(1:6)
MapMultipleFeatures(ukf269, features = paste0("nmf_", selected_factors), 
                    colors = factor_colors, image_use = "raw", 
                    override_plot_dims = TRUE, pt_size = 4)

PlotFeatureLoadings(ukf269, dims = 1:6, reduction = "nmf", nfeatures = 20, 
                    mode = "dotplot", fill = "darkmagenta", pt_size = 3)

## Gene contribution heatmap
PlotFeatureLoadings(ukf269, dims = selected_factors, reduction = "nmf", 
                    nfeatures = 10, mode = "heatmap", 
                    gradient_colors = viridis::magma(n = 11, direction = -1))
