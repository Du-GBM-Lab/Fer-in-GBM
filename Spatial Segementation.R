# load required packages
library(tidyverse)
library(SPATA2)
library(SPATAData)
library(dplyr)
library(stringr)


# the output is an object fit for the version of SPATA2 you have installed
#my_object <- updateSpataObject(object = my_object)


# load SPATA2 inbuilt data 
object_t269 <- loadExampleObject("UKF269T", meta = TRUE, process = TRUE)
object_t269 <- updateSpataObject(object = object_t269)
# variables that have been created using spatial segmentation
# curently only one: "histology"
getSpatSegmVarNames(object_t269)
# plot the histology image in coordinates frame
plotImage(object_t269) + 
  ggpLayerFrameByCoords(object_t269)

# grouping variable created with spatial segmentation
plotSurface(object_t269, color_by = "histology") 

object_t269 <- createSpatialSegmentation(object_t269)
# obtain only names from variables created with spatial segmentation
ferroptosis <- getSpatSegmVarNames(object_t269)

ferroptosis
# show only variables that have been created by spatial segmentation
getMetaDf(object_t269) %>% 
  select(barcodes, all_of(ferroptosis))
plotSurface(object_t269, color_by = "ferroptosis", clrp_adjust = c("unnamed" = "lightgrey"))
# histology variable is a grouping option ...
getGroupingOptions(object_t269)
# only histology
plotSurface(object = object_t269, pt_alpha = 0)

# histological grouping
plotSurface(object = object_t269, color_by = "ferroptosis")


object_t269 <-
  runCNV(
    object = object_t269,
    # the directory must exist in your file system
    directory_cnv_folder = "data-gbm269/cnv-results", 
    cnv_prefix = "Chr"
  )
cnv_results <- getCnvResults(object_t269)

names(cnv_results)
# if you want the complete object
plotCnvHeatmap(object = object_t269, across = "ferroptosis", clrp = "npg")


# a more complex set up
plotCnvHeatmap(
  object = object_t269, 
  across = "ferroptosis", 
  across_subset = c("fffffff", "unnamed"), # dont show the transition part
  meta_vars = "SPP1",  # visualize SNAP25 expression on the left
  meta_vars_clrs = c("SNAP25" = "inferno"), # with the inferno color spectrum 
  chrom_subset = c("1", "21", "7", "10", "14", "15", "18", "22"), # only show these chromosomes
  ggpLayers = list(arm = list(legendNone())) # remove the chrom arm legend
)

# right plot
plotSurface(object_t269, color_by = "SNAP25", pt_clrsp = "inferno")
plotCnvLineplot(
  object = object_t269,
  across = "ferroptosis", 
  n_bins_bcsp = 1000,
  n_bins_genes = 1000,
  nrow = 3
)

# cnv feature names
getCnvFeatureNames(object = object_t269) %>% head()
# are part of all feature names
getFeatureNames(object = object_t269) %>% head()
plotSurface(
  object = object_t269, 
  color_by = "Chr7", 
  pt_clrsp = "Reds"
)

plotSurface(
  object = object_t269, 
  color_by = "Chr10", 
  pt_clrsp = "Oslo" 
)


plotSpatialTrajectories(
  object = object_t269, 
  ids = "horizontal_mid", 
  color_by = "Chr7", 
  pt_clrsp = "Reds 3"
)

plotSpatialTrajectories(
  object = object_t269, 
  ids = "horizontal_mid", 
  color_by = "Chr10", 
  pt_clrsp = "Oslo"
)

plotStsLineplot(
  object = object_t269, 
  id = "horizontal_mid", 
  variables = "Chr7" 
)

plotStsLineplot(
  object = object_t269, 
  id = "horizontal_mid", 
  variables = "Chr10", 
  line_color = "blue" 
)



getGroupingOptions(object_t269)
object_t269 <- runDEA(object = object_t269, across = "ferroptosis", method = "wilcox")
# extract the complete data.frame
dea_df <- 
  getDeaResultsDf(
    object = object_t269, 
    across = "ferroptosis"
  )

nrow(dea_df)
head(dea_df)
# e.g. top 10 genes for histology area 'tumor' 
getDeaResultsDf(
  object = object_t269, 
  across = "ferroptosis",
  across_subset = "fffffff", # the group name(s) of interest,
  n_highest_lfc = 10, # top ten genes
  max_adj_pval = 0.01 # pval must be lower or equal than 0.01
)
hm <- 
  plotDeaHeatmap(
    object = object_t269, 
    across = "ferroptosis",
    clrp = "npg",
    n_highest_lfc = 40, # subset genes
    n_bcs = 100
  )
hm

# top 9 markers for transition area
transition_markers <- 
  getDeaGenes(object_t269, across = "ferroptosis", across_subset = "fffffff", n_lowest_pval = 9)

plotSurfaceComparison(
  object = object_t269, 
  color_by = transition_markers,
  #pt_clrsp = color_vector("npg")[2], # plot cluster color against white
  outline = TRUE,
  nrow = 3
) + 
  legendNone()
