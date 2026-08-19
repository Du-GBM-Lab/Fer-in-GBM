options(timeout=10000) # Set timeout for operations

# Load R packages
library(cellchat)
library(Seurat)
library(tidyverse)
library(patchwork)

# Load results saved from practical 2
load("ukf269ts2.rdata")

p1 <- SpatialPlot(ukf269t, label = TRUE, label.size = 5,group.by = "seurat_clusters")
p2 <- SpatialPlot(ukf269t,pt.size.factor = 0.6)+NoLegend()
p1+p2

# Define regions
ukf269t@meta.data$region<-NA
ukf269t@meta.data$region[ukf269t@meta.data$seurat_clusters %in% c('0','3','4','5')] <- "high"
ukf269t@meta.data$region[ukf269t@meta.data$seurat_clusters %in% c('1','6','2','7')] <- "low"

SpatialPlot(ukf269t, label = TRUE, label.size = 5,group.by = 'region',cols = c('low'='#4b5cc4','high'='#aa0000'))
Idents(ukf269t)<-ukf269t$region
p1<-SpatialPlot(ukf269t, label = TRUE, label.size = 5)
p2<-SpatialPlot(ukf269t,pt.size.factor = 0.6)+NoLegend()
p1+p2

# Since each spot contains multiple cells, use region-based analysis for cell communication
# Get spatial expression matrix
data.input = Seurat::GetAssayData(ukf269t, slot = "data", assay = "SCT")

# Get metadata
meta = data.frame(labels = Idents(ukf269t), row.names = names(Idents(ukf269t)))
unique(meta$labels)

# Get spatial coordinates
spatial.locs = Seurat::GetTissueCoordinates(ukf269t, scale = NULL,cols = c("imagerow", "imagecol"))

# Load scalefactors from spaceranger output
scalefactors = jsonlite::fromJSON(txt = file.path("../ukf269_t_st/outs/spatial", 'scalefactors_json.json'))
spot.size = 65 # 10x Visium spot size: 55μm with 10μm gap
conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)

# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", 
                           coordinates = spatial.locs,
                           spatial.factors = spatial.factors)

# Set reference database
CellChatDB <- CellChatDB.human # Use CellChatDB.mouse for mouse data
showDatabaseCategory(CellChatDB)

# Use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
cellchat@DB <- CellChatDB.use

# Preprocess CellChat object
cellchat <- subsetData(cellchat) # Necessary even when using full database
future::plan("multisession", workers = 24) # Enable parallel processing

# Identify overexpressed genes
cellchat <- identifyOverExpressedGenes(cellchat)

# Identify overexpressed ligand-receptor pairs
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)

# Infer cell-cell communication network
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE, interaction.range = 250,
                              scale.distance = 0.01, contact.dependent = TRUE,
                              contact.range = 100)

# Filter communications (minimum 10 cells per group)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer communication at pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Visualize interaction counts/weights
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count),
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight),
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Heatmap visualization
p1 <- netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
p2 <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")
p1 + p2

# Display significant pathway results
cellchat@netP$pathways
par(mfrow=c(1,1), xpd = TRUE)
pathways.show <- c("SPP1") 
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Visualize signaling network in spatial context
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial",
                    edge.width.max = 2, vertex.size.max = 1,
                    alpha.image = 0.2, vertex.label.cex = 3.5)

# Compute network centrality
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show,
                                  width = 8, height = 2.5, font.size = 10)

# Visualize incoming signals in spatial context
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial",
                    edge.width.max = 2, alpha.image = 0.2,
                    vertex.weight = "incoming", vertex.size.max = 4,
                    vertex.label.cex = 3.5)

# Visualize ligand-receptor pairs
spatialFeaturePlot(cellchat, pairLR.use = "PTN_PTPRZ1", point.size = 1,
                   do.binary = TRUE, cutoff = 0.05, enriched.only = F,
                   color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, pairLR.use = "PPIA_BSG", point.size = 1,
                   do.binary = TRUE, cutoff = 0.05, enriched.only = F,
                   color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, pairLR.use = "MDK_PTPRZ1", point.size = 1,
                   do.binary = TRUE, cutoff = 0.05, enriched.only = F,
                   color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, pairLR.use = "SPP1_ITGAV_ITGB1", point.size = 1,
                   do.binary = TRUE, cutoff = 0.05, enriched.only = F,
                   color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, features = c("SPP1","ITGA5"), point.size = 0.8,
                   color.heatmap = "Reds", direction = 1)

# Bubble plot for ligand-receptor pairs
netVisual_bubble(cellchat, sources.use = NULL, signaling=cellchat@netP$pathways[1:6],
                 targets.use = NULL, remove.isolate = FALSE)

# Chord diagram for top pathways
netVisual_chord_gene(cellchat, sources.use = NULL, targets.use = NULL,
                     signaling=cellchat@netP$pathways[1:6], lab.cex = 0.5,
                     legend.pos.y = 30)

# Gene expression visualization
p = plotGeneExpression(cellchat, signaling = "SPP1")
p
saveRDS(cellchat,file = 'ukf269t-cellchat.rds')

############################## Secondary plotting
# PTN pathway analysis
par(mfrow=c(1,1), xpd = TRUE)
pathways.show <- c("PTN")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial",
                    edge.width.max = 2, vertex.size.max = 1,
                    alpha.image = 0.2, vertex.label.cex = 3.5)

# Compute centrality for PTN
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show,
                                  width = 8, height = 2.5, font.size = 10)

# Visualize PTN signaling in spatial context
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial",
                    edge.width.max = 2, alpha.image = 0.2,
                    vertex.weight = "incoming", vertex.size.max = 4,
                    vertex.label.cex = 3.5)

# PTN ligand-receptor pairs visualization
spatialFeaturePlot(cellchat, pairLR.use = "PTN_PTPRZ1", point.size = 1,
                   do.binary = TRUE, cutoff = 0.05, enriched.only = F,
                   color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, pairLR.use = "PTN_NCL", point.size = 1,
                   do.binary = TRUE, cutoff = 0.05, enriched.only = F,
                   color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, pairLR.use = "PTN_SDC3", point.size = 1,
                   do.binary = TRUE, cutoff = 0.05, enriched.only = F,
                   color.heatmap = "Reds", direction = 1)
