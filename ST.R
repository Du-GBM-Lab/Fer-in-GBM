library(Seurat) library(ggplot2) library(stringr) library(harmony) library(patchwork) library(dplyr) library(patchwork) library(hgnchelper) library(future) library(ggsci) library(gsva) library(tidyverse) library(pheatmap) library(clusterProfiler) plan("multisession", workers = 16) options(future.globals.maxSize = 1000 * 1024^5) options(stringsAsFactors = FALSE) 

# Read data - spatial transcriptomics part
## Batch read data
assays <- dir("./10xVisium 2/")
dir <- paste0("10xVisium 2/", assays, "/outs")
samples_name = c('ukf242_t_st','ukf243_t_st','ukf248_t_st','ukf251_t_st', 'ukf255_t_st','ukf259_t_st','ukf260_t_st', "ukf262_t_st","ukf265_t_st", "ukf266_t_st", "ukf269_t_st", "ukf275_t_st", "ukf296_t_st", "ukf304_t_st", "ukf313_t_st","ukf334_t_st")
spatial.list = list()
for(i in 1:length(dir)){
  spatial <- Load10X_Spatial(data.dir = dir[i],slice = samples_name[i])
  spatial$orig.ident <- samples_name[i]
  spatial$group <- str_split(samples_name[i], "_",simplify = T)[,2]
  spatial <- RenameCells(spatial, add.cell.id = samples_name[i])
  spatial[["percent.mt"]] <- PercentageFeatureSet(spatial, pattern = "^mt[-]")
  spatial.list[[samples_name[i]]] <- spatial
}

p_mt.list <- lapply(spatial.list, function(x){
  p1 <- VlnPlot(x, features = "percent.mt") + ggtitle(unique(x$orig.ident)) + theme(legend.position = "none", axis.text.x = element_blank())
  p2 <- SpatialFeaturePlot(x, features = "percent.mt") + theme(legend.position = "right")
  p <- p1|p2
  p
})
wrap_plots(p_mt.list, ncol = 4)

# Integrate samples
strna <- merge(spatial.list[], spatial.list[2:length(spatial.list)])
strna <- strna[,unname(which(colSums(GetAssayData(strna,,assay="Spatial"))!=0))]

# Normalization
strna <- SCTransform(strna, assay = "Spatial", verbose = T)

## Dimensionality reduction and clustering
strna <- RunPCA(strna, assay = "SCT", verbose = T)
ElbowPlot(strna, ndims=30, reduction="pca")
strna = RunHarmony(strna, group.by.vars = "orig.ident",assay.use = "SCT")
strna <- RunUMAP(strna, reduction = "harmony", dims = 1:20)
strna <- FindNeighbors(strna, reduction = "harmony", dims = 1:20)
strna <- FindClusters(strna, verbose = T,resolution = 0.5)

## Visualization
mycol <- c(pal_d3()(7),pal_aaas()(7),pal_uchicago()(7),pal_jama()(7))
p1 <- DimPlot(strna, reduction = "umap", label = TRUE,group.by = 'SCT_snn_res.0.5',pt.size = 0.15,cols = mycol,label.size = 7)
p1
p2 <- SpatialDimPlot(strna,label = T, alpha = c(0.3, 1), facet.highlight = TRUE, ncol = 4, pt.size.factor = 2,label.size = 3.5)
p2
p1/p2

## Deconvolution - RCTD
library(spacexr) 
library(tidyverse) 
library(ggplot2) 
library(Seurat) 
library(dittoSeq) 
library(patchwork) 
library(plot1cell) 
library(reticulate) 
library(ggplot2) 
library(scater) 
library(Seurat) 
library(dplyr) 
library(tidyr) 
library(purrr) 
library(rlog) 
library(tibble) 
library(stringr) 
library(dittoSeq) 
library(patchwork) 
library(hgnchelper) 
library(future) 
plan("multisession", workers = 12) 
options(future.globals.maxSize = 1000 * 1024^5) 
options(stringsAsFactors = FALSE) 
set.seed(123)

# Load reference and process
scrna <- readRDS('data/Azimuth_core_GBmap.rds')
scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(features = VariableFeatures(object = .), verbose = FALSE) %>% 
  RunPCA(features = VariableFeatures(object = .), verbose = FALSE)
sc_counts <- as.matrix(scrna[['RNA']]@counts)
sc_numi = colSums(sc_counts)
scrna$celltype<-scrna$annotation_level_3
scrna$celltype<-gsub('CD4/CD8','CD4-CD8',scrna$celltype)
table(scrna$celltype)
celltype=data.frame(barcode=colnames(scrna), celltype=scrna$celltype)
names(celltype) = c('barcode', 'cell_type')
cell_types = celltype$cell_type; names(cell_types) <- celltype$barcode
cell_types =as.factor(cell_types)
reference = Reference(sc_counts, cell_types, sc_numi)

# Split spatial samples for deconvolution
sto.list <- SplitObject(strna, split.by = "orig.ident")
slice <- names(sto.list)
result.rctd=data.frame()
for(i in slice){
  coords <- strna@images[[i]]@coordinates[,c(4,5)]
  colnames(coords) <- c("x","y")
  query<- SpatialRNA(coords = coords, counts = sto.list[[i]]@assays$Spatial@counts, 
                     nUMI = structure(sto.list[[i]]$nCount_Spatial, names=colnames(sto.list[[i]])))
  RCTD <- create.RCTD(spatialRNA = query, reference = reference, max_cores = 10,cell_min_instance = 20)
  RCTD <- run.RCTD(RCTD, doublet_mode = 'doublet')
  result_df <- RCTD@results$results_df
  result_df <- as.data.frame(result_df)
  result.rctd=result.rctd%>%rbind(result_df)
}
strna1=strna[,rownames(result.rctd)]
strna1 <- AddMetaData(strna1,metadata = result.rctd)
Idents(strna1)=strna1$second_type
p3<-SpatialDimPlot(strna1,label = T, alpha = c(0.3, 1), facet.highlight = TRUE, ncol = 4, pt.size.factor = 2,label.size = 2.5,repel = TRUE)
p3
save(spatial.list, strna,result.rctd, file = 'my_data.RData')
load('my_data.RData')

# Spatial transcriptomics gene set GSVA analysis
# Prepare files
gmtfile <- "./data/h.all.v2024.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

## For time saving, select 2 pathways
# hallmark.list.dect=hallmark.list[c(1,2)]

# Prepare expression matrix
exp<- as.matrix(strna@assays$SCT@counts)

# Perform GSVA analysis
matrix = gsva(exp, hallmark.list, kcdf="Poisson", method="ssgsea", abs.ranking=T, parallel.sz=10)

## Add to Seurat object
strna$ferroptosis=matrix[18,]
save(matrix, file = 'matrix_gsva_ssgsea.RData')
load('matrix_gsva_ssgsea.RData')

## Plotting
rownames(matrix)
strna$DNA_repair=matrix[12,]
p4<-SpatialFeaturePlot(strna,features = 'DNA_repair', alpha = c(0.3, 1), ncol = 8,pt.size.factor = 1.5,)
p4
