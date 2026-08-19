library(Seurat)
library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(irGSEA)
library(msigdbr)
library(SeuratData)
library(RcppML)
library(tidyverse) 
library(infercnv) 
library(harmony) 
library(export) 
library(OneGene) 
library(ggpubr)
plan("multisession", workers = 20)
options(future.globals.maxSize = 1000 * 1024^6)
options(stringsAsFactors = FALSE)
braininjury <- readRDS('GSE278511_brain.injury.samples.RData')
DimPlot(braininjury, split.by = 'sample_group')
tumourigenesis <- readRDS('GSE278511_tumourigenesis.atlas.RData')
tumourigenesis[['cell_type']] <- tumourigenesis@active.ident # Add cell type annotations to metadata
colnames(tumourigenesis@meta.data)
DimPlot(tumourigenesis, split.by = 'sample_group')
dim(tumourigenesis)
table(tumourigenesis$sample_group)
# Subset data
#seurat_obj = tumourigenesis[,tumourigenesis$sample_group %in% c("Control" ,"Early_lesion","Endpoint","Mid_lesion")]
seurat_obj = tumourigenesis[,tumourigenesis$sample_group %in% c("Control","Endpoint")]
dim(seurat_obj)
table(seurat_obj$sample_group)

DimPlot(seurat_obj, split.by = 'sample_group')
DimPlot(tumourigenesis, group.by = "sample_group")
exprMatrix <- as.matrix(GetAssayData(seurat_obj, slot = 'counts'))
# Extract cell type information 
cellAnnota <- subset(seurat_obj@meta.data, select = c('sample_group')) 
write.table(cellAnnota, "./CNV/groupFiles.txt", sep = '\t', col.names = FALSE)

#--------------------- Create InferCNV object ------------------------- 
# ref_group_names parameter is set based on cell annotation file (Control cells as reference)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exprMatrix, # single-cell counts matrix
                                    annotations_file="./CNV/groupFiles.txt", # cell grouping file
                                    delim="\t",
                                    gene_order_file= "./CNV/cnv_ref/mouse_gencode.GRCm39.vM32.basic.annotation.by_gene_name.infercnv_positions", # gene position file
                                    ref_group_names= c('Control') # Use Control as reference
)
dir.create('CNV') # Create directory
gc() 

# Check memory usage
memory.size()

# Check memory limits
memory.limit()
rm(exprMatrix, seurat_obj, tumourigenesis)
infercnv_obj = infercnv::run(infercnv_obj, # inferCNV object
                             cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir= './CNV/', # output directory
                             cluster_by_groups= T, # cluster by groups first then hierarchical clustering
                             hclust_method="ward.D2",
                             plot_steps=F,
                             num_threads = 30,
                             write_expr_matrix = T, # important!!
                             output_format = "pdf", # use PDF instead of PNG (memory efficient)
                             res = 200, useRaster = TRUE
)

tumourigenesis <- readRDS('GSE278511_tumourigenesis.atlas.RData')
# Subset Early_lesion data
#seurat_obj = tumourigenesis[,tumourigenesis$sample_group %in% c("Control" ,"Early_lesion","Endpoint","Mid_lesion")]
seurat_obj_Early_lesion = tumourigenesis[,tumourigenesis$sample_group %in% c("Control","Early_lesion")]
exprMatrix <- as.matrix(GetAssayData(seurat_obj_Early_lesion, slot = 'counts'))
# Extract cell type information
dir.create('CNV_Early_lesion') # Create directory
gc() 
cellAnnota <- subset(seurat_obj_Early_lesion@meta.data, select = c('sample_group')) 
write.table(cellAnnota, "./CNV_Early_lesion/groupFiles.txt", sep = '\t', col.names = FALSE)
table(cellAnnota$sample_group)
table(seurat_obj_Early_lesion$sample_group)

#--------------------- Create InferCNV object ------------------------- 
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exprMatrix, # single-cell counts matrix
                                    annotations_file="./CNV_Early_lesion//groupFiles.txt", # cell grouping file
                                    delim="\t",
                                    gene_order_file= "./CNV/cnv_ref/mouse_gencode.GRCm39.vM32.basic.annotation.by_gene_name.infercnv_positions", # gene position file
                                    ref_group_names= c('Control') # Use Control as reference
)
dir.create('CNV_Early_lesion') # Create directory
gc() 

# Check memory usage
memory.size()

# Check memory limits
memory.limit()
rm(exprMatrix, tumourigenesis)
infercnv_obj = infercnv::run(infercnv_obj, # inferCNV object
                             cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir= './CNV_Early_lesion//', # output directory
                             cluster_by_groups= T, # cluster by groups first then hierarchical clustering
                             hclust_method="ward.D2",
                             plot_steps=F,
                             num_threads = 30,
                             write_expr_matrix = T, # important!!
                             output_format = "pdf", # use PDF instead of PNG (memory efficient)
                             res = 200, useRaster = TRUE
)
