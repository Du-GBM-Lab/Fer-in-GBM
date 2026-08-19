####### Seurat Object obtained from: Mei, Y., Wang, X., Zhang, J. et al. Siglec-9 acts as an immune-checkpoint molecule on macrophages in glioblastoma, restricting T-cell priming and immunotherapy response. Nat Cancer 4, 1273–1291 (2023). https://doi.org/10.1038/s43018-023-00598-9 #######

####### Seurat object for the study was downloaded from https://figshare.com/articles/dataset/Single-cell_and_spatial_transcriptomic_profiling_of_human_glioblastomas/22434341. The name of the file is "GBM.RNA.integrated.24.rds" #######

############Inside R####################

library(dplyr)
library(Seurat)
library(GSVA)
library(tidyverse)
library(ggplot2)
library(dittoSeq)
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
library(HGNChelper)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(future)
library(ggsci)
library(ggrastr)
library(clusterProfiler)
library(irGSEA)
plan("multisession", workers = 16)
options(future.globals.maxSize = 1000 * 1024^5)
options(stringsAsFactors = FALSE)
set.seed(123)
mycol <- c(pal_d3()(7),pal_aaas()(7),pal_uchicago()(7),pal_jama()(7))
Immunotherapy <- readRDS("GBM.RNA.integrated.24.rds")

scRNA<-readRDS('GBM.RNA.integrated.24.rds')

colnames(scRNA@meta.data)
Idents(scRNA)
Idents(scRNA_mouse)
table(scRNA@meta.data$PD1)
table(scRNA@meta.data$diagnosis)
table(scRNA@meta.data$treatment_1)
table(scRNA@meta.data$treatment_3)
table(scRNA$IDH)
table(scRNA$Pt_number,scRNA$diagnosis)
table(scRNA$Pt_number,scRNA$PD1)
p_treatment_1<-DimPlot(scRNA, reduction = "umap",group.by = "treatment_1",raster=FALSE,pt.size = 0.15,label.size = 7,cols = mycol)
p_treatment_1
p_treatment_3<-DimPlot(scRNA, reduction = "umap",group.by = "treatment_3",raster=FALSE,pt.size = 0.15,label.size = 7,cols = mycol)
p_treatment_3
p_PD1<-DimPlot(scRNA, reduction = "umap",group.by = "PD1",raster=FALSE,pt.size = 0.15,label.size = 7,cols = mycol)
p_PD1
p_diagnosis<-DimPlot(scRNA, reduction = "umap",group.by = "diagnosis",raster=FALSE,pt.size = 0.15,label.size = 7,cols = mycol)
p_diagnosis
p_IDH<-DimPlot(scRNA, reduction = "umap",group.by = "IDH",raster=FALSE,pt.size = 0.15,label.size = 7,cols = mycol)
p_IDH
p1 <- DimPlot(scRNA, reduction = "umap", label = TRUE,group.by = 'RNA_snn_res.0.8',
              pt.size = 0.15,label.size = 7,cols = mycol,raster=FALSE)
p1
p2 <- DimPlot(scRNA, reduction = "umap", label = TRUE,group.by = 'anno_ident',
              pt.size = 0.15,label.size = 7,cols = mycol,raster=FALSE)
p2
#Scoring
Local-Ferro <- read.gmt(Local-Ferro)
scRNA<- irGSEA.score(object = scRNA, assay = "RNA", slot = "data", seeds = 123, ncores = 1, min.cells = 3, min.feature = 0, custom = F, 
                     geneset = Local-Ferro, subcategory = NULL, geneid = "symbol", method = c( "AUCell", "UCell", "singscore", "ssgsea"), 
                     aucell.MaxRank = NULL, ucell.MaxRank = NULL, kcdf = 'Gaussian')
result.dge <- irGSEA.integrate(object = scRNA, group.by = "anno_ident", metadata = NULL, col.name = NULL, method = c( "AUCell", "UCell", "singscore", "ssgsea"))
save(scRNA, result.dge,file = 'scRNA_irGSEA.RData')
load('scRNA_irGSEA.RData')
p2 <- DimPlot(scRNA, reduction = "umap", label = TRUE,group.by = 'anno_ident',
              pt.size = 0.15,label.size = 7,cols = mycol,raster=FALSE)
p2
irGSEA.heatmap.plot<-irGSEA.heatmap(object= result.dge,method= "RRA", top= 50,show.geneset= NULL)
irGSEA.heatmap.plot
heatmap.plot <- irGSEA.heatmap(object = result.dge, method = "ssgsea",top = 30,show.geneset = NULL)
heatmap.plot
# ##
# Idents(scRNA)<-scRNA$PD1
# result.dge <- irGSEA.integrate(object = scRNA, group.by = "PD1", 
#                                metadata = NULL, col.name = NULL, 
#                                method = "ssgsea")
# heatmap.plot <- irGSEA.heatmap(object = result.dge, method = "ssgsea",top = 30,show.geneset = NULL)
# heatmap.plot

scatterplot<-irGSEA.density.scatterplot(object = scRNA,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-ALLOGRAFT-REJECTION",
                                          reduction = "umap")
scatterplot
scatterplot_responder<-irGSEA.density.scatterplot(object = scRNA[,scRNA@meta.data$PD1 %in% 'responder'],
                                        method = "UCell",
                                        show.geneset = "HALLMARK-ALLOGRAFT-REJECTION",
                                        reduction = "umap")
scatterplot_responder


scatterplot_nonresponder<-irGSEA.density.scatterplot(object = scRNA[,scRNA@meta.data$PD1 %in% 'nonresponder'],
                                                  method = "UCell",
                                                  show.geneset = "HALLMARK-ALLOGRAFT-REJECTION",
                                                  reduction = "umap")
scatterplot_nonresponder

vlnplot<-irGSEA.vlnplot(object = scRNA,
                          method = c("AUCell", "UCell", "singscore", "ssgsea"),
                          show.geneset = "ferroptosis")
vlnplot
densityheatmap <- irGSEA.densityheatmap(object = scRNA,
                                        method = "ssgsea",
                                        show.geneset = "ferroptosis")
densityheatmap 

Idents(scRNA) <- "PD1"
Idents(scRNA)
irGSEA.vlnplot(object = scRNA,
               method = c("AUCell", "UCell", "singscore", "ssgsea"),
               show.geneset = "ferroptosis")



DimPlot(scRNA,reduction = "umap",raster=FALSE,
        pt.size = 0.15,label.size = 7,cols = mycol
        )|FeaturePlot(scRNA,features = 'AIF1',label = T)|FeaturePlot(scRNA,features = 'SPP1',label = T)

library(scRNAtoolVis)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有74种差异还比较明显的颜色
col_vector <- unlist(mapply(brewer.pal,qual_col_pals$maxcolors, rownames(qual_col_pals))) 
pie(rep(1,length(col_vector)), col = col_vector)

cellRatio<-cellRatioPlot(object = PD1,
                         sample.name = "Pt_number",
                         celltype.name = "anno_ident",
                         flow.curve = 0.5,fill.col = col_vector )+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
cellRatio
table(PD1@meta.data$anno_ident)
plot_cell_fraction(PD1,  
                   #celltypes = c("TAM-MG","TAM-BDM","AC-like","CD4/CD8","Mast","MES-like"), 
                   groupby = "PD1", show_replicate = T, rep_colname = "Pt_number")

plot_cell_fraction(treatment_1,  
                   #celltypes = c("TAM-MG","TAM-BDM","AC-like","CD4/CD8","Mast","MES-like"), 
                   groupby = "treatment_1", show_replicate = T, rep_colname = "Pt_number")
plot_cell_fraction
