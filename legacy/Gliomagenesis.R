library(Seurat)
library(SeuratData)
library(SeuratObject)
library(SPATA2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(irGSEA)
library(msigdbr)
library(SeuratData)
library(RcppML)
library(dittoSeq)
brain_injury <- readRDS('GSE278511_brain.injury.samples.RData')
dim(braininjury)
colnames(brain_injury@meta.data)
colnames(tumourigenesis@meta.data)
tumourigenesis <- readRDS('GSE278511_tumourigenesis.atlas.RData')
dim(tumourigenesis)
table(tumourigenesis$cell_type)
malignant_cells<-readRDS('GSE278511_malignant.cells.Sox2CEPPT.and.NestinCPPT.models.RData')
dim(malignant_cells)
DimPlot(brain_injury,split.by = 'sample_group')
DimPlot(tumourigenesis,split.by = 'sample_group')
DimPlot(brain_injury)
DimPlot(malignant_cells)
table(malignant_cells$cell_type)
table(malignant_cells$sample_group)
DimPlot(malignant_cells,split.by = 'sample_group')
DimPlot(malignant_cells)
DimPlot(tumourigenesis)
table(brain_injury@meta.data$sample)
table(brain_injury@meta.data$sample_group)
DimPlot(tumourigenesis,split.by = 'sample')
table(tumourigenesis@meta.data$sample)
table(tumourigenesis@meta.data$sample_group)
table(tumourigenesis@meta.data$mouse_model)

UMAP_celltype <- DimPlot(tumourigenesis, reduction ="umap",
                         group.by="cell_type_fine_annotation",label = T);UMAP_celltype
Idents(tumourigenesis) <- tumourigenesis$cell_type_fine_annotation
scCustomize::DimPlot_scCustom(tumourigenesis, figure_plot = TRUE)
scCustomize::DimPlot_scCustom(tumourigenesis, figure_plot = TRUE, group.by='sample_group')
DimPlot(tumourigenesis,group.by='sample_group')
tumourigenesis$sample_group
#### Seurat V5对象 ####
tumourigenesis <- SeuratObject::UpdateSeuratObject(tumourigenesis)
# tumourigenesis2 <- CreateSeuratObject(counts = CreateAssay5Object(GetAssayData(tumourigenesis,
#                                                                            assay = "RNA", 
#                                                                            slot="counts")),
#                                   meta.data = tumourigenesis[[]])
# tumourigenesis2 <- NormalizeData(tumourigenesis2)
Local.Ferro <- read.gmt(Local-Ferro)
tumourigenesis <- irGSEA.score(object = tumourigenesis, assay = "RNA",
                            slot = "data", seeds = 123, 
                            #ncores = 1,
                            min.cells = 3, min.feature = 0,
                            custom = F, geneset = Local.Ferro,
                            geneid = "symbol",
                            method = c("AUCell","UCell","singscore",
                                       "ssgsea", "JASMINE", "viper"),
                            aucell.MaxRank = NULL, 
                            ucell.MaxRank = NULL,
                            kcdf = 'Gaussian')
tumourigenesis<-readRDS('./tumourigenesis-irGSEA.score.RDS')
Assays(tumourigenesis)
tumourigenesis.dge <- irGSEA.integrate(object = tumourigenesis,
                               group.by = "sample_group",
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea", "JASMINE", "viper"))
geneset.show <- tumourigenesis.dge$RRA %>% 
  dplyr::filter(pvalue <= 0.05) %>% 
  dplyr::pull(Name) %>% unique(.)

irGSEA.heatmap.plot1 <- irGSEA.heatmap(object = tumourigenesis.dge, 
                                       method = "RRA",
                                       show.geneset = geneset.show)
irGSEA.heatmap.plot1
DimPlot(tumourigenesis,split.by = 'sample_group')

dittoDimPlot(tumourigenesis, 
             reduction.use = "umap", var = "sample_group", 
             do.label = F, size = 1
             #,color.panel = rep('gray90',110)
             )




scatterplot <- irGSEA.density.scatterplot(object = tumourigenesis[,tumourigenesis@meta.data$sample_group %in% 'Control'],
                                          method = "UCell",
                                          show.geneset = "ferroptosis",
                                          reduction = "umap")
scatterplot
scatterplot <- irGSEA.density.scatterplot(object = tumourigenesis[,tumourigenesis@meta.data$sample_group %in% 'Control'],
                                          method = "UCell",
                                          show.geneset = "ferroptosis",
                                          reduction = "umap")
scatterplot


brain_injury<-readRDS('./brain_injury-irGSEA.score.RDS')
Assays(brain_injury)
brain_injury.dge <- irGSEA.integrate(object = brain_injury,
                                     group.by = "sample_group",
                                     method = c("AUCell","UCell","singscore",
                                                "ssgsea", "JASMINE", "viper"))
geneset.show <- brain_injury.dge$RRA %>% 
  dplyr::filter(pvalue <= 0.05) %>% 
  dplyr::pull(Name) %>% unique(.)

irGSEA.heatmap.plot1 <- irGSEA.heatmap(object = brain_injury.dge, 
                                       method = "RRA",
                                       show.geneset = geneset.show)
irGSEA.heatmap.plot1
DimPlot(tumourigenesis,split.by = 'sample_group')

dittoDimPlot(brain_injury.dge, 
             reduction.use = "umap", var = "sample_group", 
             do.label = F, size = 1
             #,color.panel = rep('gray90',110)
)


table(brain_injury@meta.data$sample_group)

scatterplot <- irGSEA.density.scatterplot(object = brain_injury[,brain_injury@meta.data$sample_group %in% 'Control'],
                                          method = "UCell",
                                          show.geneset = "ferroptosis",
                                          reduction = "ferroptosis")
scatterplot
scatterplot <- irGSEA.density.scatterplot(object = tumourigenesis[,tumourigenesis@meta.data$sample_group %in% 'Control'],
                                          method = "UCell",
                                          show.geneset = "ferroptosis",
                                          reduction = "umap")
scatterplot


scatterplot_Control <- irGSEA.density.scatterplot(object = tumourigenesis[,tumourigenesis$sample_group %in% 'Control'],
                                                  method = "UCell",
                                                  show.geneset = "ferroptosis",
                                                  reduction = "umap",pt.size=0.1)
scatterplot_Control
scatterplot_Preneoplastic <- irGSEA.density.scatterplot(object = tumourigenesis[,tumourigenesis$sample_group %in% 'Preneoplastic'],
                                                        method = "UCell",
                                                        show.geneset = "ferroptosis",
                                                        reduction = "umap")
scatterplot_Preneoplastic
scatterplot_Early_lesion <- irGSEA.density.scatterplot(object = tumourigenesis[,tumourigenesis$sample_group %in% 'Early_lesion'],
                                                       method = "UCell",
                                                       show.geneset = "ferroptosis",
                                                       reduction = "umap")
scatterplot_Early_lesion
scatterplot_Mid_lesion <- irGSEA.density.scatterplot(object = tumourigenesis[,tumourigenesis$sample_group %in% 'Mid_lesion'],
                                                     method = "UCell",
                                                     show.geneset = "ferroptosis",
                                                     reduction = "umap")
scatterplot_Mid_lesion
scatterplot_Endpoint <- irGSEA.density.scatterplot(object = tumourigenesis[,tumourigenesis$sample_group %in% 'Endpoint'],
                                                   method = "UCell",
                                                   show.geneset = "ferroptosis",
                                                   reduction = "umap")
scatterplot_Endpoint




scatterplot_Control+scatterplot_Preneoplastic+scatterplot_Early_lesion+scatterplot_Mid_lesion+scatterplot_Endpoint



FeaturePlot(object = tumourigenesis[,tumourigenesis$sample_group %in% 'Endpoint'],
            "ferroptosis",
            reduction = "umap") 


dittoDimPlot(tumourigenesis, 
             reduction.use = "umap", var = "sample_group", 
             do.label = F, size = 1,
             color.panel = rep('gray90',110)) + NoLegend() + #rep('gray90',110)是因为指定的列patient有110位患者，需要全部赋予颜色
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )


table(tumourigenesis$sample_group)






metadata <- malignant_cells@meta.data

malignant_cells_ncc<-malignant_cells[,malignant_cells$cell_type %in% 'NCC-like']

dim(malignant_cells_ncc)

DimPlot(malignant_cells_ncc)
DimPlot(malignant_cells_ncc,group.by = "sample_group")
DimPlot(malignant_cells_ncc,split.by = "sample_group")
library(scRNAtoolVis)
cellRatio<-cellRatioPlot(object = malignant_cells,
                         sample.name = "orig.ident",
                         celltype.name = "cell_type",
                         flow.curve = 0.5)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
cellRatio
