
knitr::opts_chunk$set(echo = TRUE)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(irlba)

work_dir <- "D:/ProgramData/Spatial Gene Expression/Macrophage-mediated myelin recycling fuels brain cancer malignancy/GEO Data processed/scrna-seq/"
setwd(work_dir)
getwd()

gbm.combined_subset<-readRDS('gbm-complete-subset.rds')
DefaultAssay(gbm.combined_subset) <- "RNA"

# Manually annotated cell types based on cell specific markers


## Assay is set to "RNA" as we will compare RNA expression between clusters
DefaultAssay(gbm.combined_subset) <- "RNA"
p3 <- DotPlot(gbm.combined_subset, group.by = "seurat_clusters", features = c("P2ry12", "Cd74", "Olig1" , "Gfap", "Pecam1", "Cd3e", "Itga4", "Mki67", "Des", "Ccr7", "Clec9a"))

p3

DimPlot(gbm.combined_subset, reduction = "umap", label = F)
# Split dataset for upcoming analyses and figures

## Dataset is split into a macrophage and tumor cell dataset
gbm.myeloid <- gbm.combined_subset[, gbm.combined_subset$celltype %in% c("MG", "MDMs")]

DefaultAssay(gbm.myeloid) <- "integrated"
gbm.myeloid <- ScaleData(gbm.myeloid, verbose = FALSE)
gbm.myeloid <- RunPCA(gbm.myeloid, npcs = 30, verbose = FALSE)

gbm.myeloid <- RunUMAP(gbm.myeloid, reduction = "pca", dims = 1:25)
gbm.myeloid <- FindNeighbors(gbm.myeloid, reduction = "pca", dims = 1:25)
gbm.myeloid <- FindClusters(gbm.myeloid, resolution = 0.26, algorithm = 1)

## Prepping tumor cells
gbm.tumor  <- gbm.combined_subset[, gbm.combined_subset$celltype %in% c("Tumor cells")]

DefaultAssay(gbm.tumor) <- "integrated"
gbm.tumor <- ScaleData(gbm.tumor, verbose = FALSE)
gbm.tumor <- RunPCA(gbm.tumor, npcs = 30, verbose = FALSE)

gbm.tumor <- RunUMAP(gbm.tumor, reduction = "pca", dims = 1:30)
gbm.tumor <- FindNeighbors(gbm.tumor, reduction = "pca", dims = 1:30)
gbm.tumor <- FindClusters(gbm.tumor, resolution = 0.3, algorithm = 1)
# Saving data for downstream analyses

save(gbm.myeloid,file="gbm-complete-myeloid.Rda")
save(gbm.tumor,file="gbm-complete-tumor.Rda")



knitr::opts_chunk$set(echo = TRUE)

# Figure 2B visium
library(Seurat)
library(nichenetr)
library(readxl)
library(dplyr)
library(readr)
library(nichenetr)
library(dplyr)
library(hdf5r)
library(ggplot2)


# Set working directory to folder "pre-processed visium" 
setwd("D:/ProgramData/Spatial Gene Expression/Macrophage-mediated myelin recycling fuels brain cancer malignancy/GEO Data - processed/visium")
load("~/surfdrive/Code Availability/Kloosterman_et_al/Kloosterman et al. 2024 - data availibility/GEO Data - processed/visium/visium_merged_complete.rds")
gbm.merged<-readRDS('../visium/visium_merged_complete.rds')


source("D:/ProgramData/Spatial Gene Expression/Macrophage-mediated myelin recycling fuels brain cancer malignancy/Kloosterman_et_al_Cell_2024-main/unique functions/assignLocation.R")
source("D:/ProgramData/Spatial Gene Expression/Macrophage-mediated myelin recycling fuels brain cancer malignancy/Kloosterman_et_al_Cell_2024-main/unique functions/assignSubtype.R")
source("D:/ProgramData/Spatial Gene Expression/Macrophage-mediated myelin recycling fuels brain cancer malignancy/Kloosterman_et_al_Cell_2024-main/unique functions/assignVerhaak.R")
source("D:/ProgramData/Spatial Gene Expression/Macrophage-mediated myelin recycling fuels brain cancer malignancy/Kloosterman_et_al_Cell_2024-main/unique functions/assignRavi.R")

## Assign Ivy GAP defined niches
DefaultAssay(gbm.merged) <- "SCT"
IvyGAP <- read_excel("D:/ProgramData/Spatial Gene Expression/Macrophage-mediated myelin recycling fuels brain cancer malignancy/Kloosterman_et_al_Cell_2024-main/signatures/IVY_gap_signatures.xlsx")

IvyGAP$Gene <- IvyGAP$NAME %>% convert_human_to_mouse_symbols() 
IvyGAP <- IvyGAP[complete.cases(IvyGAP), ]

ct.all <-  subset(IvyGAP$Gene, subset = IvyGAP$Assigned == "CT")
pan.all <-  subset(IvyGAP$Gene, subset = IvyGAP$Assigned == "CTpan")
le.all <-  subset(IvyGAP$Gene, subset = IvyGAP$Assigned == "LE")
mvp.all <-  subset(IvyGAP$Gene, subset = IvyGAP$Assigned == "CTmvp")


gbm.merged <- assignLocation(gbm.merged, ct.features = ct.all, pan.features = pan.all, le.features = le.all, mvp.features = mvp.all)
SpatialDimPlot(gbm.merged, images = c( "Ink4a_Prim_S8_2"), group.by = "Location" ,  stroke = 0,   image.alpha = 0, alpha = 1) 

# Assign most dominant cellular subtype to each spot

## Assign Neftel cellular subset defined niches
Cell_Signatures <- read_excel("D:/ProgramData/Spatial Gene Expression/Macrophage-mediated myelin recycling fuels brain cancer malignancy/Kloosterman_et_al_Cell_2024-main/signatures/Glioblastoma_meta-module_genelist.xlsx")

MES1.features <- convert_human_to_mouse_symbols (Cell_Signatures$MES1)
MES2.features <- convert_human_to_mouse_symbols(Cell_Signatures$MES2)
AC.features <- convert_human_to_mouse_symbols(Cell_Signatures$AC)
OPC.features <- convert_human_to_mouse_symbols(Cell_Signatures$OPC)
NPC1.features <- convert_human_to_mouse_symbols(Cell_Signatures$NPC1)
NPC2.features <- convert_human_to_mouse_symbols(Cell_Signatures$NPC2)

library(nichenetr)
library(dplyr)
# Set working directory to folder "signatures" 
Idents(gbm.combined_subset) <- gbm.combined_subset$celltype
levels(gbm.combined_subset)
gbm.combined_subset <- RenameIdents(object = gbm.combined_subset,  "DCs" = "rest",
                                    "T cells" = "rest",
                                    "Endothelial cells" = "rest",
                                    "Astrocytes" = "rest",
                                    "Pericytes" = "rest",
                                    "MDMs" = "rest",
                                    "MG"= "rest")   
gbm.combined_subset$macrovsrest <- Idents(gbm.combined_subset)

MES1.expression <- AverageExpression(gbm.combined_subset, assays = "RNA", features = MES1.features, group.by = "macrovsrest")
MES1.expression <- MES1.expression$RNA
MES1.expression <- as.data.frame(MES1.expression)

MES2.expression <- AverageExpression(gbm.combined_subset, assays = "RNA", features = MES2.features, group.by = "macrovsrest")
MES2.expression <- MES2.expression$RNA
MES2.expression <- as.data.frame(MES2.expression)

AC.expression <- AverageExpression(gbm.combined_subset, assays = "RNA", features = AC.features, group.by = "macrovsrest")
AC.expression <- AC.expression$RNA
AC.expression <- as.data.frame(AC.expression)

OPC.expression <- AverageExpression(gbm.combined_subset, assays = "RNA", features = OPC.features, group.by = "macrovsrest")
OPC.expression <- OPC.expression$RNA
OPC.expression <- as.data.frame(OPC.expression)

NPC1.expression <- AverageExpression(gbm.combined_subset, assays = "RNA", features = NPC1.features, group.by = "macrovsrest")
NPC1.expression <- NPC1.expression$RNA
NPC1.expression <- as.data.frame(NPC1.expression)

NPC2.expression <- AverageExpression(gbm.combined_subset, assays = "RNA", features = NPC2.features, group.by = "macrovsrest")
NPC2.expression <- NPC2.expression$RNA
NPC2.expression <- as.data.frame(NPC2.expression)

MES1.expression$FC <- MES1.expression$`Tumor cells`/MES1.expression$rest
MES2.expression$FC <- MES2.expression$`Tumor cells`/MES2.expression$rest
AC.expression$FC <- AC.expression$`Tumor cells`/AC.expression$rest
OPC.expression$FC <- OPC.expression$`Tumor cells`/OPC.expression$rest
NPC1.expression$FC <- NPC1.expression$`Tumor cells`/NPC1.expression$rest
NPC2.expression$FC <- NPC2.expression$`Tumor cells`/NPC2.expression$rest

MES1.features <- subset(rownames(MES1.expression), subset = MES1.expression$FC > 0 & MES1.expression$`Tumor cells` > 0.5  & MES1.expression$`Tumor cells` < 5)
MES2.features<- subset(rownames(MES2.expression), subset = MES2.expression$FC > 0 & MES2.expression$`Tumor cells` > 0.1  & MES2.expression$`Tumor cells` < 5)
AC.features <- subset(rownames(AC.expression), subset = AC.expression$FC > 1  & AC.expression$`Tumor cells` > 1 & AC.expression$`Tumor cells` < 100)
OPC.features <- subset(rownames(OPC.expression), subset = OPC.expression$FC > 5 & OPC.expression$`Tumor cells` > 1  & OPC.expression$`Tumor cells` < 100)
NPC1.features <-  subset(rownames(NPC1.expression), subset = NPC1.expression$FC > 5  & NPC1.expression$`Tumor cells` > 1 & NPC1.expression$`Tumor cells` < 100)
NPC2.features <-  subset(rownames(NPC2.expression), subset = NPC2.expression$FC > 5  & NPC2.expression$`Tumor cells` > 1  & NPC2.expression$`Tumor cells` < 100)
DefaultAssay(gbm.merged) <- "SCT"

gbm.merged <- assignSubtype(gbm.merged, MES1.features = MES1.features, MES2.features = MES2.features, AC.features = AC.features, OPC.features = OPC.features, NPC1.features = NPC1.features, NPC2.features = NPC2.features)
SpatialDimPlot(gbm.merged, images = c("Ink4a_Prim_S8_2"), group.by = "Cell.Subtype",  stroke = 0,   image.alpha = 0, alpha = 1) 



# Calculate GPNMBhigh score for each spot and classify cells as GPNMBhigh enriched when score is higher than 0.1 or GPNMBlow when score is equal or lower than 0.1

```{r}
## GPNMBhigh gene signature
Gene <- c("Gpnmb",  "Fabp5",  "Hmox1",  "Spp1" ,  "Arg1")

##Add score to every cell and create new seurat object. Not necessary to create a new object but in case anything goes wrong your original data is still there.
##This part is taken from the cellcycle classification, however I only used the part where Geneset scores were assigned to the cell
##create a score for every group

gbm.merged <- NormalizeData(gbm.merged) # Needs to be normalised 
fdscores <- AddModuleScore(gbm.merged, features= list(c(Gene)), name="Gene",nbin=100)

##define the different groups in your Genelist, add '1' to every groupname. Groupnames can be checked in the metadata -> str(fdscores@meta.data)
groups <- c("Gene1")

##load function to create density values
densMode <- function(x){
  td <- density(x)
  tdx <- td$x
  tdy <- td$y
  minx <- tdx[which(diff(sign(diff(tdy)))==2)]
  peakx <- tdx[which(diff(sign(diff(tdy)))==-2)]
  return(list(minx=minx, maxy=peakx))
}

##For every group determine the thesholds and plot several plots
for (i in groups){
  ##create density plots 
  plot(density(fdscores@meta.data[,i]), main=paste("densityplot",i, sep=""))
  ##classify the cells based on thresholds of 0.1
  gbm.merged@meta.data[,paste("assignedto",i, sep="")] <- "nonclassified"
  gbm.merged@meta.data[which(fdscores@meta.data[,i]>0.1),paste("assignedto",i, sep="")] <- paste("assignedto",i, sep=" ")
  fdscores_llm <- fdscores
}

Idents(gbm.merged) <- gbm.merged@meta.data$assignedtoGene1
new.cluster.ids <- c("GPNMBlow", "GPNMBhigh")

names(new.cluster.ids) <- levels(gbm.merged)
gbm.merged <- RenameIdents(gbm.merged, new.cluster.ids)
gbm.merged$GPNMBhigh_sig <- Idents(gbm.merged)



## Load macrophage subset signatures 

#load("gbm-complete-subset.Rda")
Top30genes <- read_csv("D:/ProgramData/Spatial Gene Expression/Macrophage-mediated myelin recycling fuels brain cancer malignancy/Kloosterman_et_al_Cell_2024-main/signatures/Top30genes.csv")

M01 <-Top30genes[ Top30genes$cluster %in% c("MG1-P2RY12"), ]$gene
M02 <-Top30genes[ Top30genes$cluster %in% c("MG2-TNF"), ]$gene
M03 <-Top30genes[ Top30genes$cluster %in% c("MG3-GPNMB"), ]$gene
M04 <-Top30genes[ Top30genes$cluster %in% c("MG4-MKI67"), ]$gene

B01 <-Top30genes[ Top30genes$cluster %in% c("MDM1-CCR7"), ]$gene
B02 <-Top30genes[ Top30genes$cluster %in% c("MDM2-H2-EB1"), ]$gene
B03 <-Top30genes[ Top30genes$cluster %in% c("MDM3-GPNMB"), ]$gene
B04 <-Top30genes[ Top30genes$cluster %in% c("MDM4-MKI67"), ]$gene

# Set working directory to folder "signatures" 
levels(gbm.combined_subset)
Idents(gbm.combined_subset) <- gbm.combined_subset$celltype
gbm.combined_subset <- RenameIdents(object = gbm.combined_subset,  "DCs" = "rest",
                                    "T cells" = "rest",
                                    "Endothelial cells" = "rest",
                                    "Astrocytes" = "rest",
                                    "Pericytes" = "rest",
                                    "Tumor cells"  = "rest")   
gbm.combined_subset$macrovsrest <- Idents(gbm.combined_subset)

M01.expression <- AverageExpression(gbm.combined_subset, assays = "RNA", features = M01, group.by = "macrovsrest")
M01.expression <- M01.expression$RNA
M01.expression <- as.data.frame(M01.expression)

M02.expression <- AverageExpression(gbm.combined_subset, assays = "RNA", features = M02, group.by = "macrovsrest")
M02.expression <- M02.expression$RNA
M02.expression <- as.data.frame(M02.expression)

M03.expression <- AverageExpression(gbm.combined_subset, assays = "RNA", features = M03, group.by = "macrovsrest")
M03.expression <- M03.expression$RNA
M03.expression <- as.data.frame(M03.expression)

B01.expression <- AverageExpression(gbm.combined_subset, assays = "RNA", features = B01, group.by = "macrovsrest")
B01.expression <- B01.expression$RNA
B01.expression <- as.data.frame(B01.expression)

B02.expression <- AverageExpression(gbm.combined_subset, assays = "RNA", features = B02, group.by = "macrovsrest")
B02.expression <- B02.expression$RNA
B02.expression <- as.data.frame(B02.expression)

B03.expression <- AverageExpression(gbm.combined_subset, assays = "RNA", features = B03, group.by = "macrovsrest")
B03.expression <- B03.expression$RNA
B03.expression <- as.data.frame(B03.expression)

# Make macrophage specific features
M01 <- subset(rownames(M01.expression), subset = M01.expression$`rest` < 1)
M02 <- subset(rownames(M02.expression), subset = M02.expression$`rest` < 1)
M03 <- subset(rownames(M03.expression), subset = M03.expression$`rest` < 1)
B01 <-  subset(rownames(B01.expression), subset = B01.expression$`rest` < 1)
B02 <-  subset(rownames(B02.expression), subset = B02.expression$`rest` < 0.95) # Have to put this threshold to 0.95 as addition of Cxcl16 to the signature makes the function fail
B03 <-  subset(rownames(B03.expression), subset = B03.expression$`rest` < 1)

# Create functions to score each subset based on its ontogenty (MG or MDMs)


MGSubsetScoring <- function (object, 
                             M01.features,  M02.features, M03.features
                             ,set.ident = FALSE, ctrl_genes = 100, ...)
{
  name <- "MG.Subset"
  features <- list(  M01.Score = M01.features, M02.Score =  M02.features, M03.Score = M03.features)
  object.cc <- AddModuleScore(object = object,
                              features = features,
                              name = name,
                              ctrl = min(c(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1))), 
                                         ctrl_genes), 
                              ...)
  cc.columns <- grep(pattern = name, x = colnames(x = object.cc[[]]),
                     value = TRUE)
  cc.scores <- object.cc[[cc.columns]]
  cc.scores <- t(cc.scores)
  cc.scores <- t(apply(cc.scores, 1, function(y)(y-min(y))/(max(y)-min(y))))
  cc.scores <- t(cc.scores)
  rm(object.cc)
  CheckGC()
  assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores,
                                                                 first = "M01", second = "M02", third = "M03", 
                                                                 null = "Undecided") {
    if (all(scores < -0.0)) {
      return(null)
    }
    else {
      return(c(first, second, third)[which(x = scores == max(scores))])
    }
  }
  )
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments),
                     by = 0)
  colnames(x = cc.scores) <- c("rownames", "M01.Score", "M02.Score", "M03.Score",
                               "MG_subset")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("M01.Score", "M02.Score", "M03.Score",
                             "MG_subset")]
  object[[colnames(x = cc.scores)]] <- cc.scores
  if (set.ident) {
    object[["old.ident"]] <- Idents(object = object)
    Idents(object = object) <- "MG_subset"
  }
  return(object)
}

MDMSubsetScoring <- function (object, 
                              B01.features,  B02.features, B03.features
                              ,set.ident = FALSE, ctrl_genes = 100, ...)
{
  name <- "MDM.Subset"
  features <- list(  B01.Score = B01.features, B02.Score =  B02.features, B03.Score = B03.features)
  object.cc <- AddModuleScore(object = object,
                              features = features,
                              name = name,
                              ctrl = min(c(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1))), 
                                         ctrl_genes), 
                              ...)
  cc.columns <- grep(pattern = name, x = colnames(x = object.cc[[]]),
                     value = TRUE)
  cc.scores <- object.cc[[cc.columns]]
  cc.scores <- t(cc.scores)
  cc.scores <- t(apply(cc.scores, 1, function(y)(y-min(y))/(max(y)-min(y))))
  cc.scores <- t(cc.scores)
  rm(object.cc)
  CheckGC()
  assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores,
                                                                 first = "B01", second = "B02", third = "B03", 
                                                                 null = "Undecided") {
    if (all(scores < -0.0)) {
      return(null)
    }
    else {
      return(c(first, second, third)[which(x = scores == max(scores))])
    }
  }
  )
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments),
                     by = 0)
  colnames(x = cc.scores) <- c("rownames", "B01.Score", "B02.Score", "B03.Score",
                               "MDM_subset")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("B01.Score", "B02.Score", "B03.Score",
                             "MDM_subset")]
  object[[colnames(x = cc.scores)]] <- cc.scores
  if (set.ident) {
    object[["old.ident"]] <- Idents(object = object)
    Idents(object = object) <- "MDM_subset"
  }
  return(object)
}

# Use MacrophageSubsetScoring functions


DefaultAssay(gbm.merged) <- "SCT"

gbm.merged <- MGSubsetScoring(gbm.merged, 
                              M01.features = M01,  M02.features = M02, M03.features = M03)

gbm.merged <- MDMSubsetScoring(gbm.merged, 
                               B01.features = B01,  B02.features = B02, B03.features = B03)

Gene <- c("Aif1", "P2ry12", "Lgals3", "Itga1", "Trem2", "Siglech")

##Add score to every cell and create new seurat object. Not necessary to create a new object but in case anything goes wrong your original data is still there.
##This part is taken from the cellcycle classification, however I only used the part where Geneset scores were assigned to the cell
##create a score for every group
fdscores <- AddModuleScore(gbm.merged, features= list(c(Gene)), name="Gene")

##define the different groups in your Genelist, add '1' to every groupname. Groupnames can be checked in the metadata -> str(fdscores@meta.data)
groups <- c("Gene1")

##load function to create density values
densMode <- function(x){
  td <- density(x)
  tdx <- td$x
  tdy <- td$y
  minx <- tdx[which(diff(sign(diff(tdy)))==2)]
  peakx <- tdx[which(diff(sign(diff(tdy)))==-2)]
  return(list(minx=minx, maxy=peakx))
}

##For every group determine the thesholds and plot several plots
for (i in groups){
  ##create density plots 
  plot(density(fdscores@meta.data[,i]), main=paste("densityplot",i, sep=""))
  ##classify the cells based on thresholds
  gbm.merged@meta.data[,paste("assignedto",i, sep="")] <- "nonclassified"
  gbm.merged@meta.data[which(fdscores@meta.data[,i]>0.1),paste("assignedto",i, sep="")] <- paste("assignedto",i, sep=" ")
  fdscores_llm <- fdscores
}

Idents(gbm.merged) <- gbm.merged@meta.data$assignedtoGene1
new.cluster.ids <- c( "Rest", "Macrophage")

names(new.cluster.ids) <- levels(gbm.merged)
gbm.merged <- RenameIdents(gbm.merged, new.cluster.ids)
gbm.merged$Macrophage_score <- Idents(gbm.merged)



## Localise in pan-macrophages the MG 

## LLM 18 gene signature
Gene <- c("P2ry12", "Siglech")
##Add score to every cell and create new seurat object. Not necessary to create a new object but in case anything goes wrong your original data is still there.
##This part is taken from the cellcycle classification, however I only used the part where Geneset scores were assigned to the cell
##create a score for every group
fdscores <- AddModuleScore(gbm.merged, features= list(c(Gene)), name="Gene")

##define the different groups in your Genelist, add '1' to every groupname. Groupnames can be checked in the metadata -> str(fdscores@meta.data)
groups <- c("Gene1")

##load function to create density values
densMode <- function(x){
  td <- density(x)
  tdx <- td$x
  tdy <- td$y
  minx <- tdx[which(diff(sign(diff(tdy)))==2)]
  peakx <- tdx[which(diff(sign(diff(tdy)))==-2)]
  return(list(minx=minx, maxy=peakx))
}

##For every group determine the thesholds and plot several plots
for (i in groups){
  ##create densityplots and set cut-offs
  vl <- densMode(fdscores@meta.data[,i])[1]
  vl2 <- densMode(fdscores@meta.data[,i])[2]
  plot(density(fdscores@meta.data[,i]), main=paste("densityplot",i, sep=""))
  if(density(fdscores@meta.data[,i])$y[which(density(fdscores@meta.data[,i])$x==vl$minx[1])] > 6 || density(fdscores@meta.data[,i])$y[which(density(fdscores@meta.data[,i])$x==vl$minx[1])] <0.01){
    abline(v=vl$minx[2])
    threshold=vl$minx[2]
  }else{
    abline(v=vl$minx[1], col="red")
    threshold=vl$minx[1]
  }
  ##classify the cells based on thresholds
  gbm.merged@meta.data[,paste("assignedto",i, sep="")] <- "nonclassified"
  gbm.merged@meta.data[which(fdscores@meta.data[,i]>threshold),paste("assignedto",i, sep="")] <- paste("assignedto",i, sep=" ")
  fdscores_llm <- fdscores
}

Idents(gbm.merged) <- gbm.merged@meta.data$assignedtoGene1
new.cluster.ids <- c( "Rest", "MG")

names(new.cluster.ids) <- levels(gbm.merged)
gbm.merged <- RenameIdents(gbm.merged, new.cluster.ids)
gbm.merged$MG_score <- Idents(gbm.merged)

# Classify each cell based on the pan-macrophage and microglia marker (Macrophage - Rest = MDM, Macrophage - MG = MG)


Idents(gbm.merged) <- paste(gbm.merged$Macrophage_score, gbm.merged$MG_score)
levels(gbm.merged)
gbm.merged <- RenameIdents(object = gbm.merged, "Rest Rest" = "TME" ,
                           "Macrophage MG" = "MG",
                           "Macrophage Rest" = "MDM",
                           "Rest MG" = "TME" )
levels(gbm.merged)
gbm.merged$celltype <- Idents(gbm.merged)

Idents(gbm.merged) <- paste(gbm.merged$celltype,gbm.merged$MDM_subset, gbm.merged$MG_subset)
levels(Idents(gbm.merged))


gbm.merged <- RenameIdents(object = gbm.merged, "TME B01 M02" = "NA", 
                           "MDM B01 M01" = "MDM1-CCR7",
                           "MDM B01 M02" = "MDM1-CCR7",
                           "MDM B03 M01" = "MDM3-GPNMB",
                           "MG B01 M02" = "MG2-TNF" ,
                           "TME B01 M01" = "NA",
                           "MG B03 M01" = "MG1-P2RY12",
                           "TME B02 M01" = "NA",
                           "MDM B01 M03" = "MDM1-CCR7",
                           "TME B03 M01" = "NA",
                           "TME B02 M02" = "NA",
                           "MDM B03 M03" = "MDM3-GPNMB",
                           "MG B01 M01" = "MG1-P2RY12" ,
                           "MDM B02 M01" = "MDM2-H2-EB1",
                           "TME B03 M02" ="NA",
                           "MG B02 M01" = "MG1-P2RY12" ,
                           "MG B03 M02" = "MG2-TNF",
                           "MDM B03 M02" = "MDM3-GPNMB",
                           "TME B01 M03" = "NA",
                           "MDM B02 M02" = "MDM2-H2-EB1",
                           "TME B03 M03" = "NA",
                           "MG B02 M02" = "MG2-TNF",
                           "MG B03 M03" = "MG3-GPNMB",
                           "MG B01 M03" = "MG3-GPNMB",
                           "MDM B02 M03" = "MDM2-H2-EB1",
                           "TME B02 M03" = "NA",
                           "MG B02 M03" = "MG3-GPNMB"
)


gbm.merged$TME <- Idents(gbm.merged)
levels(Idents(gbm.merged))

my_levels <-  c("MG1-P2RY12",	"MG2-TNF", "MG3-GPNMB" ,	"MDM1-CCR7",	"MDM2-H2-EB1", "MDM3-GPNMB", "NA")

# Relevel object@ident
Idents(gbm.merged) <- factor(Idents(gbm.merged), levels = my_levels)
gbm.merged$TME <- Idents(gbm.merged)
Idents(gbm.merged) <- gbm.merged$TME

SpatialDimPlot(gbm.merged,  stroke = 0,   image.alpha = 0, alpha = 1, images = c("Ink4a_Prim_S8_2"), cells.highlight = CellsByIdentities(gbm.merged), facet.highlight = TRUE,  ncol = 4)
SpatialDimPlot(gbm.merged,  stroke = 0,   image.alpha = 0, alpha = 1, images = c("Ink4a_Rec_S2"), cells.highlight = CellsByIdentities(gbm.merged), facet.highlight = TRUE,  ncol = 4)
```


# Assign tumor area based on RNA count data (>7000 = tumor tissue, lower = non tumor bearing (Ntb))

gbm.merged@meta.data[,"tissue"] <- c("Ntb")
gbm.merged@meta.data[which(gbm.merged$nCount_Spatial>7000),"tissue"] <-  c("Tumor")
Idents(gbm.merged) <- gbm.merged$tissue
VlnPlot(gbm.merged, features = "nCount_Spatial", pt.size = 0.1, group.by = c("orig.ident"), split.by = "tissue") + NoLegend()


# Figure 2B: Classification of spots based on glioblastoma cellular subtype dominance or anatatomical niche


SpatialDimPlot(gbm.merged, group.by = "Location", images = c("Ink4a_Prim_S8_2")) + scale_fill_manual(values = c("orange", "violetred1", "darkolivegreen3","lightblue", "grey")) 
SpatialDimPlot(gbm.merged, group.by = "Cell.Subtype", images = c("Ink4a_Prim_S8_2"))


saveRDS(gbm.merged, file = "~/Desktop/visium_merged_complete.rds")








---
  title: "Kloosterman and Erbani et al., Cell 2024: Figure 2D-G"
author: "Daan J. Kloosterman"
date: "15/03/24"
output: html_document
editor_options: 
  chunk_output_type: console
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Figure 2D-G - Co-localisation analysis of macrophage subsets, glioblastoma cellular subtype and anatamical niches

Using the modules applied on the VISIUM 10x dataset in the script for Figure 2B, we will perform a co-localisation analysis to understand which cell types are more probable to co-localise within the tumor. In addition, we will quantify the composition of cell types and compare their dynamics between primary and recurrent tumors.

# Loading packages/libraries required for the data pre-processing

```{r}
library(Seurat)
library(nichenetr)
library(readxl)
library(dplyr)
library(readr)
library(nichenetr)
library(dplyr)
library(corrplot)
library(PerformanceAnalytics)
library(RColorBrewer)
library(ggplot2)
```


# Load data

# Figure 2D (1): Data generated to create correlation plot between different spot assignments. Heatmap generated in this chuck was remodelled in illustrator for aesthetic purposes.

gbm.merged_tumor <- gbm.merged[, gbm.merged$tissue %in% c("Tumor")]
gbm.merged_tumormac <- gbm.merged_tumor[, gbm.merged_tumor$Macrophage_score %in% c("Macrophage")]

Idents(gbm.merged_tumormac) <- gbm.merged_tumormac$orig.ident
gbm.merged_tumormac <- RenameIdents(object = gbm.merged_tumormac, "S1_2" = "S9",
                                    "S2_2" = "S10",
                                    "S3_2" = "S11",
                                    "S4_2" = "S12",
                                    "S8_2" = "S16"
)
gbm.merged_tumormac$orig.ident <- Idents(gbm.merged_tumormac)
correlation <- gbm.merged_tumormac@meta.data


## Make LLM score
DefaultAssay(gbm.merged_tumormac) <- "SCT"

## GPNMBhigh gene signature
Gene <- c("Gpnmb",  "Fabp5",  "Hmox1",  "Spp1" ,  "Arg1")

##Add score to every cell and create new seurat object. Not necessary to create a new object but in case anything goes wrong your original data is still there.
##This part is taken from the cellcycle classification, however I only used the part where Geneset scores were assigned to the cell
##create a score for every group

gbm.merged_tumormac <- NormalizeData(gbm.merged_tumormac) # Needs to be normalised 
fdscores <- AddModuleScore(gbm.merged_tumormac, features= list(c(Gene)), name="Gene",nbin=100)

##define the different groups in your Genelist, add '1' to every groupname. Groupnames can be checked in the metadata -> str(fdscores@meta.data)
groups <- c("Gene1")

##load function to create density values
densMode <- function(x){
  td <- density(x)
  tdx <- td$x
  tdy <- td$y
  minx <- tdx[which(diff(sign(diff(tdy)))==2)]
  peakx <- tdx[which(diff(sign(diff(tdy)))==-2)]
  return(list(minx=minx, maxy=peakx))
}

##For every group determine the thesholds and plot several plots
for (i in groups){
  ##create densityplots
  plot(density(fdscores@meta.data[,i]), main=paste("densityplot"))
  ##classify the cells based on thresholds of 0.1
  gbm.merged@meta.data[,paste("assignedto",i, sep="")] <- "nonclassified"
  gbm.merged@meta.data[which(fdscores@meta.data[,i]>0.1),paste("assignedto",i, sep="")] <- paste("assignedto",i, sep=" ")
  fdscores_gpnmbhigh <- fdscores
}

correlation_gpnmbhigh <- fdscores_gpnmbhigh@meta.data # acquired from Figure 1k

correlation_perwell <- data.frame(correlation$MES1.Score, correlation$MES2.Score, correlation$AC.Score, correlation$NPC1.Score, correlation$NPC2.Score, correlation$OPC.Score,  correlation$PAN.Score, correlation$LE.Score, correlation$MVP.Score, correlation$CT.Score,  correlation$M01.Score, correlation$M02.Score, correlation$M03.Score,correlation$B01.Score, correlation$B02.Score, correlation$B03.Score, correlation_gpnmbhigh$Gene1)


corrplot(cor(correlation_perwell[-c(631),]),  
         tl.col = "black", # Labels color
         bg = "white",     # Background color
         col=brewer.pal(n=8, name="PRGn"))       # Color palette

chart.Correlation(cor(correlation_perwell[-c(631),]), histogram=TRUE, pch=20)

write.csv2(cor(correlation_perwell[-c(631),]), "~/Desktop/Correlation_Figure2D.csv")
rm(gbm.merged_tumor)
rm(gbm.merged_tumormac)
```


# Figure 2D (2): Assessing the fold-change of the various cell types and tumor areas between primary and recurrent PDG-Ink4a and PDG-shp53 mice.

```{r}
gbm.merged_quantify <- gbm.merged[, gbm.merged$tissue %in% c("Tumor")]
Idents(gbm.merged_quantify) <- gbm.merged_quantify$orig.ident
gbm.merged_quantify <- RenameIdents(object = gbm.merged_quantify, "S1_2" = "S9",
                                    "S2_2" = "S10",
                                    "S3_2" = "S11",
                                    "S4_2" = "S12",
                                    "S8_2" = "S16"
)
gbm.merged_quantify$orig.ident <- Idents(gbm.merged_quantify)

table_TME <- table(gbm.merged_quantify$TME, gbm.merged_quantify$orig.ident)
table_CellSubtype <- table(gbm.merged_quantify$Cell.Subtype, gbm.merged_quantify$orig.ident)
table_IVY <- table(gbm.merged_quantify$Location, gbm.merged_quantify$orig.ident)
table_GPNMBhigh <- table(gbm.merged_quantify$GPNMBhigh_sig, gbm.merged_quantify$orig.ident)

library(openxlsx)
write.xlsx(table_TME, file="~/Desktop/table_fractions-TME.xlsx", sheetName="table_TME", append=TRUE, rowNames = F)
write.xlsx(table_CellSubtype, file="~/Desktop/table_fractions-Neftel.xlsx", sheetName="table_CellSubtype", append=TRUE, rowNames = F)
write.xlsx(table_IVY, file="~/Desktop/table_fractions-IVY.xlsx", sheetName="table_IVY", append=TRUE, rowNames = F)
write.xlsx(table_GPNMBhigh, file="~/Desktop/table_fractions-GPNMBhigh.xlsx", sheetName="table_GPNMBhigh", append=TRUE, rowNames = F)
```


# Figure 2E-G Representative visualization of VISIUM 10X spatial transcriptomic analyses in recurrent PDG-Ink4a glioblastoma, highlighting GPNMBhigh deserted and enriched areas

```{r}
Idents(gbm.merged) <- gbm.merged$TME
my_levels <-  c("MG1-P2RY12",	"MG2-TNF", "MG3-GPNMB" ,	"MDM1-CCR7",	"MDM2-H2-EB1", "MDM3-GPNMB", "NA")

# Relevel object@ident
Idents(gbm.merged) <- factor(Idents(gbm.merged), levels = my_levels)
gbm.merged$TME <- Idents(gbm.merged)

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=6)
color_list <- c(color_list, "grey")

Idents(gbm.merged) <- gbm.merged$TME

SpatialDimPlot(gbm.merged, group.by = "TME",   stroke = 0,   image.alpha = 0, alpha = 1, 
               images = c("Ink4a_Prim_S8_2")) + scale_fill_manual(values = color_list)

SpatialDimPlot(gbm.merged,  stroke = 0,   image.alpha = 0, alpha = 1, images = c("Ink4a_Prim_S8_2"), cells.highlight = CellsByIdentities(gbm.merged), facet.highlight = TRUE,  ncol = 4)

## Figure 2E-G
SpatialDimPlot(gbm.merged[, gbm.merged$TME %in% c("MG1-P2RY12",	"MG2-TNF", "MG3-GPNMB" ,	"MDM1-CCR7",	"MDM2-H2-EB1", "MDM3-GPNMB")], group.by = "GPNMBhigh_sig",   stroke = 0,   image.alpha = 1, alpha = 1, 
               images = c("Ink4a_Rec_S2")) + scale_fill_manual(values =  c("#FFA07A", "#AB82FF"))
SpatialDimPlot(gbm.merged, group.by = "Cell.Subtype",stroke = 0,  images = c("Ink4a_Rec_S2"), image.alpha = 0) + scale_fill_manual(values = c("darkolivegreen3", "pink", "tomato3","darkgoldenrod1", "darkgoldenrod3", "royalblue3", "NA" ))
SpatialDimPlot(gbm.merged, group.by = "Location",stroke = 0,  images = c("Ink4a_Rec_S2"), image.alpha =  0) + scale_fill_manual(values = c("orange", "violetred1", "darkolivegreen3","lightblue", "grey")) 

## Supplementary Figure 7C-E (additional representitive images of primary tumor)
SpatialDimPlot(gbm.merged[, gbm.merged$TME %in% c("MG1-P2RY12",	"MG2-TNF", "MG3-GPNMB" ,	"MDM1-CCR7",	"MDM2-H2-EB1", "MDM3-GPNMB")], group.by = "GPNMBhigh_sig",   stroke = 0,   image.alpha = 1, alpha = 1, 
               images = c("Ink4a_Prim_S8_2")) + scale_fill_manual(values =  c("#FFA07A", "#AB82FF"))
SpatialDimPlot(gbm.merged, group.by = "Cell.Subtype",stroke = 0,  images = c("Ink4a_Prim_S8_2"), image.alpha = 0) + scale_fill_manual(values = c("darkolivegreen3", "pink", "tomato3","darkgoldenrod1", "darkgoldenrod3", "royalblue3", "NA" ))
SpatialDimPlot(gbm.merged, group.by = "Location",stroke = 0,  images = c("Ink4a_Prim_S8_2"), image.alpha =  0) + scale_fill_manual(values = c("orange", "violetred1", "darkolivegreen3","lightblue", "grey")) 
```
















