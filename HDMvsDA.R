# DOG ALLERGEN vs HDM extracts
## Single-cell RNA_Seq datasets in R using Seurat
# JEBUNNAHAR MISHU
setwd("/Users/cbm737")

# Install the packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")

install.packages("devtools", mirror = 'usethis')
devtools::install_github("Bioconductor/BiocManager", ref="ghost-binary-repo")
install.packages("remotes")
remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')
#install.packages('Seurat')

BiocManager::install("targets")
BiocManager::install("tidyverse") # declaratively creating graphics
BiocManager::install("ggplot2") # declaratively creating graphics
BiocManager::install("gridExtra") # provides low-level functions for drawing graphical elements.

BiocManager::install("DESeq2")
install.packages("pheatmap")
install.packages("harmony")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
BiocManager::install("ensembldb")


## Loading the libraries
library(BiocManager)
library(devtools)
library(Seurat)
library(SeuratObject)
library(remotes)

library(targets)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(Matrix)
library(pheatmap)
library(harmony)
#library(SingleCellExperiment)




getwd()

# Get or read the file
base_path <- "/Volumes/Groupdir/SUN-ISIM-skin-allergy/Mishu/DAvsHDM"
# Set the directory
setwd(base_path)

Dirs <- list.dirs(path = base_path, full.names = FALSE, recursive = FALSE)
Dirs
list.files(Dirs)


for (x in Dirs) {
  name <- gsub('_filtered_feature_bc_matrix','', x)
  
  # Construct file paths
  mtx_path <- file.path(base_path, x, "matrix.mtx.gz")
  features_path <- file.path(base_path, x, "features.tsv.gz")
  barcodes_path <- file.path(base_path, x, "barcodes.tsv.gz")
  
  # Read the data
  DAHDM_count <- ReadMtx(mtx = mtx_path, 
                         features = features_path, 
                         cells = barcodes_path)
  # create seurat object
  assign(name, CreateSeuratObject(count = DAHDM_count))
  
}
#ls()

### merge datasets
DAHDM_Obj <- merge(HDM2WK, y = DA2WK, add.cell.ids = c("HDM_BAL", "DA_BAL"), project = "HDM_DA") 
str(DAHDM_Obj)
View(DAHDM_Obj@meta.data)
head(DAHDM_Obj@meta.data)


# create the column
DAHDM_Obj$Sample <- rownames(DAHDM_Obj@meta.data)

# splite the column
DAHDM_Obj@meta.data <- separate(DAHDM_Obj@meta.data, col = 'Sample', into = c('Id', 'Type', 'Barcode'), sep = '_')
unique(DAHDM_Obj@meta.data$Id)


# percentage mitocondrial 
DAHDM_Obj$PercentMT <- PercentageFeatureSet(DAHDM_Obj, pattern = "^mt-")

## Visualize data before filtering as a violin plot
VlnPlot(DAHDM_Obj, features = c("nCount_RNA", "nFeature_RNA","PercentMT"), ncol = 3)


#filtering
DAHDM_Obj_filtered <- subset(DAHDM_Obj, subset = nCount_RNA < 15000 & # in filtering, keeping cells which have less then 15000 transcript counts (removing doublet or multiplets cells)and
                               nFeature_RNA > 500 & # cells have greater then 500 genes.(removing low quality or empty cells)
                               PercentMT < 10)# cells have smaller then 10% mitocondrial percentage.(removing low quality or dying cells)


DAHDM_Obj_filtered

## Visualize data before normalization as a violin plot
VlnPlot(DAHDM_Obj_filtered, features = c("nCount_RNA", "nFeature_RNA","PercentMT"), ncol = 3)

# Feature Scatter is used to visualize feature to feature relations
FeatureScatter(DAHDM_Obj_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = lm)



#Split objects
# split the object into a list of multiple objects based on a metadata, column creates a list of two objects
DAHDM_Obj_filtered_list <- SplitObject(DAHDM_Obj_filtered, split.by = "Id")
DAHDM_Obj_filtered_list$HDM
DAHDM_Obj_filtered_list$DA


# Identify common barcodes in the metadata of both Seurat objects
common_barcodes <- intersect(DAHDM_Obj_filtered_list$HDM@meta.data$Barcode, DAHDM_Obj_filtered_list$DA@meta.data$Barcode)
view(common_barcodes)

# Change the default assay to RNA before subsetting
DefaultAssay(DAHDM_Obj_filtered_list$HDM) <- "RNA"
DefaultAssay(DAHDM_Obj_filtered_list$DA) <- "RNA"

# Verify the common barcodes exist in the RNA assay
common_barcodes_in_HDM <- common_barcodes[common_barcodes %in% DAHDM_Obj_filtered_list$HDM@meta.data$Barcode]
common_barcodes_in_DA <- common_barcodes[common_barcodes %in% DAHDM_Obj_filtered_list$DA@meta.data$Barcode]


# Check the lengths of common barcodes found in each object
length(common_barcodes_in_HDM)
length(common_barcodes_in_DA)

# Subset the objects using the verified common barcodes
HDM_common <- subset(DAHDM_Obj_filtered_list$HDM, subset = Barcode == common_barcodes_in_HDM)
DA_common <- subset(DAHDM_Obj_filtered_list$DA, idents = common_barcodes_in_DA)
DAHDM_Obj_filtered_list$HDM@meta.data <- merge(DAHDM_Obj_filtered_list$HDM@meta.data, common_barcodes_in_HDM, by.x = "Barcode", by.y = "emp_id")


# Read file with the information of cells origin
Hashing <- read.csv("Hashing.csv")
Hashing[1:5, ]
View(Hashing)

## Prepare HDM object
## Save the rownames
DAHDM_Obj_filtered_list$HDM@meta.data$Rownames <- row.names(DAHDM_Obj_filtered_list$HDM@meta.data)
#Adding the information about the source of the cells
DAHDM_Obj_filtered_list$HDM@meta.data <- merge(DAHDM_Obj_filtered_list$HDM@meta.data, Hashing, by = 'Barcode')
head(DAHDM_Obj_filtered_list$HDM@meta.data)
rownames(DAHDM_Obj_filtered_list$HDM@meta.data) <- DAHDM_Obj_filtered_list$HDM@meta.data$Rownames

# prepare the DA Object
DAHDM_Obj_filtered_list$DA@meta.data$Hashing <- rep("BAL-1")
head(DAHDM_Obj_filtered_list$DA@meta.data)

### Keep only BAL (bronchoalveolar lavage) type of cells
#DAHDM_Obj_filtered_list$HDM@meta.data <- subset(DAHDM_Obj_filtered_list$HDM@meta.data, subset = Hashing == "BAL-1")
#view(DAHDM_Obj_filtered_list$HDM@meta.data)

DAHDM_Obj_filtered_list$HDM@meta.data[1:5,]


# Check validity of each Seurat object
validObject(DAHDM_Obj_filtered_list$HDM)
validObject(DAHDM_Obj_filtered_list$DA)


## Merge two object where only BAL (bronchoalveolar lavage) type of cells are present
mergedDA_HDM <- merge(DAHDM_Obj_filtered_list$HDM, y = DAHDM_Obj_filtered_list$DA)
view(mergedDA_HDM@meta.data)

### Keep only BAL (bronchoalveolar lavage) type of cells
mergedDA_HDM@meta.data <- subset(mergedDA_HDM@meta.data, subset = Hashing == "BAL-1")
mergedDA_HDM@meta.data <- subset(mergedDA_HDM@meta.data, select = -c(Type, Rownames))


## Removing Tcr  
mergedDA_HDM <- mergedDA_HDM[!grepl('^Tr[abdg][vjc]', rownames(mergedDA_HDM))] # mouse

## Join the layers 
Layers(mergedDA_HDM)
mergedDA_HDM[["RNA"]] <- JoinLayers(mergedDA_HDM)

# Normalize the data
## Find the variables features
## Scale the data to avoid technical noise (batch effect) or biological sources ( different cell cycle)
## perform linear dimension reduction
# Find neighbors, clusters and Run UMAP to visualization
DAHDM_Obj_Norm <- NormalizeData(object = mergedDA_HDM) %>%
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters(resolution = c(0.3, 0.5, 0.8)) %>%
  RunUMAP(dims = 1:15)

str(DAHDM_Obj_Norm)
view(DAHDM_Obj_Norm@meta.data)

# With the Harmony Integration
DAHDM_Obj_HarI <- NormalizeData(object = mergedDA_HDM) %>%
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(group.by.vars = "Id") %>%
  FindNeighbors(dims = 1:15,reduction = "harmony") %>%
  FindClusters(resolution = c(0.3, 0.5, 0.8)) %>%
  RunUMAP(dims = 1:15)

str(DAHDM_Obj_HarI)


## Visualization the clustering without Integration
# UMAP plot
WoHI <- DimPlot(DAHDM_Obj_Norm, reduction = "umap", group.by = "Id",  label = TRUE, label.box = TRUE) + labs(title = "Without Integration")
WoHI
# split.by = "Id", "RNA_snn_res.0.5"
## Visualization the clustering with Integration
# UMAP plot
WHI <- DimPlot(DAHDM_Obj_HarI,reduction = "umap",group.by = "Id", label = TRUE, label.box = TRUE) + labs(title = "With Integration")
WHI
# split.by = "Id", "RNA_snn_res.0.5"
WoHI|WHI

DimPlot(DAHDM_Obj_HarI,reduction = "umap",group.by = "RNA_snn_res.0.5", label = TRUE, label.box = TRUE) + labs(title = "With Integration")


## Visualize specific gene 

## Feature plot - visualize feature expression
FeaturePlot(DAHDM_Obj_HarI, 'Tcf7')

## Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(DAHDM_Obj_HarI, features = c('Thy1', 'Tcf7'))

VlnPlot(DAHDM_Obj_HarI, features = c('Birc5', 'Pclaf'))

# Single cell heatmap of feature expression
DoHeatmap(subset(DAHDM_Obj_HarI, downsample = 100), size = 3)


# DoHeatmap now shows a grouping bar, splitting the heatmap into groups or clusters. This can # #### Try it
# be changed with the `group.by` parameter
DoHeatmap(DAHDM_Obj_HarI, group.by = 'Id', features = VariableFeatures(DAHDM_Obj_HarI),  size = 4,
          angle = 90) + NoLegend()

#cells = 1:500,#

## Find differential expressed feature
clustmarks_DAHDM <- FindAllMarkers(DAHDM_Obj_HarI, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(clustmarks_DAHDM)
dim(clustmarks_DAHDM)


#View top markers for the specific cluster
top10_DAHDM <- clustmarks_DAHDM %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC)
view(top10_DAHDM)




# get the number of cells in each cluster
Cell_Count <- table(Idents(DAHDM_Obj_HarI))
view(Cell_Count)

## write files
write.csv(top10_DAHDM, file = "top10_DAHDM.csv")
write.csv(Cell_Count, file = "Cellcount_DAHDM.csv")
#write.csv(MarkerDAHDM, file = "clustmarks_DAHDM.csv",col.names = NA,row.names = FALSE)

# cluster 0
cluster0.DAHDM <- FindMarkers(DAHDM_Obj_HarI, ident.1 = 0, test.use = "roc", only.pos = TRUE )
VlnPlot(seurat_obj_NormUMAP, features = c("Ramp1", "Ckb"))
view(cluster0.DAHDM)

cluster4.DAHDM <- FindMarkers(DAHDM_Obj_HarI, ident.1 = 4, test.use = "roc", only.pos = TRUE )


FeaturePlot(seurat_obj_NormUMAP, features = c("Ramp1", "Capg"),min.cutoff = 'q10')

# cluster 1
cluster1.DAHDM <- FindMarkers(DAHDM_Obj_HarI, ident.1 = 1, test.use = "roc", only.pos = TRUE)





########### Annotation

# Biomart 
# create mart
library(biomaRt)
listEnsembl()

Mus <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", mirror = "useast")

Mus_DAHDM <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id', 'ensembl_transcript_id', 'gene_biotype'),
                   filters = 'external_gene_name',
                   values = clustmarks_DAHDM$gene,
                   mart = Mus)

View(Mus_DAHDM)

Mus_DAHDM[1:5, ]

## match and keep the genes which are common in mus_mart and cluster_marker

DAHDM_mart <- clustmarks_DAHDM[clustmarks_DAHDM$gene %in%
                                 Mus_DAHDM$external_gene_name, ]


DAHDM_mart[1:5, ]
dim(DAHDM_mart)
view(DAHDM_mart)

write.csv(DAHDM_mart, file = "DA_HDM.csv", row.names = TRUE)

## Cell type identification
DAHDM_mart[1:10, 7]


# Single cell heatmap of feature expression
DoHeatmap(subset(DAHDM_mart, downsample = 100), features = gene, size = 3)


# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
DotPlot(DAHDM_mart, features = gene, split.by = "Id") + RotatedAxis()





