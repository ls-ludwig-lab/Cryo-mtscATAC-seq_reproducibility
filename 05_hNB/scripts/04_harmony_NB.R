# Load Libraries
library(Signac)             
library(Seurat)             
library(ggplot2)            
library(tidyverse)          
library(dplyr)              
library(devtools)           
library(patchwork)          
library(SummarizedExperiment) 
library(GenomeInfoDb)       
library(GenomicRanges)      
library(harmony)            

# Load Neuroblastoma Objects
NB01 <- readRDS("path_to/NB01.rds")
NB02 <- readRDS("path_to/NB02.rds")

# ============================
# Merge Seurat Objects
# ============================
seurat.list = list(NB01, NB02)
sample.names <- c("NB01", "NB02")
seurat_objects <- seurat.list

## Merge the Seurat objects into a single object
seurat = Reduce(function(x, y) merge(x, y, add.cell.ids = sample.names), seurat_objects)

# ============================
# Normalize and Scale RNA Assay
# ============================
DefaultAssay(seurat) <- "RNA"
seurat <- NormalizeData(seurat)
seurat <- ScaleData(seurat)

# ============================
# Dimensionality Reduction after Merging
# ============================
selected_resolutions <- c(0.5)

## TF-IDF and Latent Semantic Indexing (LSI)
seurat <- RunTFIDF(seurat)                         # Standard preprocessing for scATAC-seq
seurat <- FindTopFeatures(seurat, min.cutoff = "q0") # Select all features (peaks)
seurat <- RunSVD(seurat)                            # SVD for LSI dimensionality
seurat <- RunUMAP(seurat, reduction = 'lsi', dims = 2:30) # UMAP embedding
seurat <- FindNeighbors(seurat, reduction = 'lsi', dims = 2:30)
seurat <- FindClusters(seurat, resolution = 0.3)

## Normalize RNA assay again with scaling by median counts
seurat <- NormalizeData(
  object = seurat,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat$nCount_RNA)
)

# ============================
# Harmony Integration
# ============================
## Define parameters for Harmony
ResolutionClustering = 0.5
var.batch = "Sample"        # Batch variable
lambdaHarmony = 1           # Harmony lambda parameter
AlgoClustering = 3          # Seurat clustering algorithm (Louvain = 1, SLM = 3)

## Convert batch variable to factor
seurat@meta.data[[var.batch]] = as.factor(seurat@meta.data[[var.batch]])

## Run Harmony integration
seurat = RunHarmony(
  object = seurat,
  group.by.vars = "Sample", 
  reduction = "lsi", 
  assay.use = "peaks", 
  project.dim = FALSE, 
  lambda = lambdaHarmony, 
  reduction.save = "harmony"
)

## UMAP, neighbors, and clustering on Harmony embeddings
seurat = RunUMAP(object = seurat, reduction = 'harmony', dims = 2:30)
seurat = FindNeighbors(object = seurat, reduction = 'harmony', dims = 2:30)
seurat = FindClusters(object = seurat, verbose = FALSE, resolution = ResolutionClustering, algorithm = AlgoClustering)
