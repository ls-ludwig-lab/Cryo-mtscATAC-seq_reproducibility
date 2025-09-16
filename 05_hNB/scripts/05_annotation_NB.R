# Load Libraries
library(Seurat)       
library(ggplot2)      
library(tidyverse)    
library(dplyr)        
library(scCustomize)  

# Load Seurat Object
seurat <- readRDS("path_to/Harmony_Integrated_Object.rds")  # Neuroblastoma dataset

# ============================
# Manual Annotation - Renaming Cluster Identities
# ============================
## Copy cluster IDs to a new metadata column for manual annotation
seurat$Manual_Annotation <- seurat$seurat_clusters

## Set active identities to the manual annotation column
Idents(seurat) <- "Manual_Annotation"

## Rename clusters manually based on reference integration
seurat <- RenameIdents(object = seurat,
                       "0" = "Neuroendocrine",
                       "1" = "Neuroendocrine",
                       "2" = "Neuroendocrine",
                       "3" = "Neuroendocrine",
                       "4" = "Endothelial", 
                       "5" = "Fibroblast",
                       "6" = "Neuroendocrine",
                       "7" = "Immune")

## Update metadata to reflect new identities
seurat$Manual_Annotation <- Idents(seurat)

# ============================
# Define Marker Genes
# ============================
## Manual list of marker genes for each cell type
Neuroendocrine <- c("PHOX2B", "NXPH1", "CHGA", "DBH", "ISL1", "MYCN", "CHGB", "TH", "HAND2", "PHOX2A")
Immune <- c('PTPRC', 'MS4A7', 'CD68', "IL1B", "CD14", "FCN1")
Endothelial <- c("EGFL7", "EMCN", "PLVAP", "FLT1", "PTPRB", "PECAM1")
Fibroblast <- c("COL1A1", "COL1A2", "COL3A1", "PDGFRB", "DCN")

# Set RNA as default assay for gene expression
DefaultAssay(seurat) <- "RNA"
Idents(seurat) <- "Manual_Annotation"

# ============================
# DotPlot Marker Genes for Neuroblastoma
# ============================
DotPlot(
  seurat,
  assay = "RNA",
  features = unique(c(Neuroendocrine, Immune, Endothelial, Fibroblast)),  # Combine all marker genes
  cols = c("lightblue", "red4"),        # Color gradient from low to high expression
  col.min = -2.5,                       # Minimum scaled expression
  col.max = 2.5,                        # Maximum scaled expression
  dot.min = 0,                           # Minimum fraction of cells expressing gene to display dot
  dot.scale = 8,                         # Scaling factor for dot size
  group.by = NULL,                       # Group by active identities
  split.by = NULL,                       # No splitting
  cluster.idents = FALSE,                # Do not reorder cell types
  scale = TRUE,                          # Scale expression per gene
  scale.by = "radius"                     # Dot radius reflects fraction of expressing cells
) + theme(axis.text.x = element_text(color="black", hjust = 1, angle = 45, size=14),
          axis.text.y = element_text(color="black", size=14))
