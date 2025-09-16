# Load Libraries
library(Seurat)       
library(ggplot2)      
library(tidyverse)    
library(dplyr)        
library(scCustomize)  

# Load Seurat Object
seurat <- readRDS("path_to/GN01.rds")  # Load Ganglioneuroma dataset

# ============================
# Manual Annotation - Renaming Cluster Identities
# ============================
## Copy cluster IDs to a new metadata column for manual annotation
seurat$Annotation_Manual <- seurat$seurat_clusters

## Set active identities to the manual annotation column
Idents(seurat) <- "Annotation_Manual"

## Rename clusters manually based on integration with reference datasets
seurat <- RenameIdents(object = seurat,
                       "0" = "Monocyte_Macrophage",
                       "1" = "Fibroblast",
                       "2" = "Fibroblast",
                       "3" = "Schwann_cells",
                       "4" = "T_cells", 
                       "5" = "Other",
                       "6" = "B_cells",
                       "7" = "Endothelial",
                       "8" = "Fibroblast")

## Update metadata to reflect new identities
seurat$Annotation_Manual <- Idents(seurat)

# ============================
# Dotplot for Marker Gene Expression
# ============================
## Define marker genes manually for each cell type
Endothelial <- c("EMCN", "FLT1", "PTPRB", "PECAM1")
Fibroblast <- c("COL1A1", "COL1A2", "COL3A1", "PDGFRB", "DCN")
Schwann_cells <- c("SOX10", "S100B", "PLP1", "ERBB3")
Bcell <- c("CD22", "CD20", "CCR6", "FCRL1")
Tcell <- c("CD3E", "CD8A", "TIGIT")
Monocyte_Macrophage <- c("CD68", "CD163", "MRC1", "CD14")

## Use RNA assay for gene expression dotplot
DefaultAssay(seurat) <- "RNA"
Idents(seurat) <- "Annotation_Manual"

## Reorder identities in the desired order for plotting
desired_order <- c("T_cells", "B_cells","Endothelial", "Monocyte_Macrophage", "Other", "Fibroblast", "Schwann_cells")  
Idents(seurat) <- factor(Idents(seurat), levels = desired_order)

# ============================
# DotPlot Marker Genes
# ============================
DotPlot(
  seurat,  
  assay = "RNA",
  features = unique(c(Endothelial, Fibroblast, Schwann_cells, Bcell, Tcell, Monocyte_Macrophage)),
  cols = c("lightblue", "red4"),   # Gradient from low to high expression
  col.min = -2.5,                  # Minimum expression cutoff
  col.max = 2.5,                   # Maximum expression cutoff
  dot.min = 0,                      # Minimum fraction of cells expressing a gene to show a dot
  dot.scale = 8,                    # Dot size scaling
  group.by = NULL,                  # Default grouping by active identities
  split.by = NULL,                  # No splitting
  cluster.idents = FALSE,           # Do not cluster cell types on y-axis
  scale = TRUE,                     # Scale gene expression across features
  scale.by = "radius"               # Scale dot size by fraction of expressing cells
) + theme(axis.text.x = element_text(color="black", hjust = 1, angle = 45,
                                     size=14),
          axis.text.y = element_text(color="black",
                                     size=14))
