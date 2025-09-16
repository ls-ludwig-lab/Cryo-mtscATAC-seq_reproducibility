library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(readxl)
library(RColorBrewer)
library(dplyr)
library(tidyr)

# Read in th e
setwd(".../output")

df <- readRDS(file="integrated.rds")


DefaultAssay(df) <- "RNA"

# Define marker genes by cell type
markers <- list(
  Oligodendrocytes = c("MBP", "PLP1", "MOG", "CLDN11"),
  OPCs = c("PDGFRA", "CSPG4", "SOX10", "OLIG2"),
  Committed_OPCs = c("SOX10", "OLIG1", "ASPA", "TNR"),
  Astrocytes = c("GFAP", "AQP4", "SLC1A3", "ALDH1L1"),
  Excitatory_Neurons = c("SLC17A7", "CAMK2A", "NEUROD2", "SATB2"),
  Inhibitory_Neurons = c("GAD1", "GAD2", "SST", "PVALB"),
  Endothelial_Pericytes = c("CLDN5", "PECAM1", "PDGFRB", "RGS5"),
  TCells = c("CD3D", "IL7R", "FOXP3", "CD8A"),
  Microglia_TAMs = c("CX3CR1", "P2RY12", "CD68", "TMEM119"),
  GBM_Cells = c("SOX2", "EGFR", "CD44", "VIM")
)

# Flatten the list into a single vector of features
all_features <- unlist(markers)

plot_feature <- function(feature) {
  p <- FeaturePlot(object = df, 
                   features = feature,  
                   pt.size = 0.1, order = TRUE) &
    scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu"))) &
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          axis.line = element_blank())
  
  return(p)
}

# Loop through each cell type
for (cell_type in names(markers)) {
  features <- markers[[cell_type]]  # Get markers for this cell type
  
  # Generate FeaturePlots for each marker
  plots <- lapply(features, plot_feature)
  
  # Arrange them in a grid
  grid_plot <- wrap_plots(plots, ncol = 2) + plot_annotation(title = paste0(cell_type, " Marker Expression"))
  
  # Save each cell type as an SVG
  file_name <- paste0(cell_type, "_markers.pdf")
  ggsave(file_name, plot = grid_plot, width = 6, height = 6, dpi = 300)
  
  print(paste("Saved:", file_name))
}

