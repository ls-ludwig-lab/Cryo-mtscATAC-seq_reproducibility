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


setwd(".../integration/output")

d <- readRDS(file="integrated_annotated.rds")
table(d$dataset)
meta.data <- d@meta.data

# List of y variables
y_vars <- c("nCount_ATAC", "mtDNA_depth", "TSS.enrichment")  # Replace with actual variable names

dataset_color <- c("FF_10min" = "palevioletred2",
                   "FF30_min" = "darkorange",
                   "Cells" = "darkslategray3" )


cell_type_order <- c(
  # Myeloid lineage
  "Myeloid Cell",
  "Monocyte",                      
  "Monocyte/ Macrophage",         
  "Macrophage",                   
  "Dendritic Cell",               
  "Plasmacytoid Dendritic Cell",  
  "Neutrophil / Monocytes",       
  
  # Erythroid lineage
  "Erythroid Progenitor Cell",               
  
  # Lymphoid lineage
  "NK Cell",                      
  "NK/ T Cell",                   
  "CD8 T Cell",                   
  "CD4 T Cell",                   
  
  # B cell lineage
  "B Cell",                       
  "Marginal Zone B Cell",        
  
  # Unresolved clusters
  "other"
)

# Set the factor levels in the biologically-informed order
meta.data$possible_Celltype <- factor(meta.data$possible_Celltype, levels = cell_type_order)
cell_type_colors <- c(              # Blue
  "Myeloid Cell" = "lightpink2",      # Orange
  "CD8 T Cell" = "dodgerblue3",                  # Green
  "Monocyte/ Macrophage" = "plum2",        # Red
  "Marginal Zone B Cell" = "darkseagreen3",        # Purple
  "Erythroid Progenitor Cell" = "coral2",                     # Pink
  "other" = "wheat2",                  # Gray
  "Macrophage" = "plum3",                  # Olive
  "unresolved II" = "navajowhite",               # Cyan
  "B Cell" = "darkseagreen2",                      # Light Blue
  "unresolved I" = "wheat3",                # Peach
  "NK/ T Cell" = "dodgerblue",                  # Light Green
  "NK Cell" = "cornflowerblue",                     # Light Red
  "CD4 T Cell" = "deepskyblue3",                  # Lavender
  "Dendritic Cell" = "lightpink", # Light Pink
  "Monocyte" = "plum"                     # Light Brown
)


# Loop through y_vars to create and save plots
for (y_var in y_vars) {
  # Create the plot for the current y_var
  p <- ggplot(meta.data, aes_string(x = "dataset", y = y_var, fill = "dataset")) +
    geom_boxplot(trim = TRUE, color = "grey50") +
    theme_classic() +
    scale_fill_manual(values=dataset_color) +
    scale_y_log10() +  # Apply logarithmic scale to y-axis
    theme(
      axis.title.y = element_text(size = 20),
      axis.line.x = element_line(),
      axis.text.y = element_text(size = 15),
      axis.title.x = element_blank()
    )
  print(p)
  # Save the plot
  ggsave(paste0("boxplot_", y_var, ".pdf"), plot = p, width = 7, height = 5)
}


# Loop through y_vars to create and save plots
for (y_var in y_vars) {
  # Create the plot for the current y_var
  p <- ggplot(meta.data, aes_string(x = "possible_Celltype", y = y_var, fill = "possible_Celltype")) +
    geom_boxplot(trim = TRUE, color = "grey50") +
    theme_classic() +
    scale_fill_manual(values=cell_type_colors) +
    scale_y_log10() +  # Apply logarithmic scale to y-axis
    theme(
      axis.title.y = element_text(size = 20),
      axis.line.x = element_line(),
      axis.text.y = element_text(size = 15),
      axis.title.x = element_blank()
    )
  print(p)
  # Save the plot
  ggsave(paste0("boxplot_", y_var, ".pdf"), plot = p, width = 7, height = 5)
}

