library(Signac)
library(Seurat)
library(scCustomize)
library(ggplot2)
library(patchwork)
library(readxl)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(readr)


setwd("../output")

df <- readRDS(file="integrated.rds")

feature <- as.data.frame(read_excel(".../BlancoCarmona_modules.xlsx"))

feature <- feature %>%
  pivot_longer(
    cols = everything(),   # Pivot all columns including BAMs
    names_to = "Module",   # Column names go to 'Module'
    values_to = "Gene"    # Values go to 'Genes'
  ) %>%
  drop_na(Gene) 


unique_values <- unique(feature$Module)
# Create an empty list to store the individual data frames
separated_dfs <- list()

# Iterate over each unique value
for (value in unique_values) {
  # Filter the data frame based on the current value
  separated_dfs[[value]] <- subset(feature, Module == value)
}
# Access individual data frames
for (value in unique_values) {
  cat(paste("Data frame for '", value, "':\n", sep = ""))
  print(separated_dfs[[value]])
  cat("\n")
}

#The add_module_score function takes a module value, extracts the gene list, and uses AddModuleScore to add the module score to the Seurat object.
#lapply is used to apply this function to each unique module value.
#The <<- operator is used to modify the df object within the function.

# Function to add module scores to the Seurat object
add_module_score <- function(value) {
  cell_marker_gene_list <- separated_dfs[[value]]$Gene
  DefaultAssay(df) <- "RNA"
  df <<- AddModuleScore(object = df, features = list(cell_marker_gene_list), name = paste0(value, "_Score"))
}

# Apply the function to each unique value using lapply
lapply(unique_values, add_module_score)
head(df@meta.data)


meta_columns <- c("Astro.like_Score1",
"Astrocytes_Score1", 
"Cycling_Score1", 
"Endothelial_Score1",  
"Gradient_Score1" , 
"Microglia_Score1",
"Mixed_Score1",
"Neurons_Score1",
"Oligodendrocytes_Score1", 
"OPC.like_Score1",
"Pericytes_Score1",
"RE_Score1",
"T.cells_Score1")
# Loop through the metadata columns and create DimPlots


# Create DotPlot using the module score columns
dot_plot <- DotPlot_scCustom(
  seurat_object = df,
  features = meta_columns,
  group.by = "seurat_clusters"
) +
  scale_color_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu"))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(x = NULL, y = "Cluster")
dot_plot
ggsave(filename = "Blanco_Camona_dotPlot.pdf", plot = dot_plot, width = 8, height = 8)


# Neftle et al.
feature <- read_excel(".../meta_module_list.xlsx")

unique_values <- unique(feature$Module)
# Create an empty list to store the individual data frames
separated_dfs <- list()

# Iterate over each unique value
for (value in unique_values) {
  # Filter the data frame based on the current value
  separated_dfs[[value]] <- subset(feature, Module == value)
}
# Access individual data frames
for (value in unique_values) {
  cat(paste("Data frame for '", value, "':\n", sep = ""))
  print(separated_dfs[[value]])
  cat("\n")
}

#The add_module_score function takes a module value, extracts the gene list, and uses AddModuleScore to add the module score to the Seurat object.
#lapply is used to apply this function to each unique module value.
#The <<- operator is used to modify the df object within the function.

# Function to add module scores to the Seurat object
add_module_score <- function(value) {
  cell_marker_gene_list <- separated_dfs[[value]]$Gene
  DefaultAssay(df) <- "RNA"
  df <<- AddModuleScore(object = df, features = list(cell_marker_gene_list), name = paste0(value, "_Score"))
}

# Apply the function to each unique value using lapply
lapply(unique_values, add_module_score)
head(df@meta.data)


meta_columns <- c("MES2_Score1", "MES1_Score1", "AC_Score1", "OPC_Score1", "NPC1_Score1", "NPC2_Score1", "G1_S_Score1", "G2_M_Score1")
# Create DotPlot using the module score columns
dot_plot <- DotPlot_scCustom(
  seurat_object = df,
  features = meta_columns,
  group.by = "seurat_clusters"
) +
  scale_color_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu"))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(x = NULL, y = "Cluster")
dot_plot
ggsave(filename = "Neftel_etal_dotplot.pdf", plot = dot_plot, width = 8, height = 8)




