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


setwd("...")

df <- readRDS(file="hm.integrated.rds")
p <- DimPlot(df, label = TRUE)
p


Idents(df)
DefaultAssay(df) <- "RNA"
marker<- FindAllMarkers(df,  min.pct = 0.25)

# Identify unique values
unique_values <- unique(marker$cluster)

# Create an empty list to store the individual data frames
separated_dfs <- list()

# Iterate over each unique value
for (value in unique_values) {
  # Filter the data frame based on the current value
  separated_dfs[[value]] <- subset(marker, cluster == value)
}

# Access individual data frames
for (value in unique_values) {
  cat(paste("Data frame for '", value, "':\n", sep = ""))
  print(separated_dfs[[value]])
  cat("\n")
}

# Define the path where CSV files will be saved
output_directory <- "../output/"

# Ensure the output directory exists
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# Save individual data frames as CSV files
for (value in unique_values) {
  cat(paste("Data frame for '", value, "':\n", sep = ""))
  print(separated_dfs[[value]])
  cat("\n")
  
  # Define the file name for the CSV
  file_name <- file.path(output_directory, paste("data_frame_", value, ".csv", sep = ""))
  
  # Save the data frame as a CSV
  write.csv(separated_dfs[[value]], file_name, row.names = FALSE)
  cat(paste("Data frame for '", value, "' saved as '", file_name, "'.\n", sep = ""))
}

DimPlot(df, label = TRUE)
 
df <- RenameIdents(df,"0" = "Oligodendrocytes",
                   "1" = "Oligodendrocytes I",
                   "2" = "Microglia",
                   "3" = "Glia Cells",
                   "4" = "Astrocytes",
                   
                   "5" = "Glia Cells I",
                   "6" = "Glia Cells II",
                   
                   "7" = "exitatory neurons",
                   
                   "8" = "Glia Cells III",
                   
                   "9" = "GABAergic interneurons",
                   "10" = "excitatory neurons I",
                   
                   "11" = "Endothelial Cells",
                   
                   "12" = "OPCs",
                   "13" = "OPCs I",
                   "14" = "OPCs II",
                   "15" = "proinflammatory Microglia",
                   
                   "16" = "CD8 T-Cells",
                   "17" = "Astrocytes I")
DimPlot(df, label = TRUE)


Idents(df)
# Save cluster identities as a new column in meta.data
df$possible_Celltype <- df@active.ident

saveRDS(df, "integrated_annotated_1.rds")
