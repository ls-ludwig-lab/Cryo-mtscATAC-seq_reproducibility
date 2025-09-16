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

setwd(".../output")

df <- readRDS(file="integrated.rds")
df
p1 <- DimPlot(df, label = TRUE)
p1
ggsave(filename = "DimPlot_Cluster_resolution.pdf", plot = p1, width = 7, height = 6, dpi=300)

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
output_directory <- ".../output/"

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

