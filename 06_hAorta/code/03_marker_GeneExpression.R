library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)


df <- readRDS(file="preprocessed.rds")

DefaultAssay(df) <- "peaks"
 
df <- FindClusters(df, resolution = 0.7, algorithm = 3)

Idents(df) <-"peaks_snn_res.0.5"

p1 <- DimPlot(df, label.size = 4, label=TRUE) +
  # -- Here you tell ggplot to eliminate all the axis info from your original plot -- #
  theme(aspect.ratio=1/1,
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),#legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
p1
ggsave(plot = p1, filename = "UMAP_peaks_snn_res.0.5.pdf", width = 6, height = 5, dpi=300)



DefaultAssay(df) <- 'RNA'

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
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


