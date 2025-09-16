library(Seurat)
library(Signac)           
library(SummarizedExperiment)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tibble)


setwd(".../output/")

df <- readRDS(file="preprocessed.rds")
df
variants <- readRDS("variants_ncof5.rds")

#adding the variants as assay data to the seurat objects
variants <- variants[,colnames(df)]

assay_data <- assay(variants, "allele_frequency") 
assay_data
# Add the new assay to the Seurat object
df[["alleles"]] <- CreateAssayObject(counts = assay_data)
DefaultAssay(df) <- "alleles"
misc_df <- data.frame(rowData(variants))
variant <- misc_df$variant

# Function to process each variant and update metadata
kept_mat <- as.matrix(assay_data[variant, ])
updated_labels <- ifelse(kept_mat > 0.20, paste(variant), "wt")
rownames(updated_labels) <- gsub(">", "_", rownames(updated_labels))
updated_labels1 <- t(updated_labels)

new_column_names <- colnames(updated_labels1)
new_column_names <- ifelse(grepl("^[0-9]", new_column_names), paste0("X", new_column_names), new_column_names)
colnames(updated_labels1) <- new_column_names
umap_coords <- Embeddings(df, "umap")
umap_df <- as.data.frame(umap_coords)
df@meta.data <- cbind(df@meta.data, umap_df)
# View the updated metadata
meta.data <- df@meta.data
meta.data

meta.data <- cbind(meta.data,updated_labels1)

#variants manual curated from other scripts
mutation_columns <- c("X16390G_A", "X16292C_T", "X10310G_A")
mutation_columns                      


# Loop through each mutation column
for (column in mutation_columns) {
  
  # Check column exists in meta.data
  if (!(column %in% colnames(meta.data))) {
    warning(paste("Column", column, "not found in meta.data. Skipping."))
    next
  }
  
  # Subset cells with mutation (anything != "wt" and not NA)
  mutated_cells <- meta.data[meta.data[[column]] != "wt" & !is.na(meta.data[[column]]), ]
  
  # Create the plot
  p1 <- ggplot() +
    geom_point(
      data = meta.data,
      aes(x = UMAP_1, y = UMAP_2),
      color = "grey", size = 1.5, alpha = 0.6
    ) +
    geom_point(
      data = mutated_cells,
      aes(x = UMAP_1, y = UMAP_2),
      color = "firebrick", size = 1.5
    ) +
    ggtitle(column) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 10)
    )
  
  # Save the plot
  filename <- paste0("UMAP_", column, ".pdf")
  ggsave(filename, plot = p1, width = 5, height = 5, dpi = 300)
}

