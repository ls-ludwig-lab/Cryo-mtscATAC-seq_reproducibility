library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(readxl)
library(RColorBrewer)
library(SummarizedExperiment)



setwd("...")

d <- readRDS(file="integrated_annotated.rds")

#change to the Region of Interest
df <- subset(d, subset = dataset == "Region_5")

# Check the current cell names
head(Cells(df))
# Modify the cell names to remove the  prefix
new_cell_names <- gsub(pattern = "^Region_5_", replacement = "", x = Cells(df))
# Set the modified cell names back to the Seurat object
df <- RenameCells(object = df, new.names = new_cell_names)
# Check the modified cell names
head(Cells(df))

variantsR5 <- readRDS(".../Region/R5.rds")

#adding the variants as assay data to the seurat objects
variantsR5 <- variantsR5[,colnames(df)]

assay_data <- assay(variantsR5, "allele_frequency") 
assay_data
# Add the new assay to the Seurat object
df[["alleles"]] <- CreateAssayObject(counts = assay_data)
DefaultAssay(df) <- "alleles"
misc_df_R10 <- data.frame(rowData(variantsR5))
variant <- misc_df_R10$variant

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

R5_meta <- cbind(meta.data,updated_labels1)

colnames <- as.data.frame(colnames(R5_meta))

# Add a new column with the desired transformation
colnames <- colnames %>%
  mutate(
    formatted_variants = gsub("^X", "", gsub("_", ">", `colnames(R5_meta)`))
  )

variants <- colnames %>% dplyr::filter(formatted_variants %in% variant)
head(variants)
# Rename a column
variants <- variants %>%
  dplyr::rename(column = "colnames(R5_meta)")

variants <- variants %>%
  dplyr::rename(value = formatted_variants)


write.csv(variants, "variants_to_plot.csv")

to_subset <- c("16147C>T", "1206G>A", "1345G>A")

variants1 <- variants %>%
  dplyr::filter(value %in% to_subset)


# Set default assay
DefaultAssay(df) <- "alleles"


# Loop through each variant
for (i in 1:nrow(variants1)) {
  column <- variants1$column[i]
  value <- variants1$value[i]
  
  # Create the first plot (p1)
  p1 <- ggplot() +
    geom_point(
      data = R5_meta,
      aes(x = UMAP_1, y = UMAP_2, color = "orig.ident"),
      size = 1, alpha = 0.8
    ) +
    geom_point(
      data = subset(R5_meta, get(column) == value),
      aes(x = UMAP_1, y = UMAP_2, color = value),
      size = 1
    ) +
    scale_color_manual(values = c("firebrick", "grey")) +   
    theme_void() + 
    theme(legend.position = "none")
  
  # Define filenames
  p1_filename <- paste0("Region5_", value, "_p1.pdf")
  
  # Save plots separately
  ggsave(plot = p1, filename = p1_filename, width = 5, height = 5, dpi = 300, limitsize = FALSE)
}

