library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(circlize)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(readxl)
library(RColorBrewer)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(BuenColors)

setwd(".../output/")

df <- readRDS(file="preprocessed.rds")

variants1 <- readRDS("variants_ncof5.rds")

# Subset to HQ cells that exist so far
variants1 <- variants1[,colnames(df)]
misc_df <- data.frame(rowData(variants1))

assay_data <- assay(variants1, "allele_frequency")  # Assuming "counts" is the assay name
assay_data

# Add the new assay to the Seurat object
df[["alleles"]] <- CreateAssayObject(counts = assay_data)
DefaultAssay(df) <- "alleles"

misc_df1_R10 <- data.frame(rowData(variants1))
variant <- misc_df1_R10$variant

# Function to process each variant and update metadata
kept_mat <- as.matrix(assay_data[variant, ])
updated_labels <- ifelse(kept_mat > 0.2, paste(variant), "wt")
rownames(updated_labels) <- gsub(">", "_", rownames(updated_labels))
updated_labels1 <- t(updated_labels)

new_column_names <- colnames(updated_labels1)
new_column_names <- ifelse(grepl("^[0-9]", new_column_names), paste0("X", new_column_names), new_column_names)
colnames(updated_labels1) <- new_column_names

umap_coords <- Embeddings(df, "umap")
umap_df1 <- as.data.frame(umap_coords)
df@meta.data <- cbind(df@meta.data, umap_df1)
# View the updated metadata
meta.data <- df@meta.data
meta.data

R10_meta <- cbind(meta.data,updated_labels1)


mean_df <- read.csv("mean_abundance.csv")

head(mean_df)

mean_df <- mean_df %>%
  arrange(desc(primary)) %>%
  mutate(v2 = factor(v2, levels = v2))
head(mean_df)

top_v2_vector <- mean_df %>% # Order by 'primary' from highest to lowest
  slice_head(n = 50) %>%  # Select top 20 rows
  pull(v2)  # Extract the 'v2' column as a vector

# Check the vector of top 20 v2 values
print(top_v2_vector)


mutation_columns <- c("16390G>A", "2221C>T", "16391G>A", "5540G>A", "16291C>T", "13676A>G", "7362G>A", "11711G>A", "5703G>A", "15755T>C", "16292C>T", "189A>G", "709G>A", "2992G>A", "5112G>A",
                      "2702G>A", "5814T>C", "5329T>C", "564G>A", "10310G>A", "66G>A", "7950T>C", "2571G>A", "3591G>A", "16366C>T", "16261C>T", "879T>C", "385A>G", "16293A>G", "195T>C",
                      "16327C>T", "14384G>A", "368A>G", "199T>C", "3380G>A", "1485G>A", "249A>G", "16223C>T", "16166A>G", "294T>C", "2275T>C", "6943T>C", "5140G>A", "16465C>T", "3243A>G",
                      "16389G>A", "16357T>C", "16355C>T", "13042G>A", "279T>C")


# Subset the matrix to keep only the rows of interest
subset_mat <- kept_mat[rownames(kept_mat) %in% mutation_columns, ]

# Identify columns with at least one value > 0.2
cols_to_keep <- colSums(subset_mat > 0.2) > 0

# Subset the matrix to keep only relevant columns
filtered_mat <- subset_mat[, cols_to_keep]

# View result
filtered_mat

filtered_metadata <- R10_meta[rownames(R10_meta) %in% colnames(filtered_mat), ]

head(filtered_metadata)

# Ensure rownames of metadata match colnames of matrix
filtered_metadata <- filtered_metadata[rownames(filtered_metadata) %in% colnames(filtered_mat), ]
filtered_mat <- filtered_mat[, colnames(filtered_mat) %in% rownames(filtered_metadata)]

brewer_palette <- brewer.pal(5, "YlGnBu")

# Reorder metadata by Celltype
filtered_metadata <- filtered_metadata[order(filtered_metadata$peaks_snn_res.0.5), ]

# Reorder matrix columns accordingly
filtered_mat <- filtered_mat[, rownames(filtered_metadata)]

ordered_indices <- order(filtered_metadata$peaks_snn_res.0.5)  # Sort by Celltype

filtered_metadata <- filtered_metadata[ordered_indices, ]  # Keep annotation in sync

head(filtered_metadata$peaks_snn_res.0.5)
head(filtered_mat)

filtered_metadata <- filtered_metadata %>%
  mutate(Celltype = case_when(
    peaks_snn_res.0.5 == "0" ~ "SMC",
    peaks_snn_res.0.5 == "1" ~ "SMC I",
    peaks_snn_res.0.5 == "2" ~ "Fibroblasts",
    peaks_snn_res.0.5 == "3" ~ "Macrophages",
    peaks_snn_res.0.5 == "4" ~ "NK/T Cells",
    peaks_snn_res.0.5 == "5" ~ "Endothelial Cells",
    TRUE ~ "Unknown"
  ))

colnames(filtered_metadata)

head(filtered_metadata)
write.csv(filtered_metadata, "filtered_metadataTop50.csv", row.names = TRUE)

col_fun_mtDNA = colorRamp2(c(5, 50, 10), hcl_palette = "Viridis")

column_ha = HeatmapAnnotation(
  "Celltype" = filtered_metadata$Celltype,
  "mtDNA_depth" = filtered_metadata$mtDNA_depth,
  show_legend = c(FALSE),
  col = list(Celltype = c(
    "SMC" = "deepskyblue3",
    "SMC I" = "skyblue2",
    "Fibroblasts" = "darkseagreen3",
    "Macrophages" = "plum2",
    "NK/T Cells" = "lightsalmon",
    "Endothelial Cells" = "brown2"
  ),
  "mtDNA_depth" = col_fun_mtDNA))

# 2. Count number of cells where each variant is present
variant_abundance <- rowSums(filtered_mat)

# 3. Order rows by abundance (descending)
ordered_rows <- names(sort(variant_abundance, decreasing = TRUE))

# 4. Reorder matrix
filtered_mat_ordered <- filtered_mat[ordered_rows, ]


# 1. Binarize heteroplasmy matrix
# You can adjust the threshold (e.g., > 0.05)
binary_mat <- (filtered_mat > 0.05) * 1

# 2. Compute binary distance and clustering
binary_dist <- dist(t(binary_mat), method = "binary")  # Binary (Jaccard) distance
binary_clust <- hclust(binary_dist, method = "ward.D")



# 3. Heatmap with column clustering based on binary shared variant presence
p1 <- Heatmap(
  filtered_mat_ordered,
  name = "Heteroplasmy",
  top_annotation = column_ha,
  col = as.character(jdb_palette("solar_rojos", type = "continuous")),
  cluster_rows = TRUE,
  cluster_columns = binary_clust,
  show_row_names = TRUE,
  show_column_names = FALSE,
  show_column_dend = TRUE,
  show_row_dend = FALSE,
  border = TRUE,
  height = unit(200, "mm"),
  width = unit(200, "mm"),
  row_gap = unit(0, "mm"), 
  column_gap = unit(0, "mm")
)


p1

pdf(file="ComplexHeat_matrix_Top50_hclust.pdf", width = 20, height = 10)
p1
dev.off()

