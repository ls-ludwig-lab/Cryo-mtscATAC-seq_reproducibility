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


meta.data <- meta.data %>%
  mutate(Celltype = case_when(
    peaks_snn_res.0.5 == "0" ~ "SMC",
    peaks_snn_res.0.5 == "1" ~ "SMC I",
    peaks_snn_res.0.5 == "2" ~ "Fibroblasts",
    peaks_snn_res.0.5 == "3" ~ "Macrophages",
    peaks_snn_res.0.5 == "4" ~ "NK/T Cells",
    peaks_snn_res.0.5 == "5" ~ "Endothelial Cells",
    TRUE ~ "Unknown"
  ))


R10_meta <- cbind(meta.data,updated_labels1)

mean_df <- read.csv("mean_df_variants.csv")

top_v2_vector <- mean_df %>%
  arrange(desc(primary)) %>%  # Order by 'primary' from highest to lowest
  slice_head(n = 50) %>%  # Select top 20 rows
  pull(v2)  # Extract the 'v2' column as a vector

# Check the vector of top 20 v2 values
print(top_v2_vector)

mutation_columns <- c( "16390G>A", "2221C>T" , "16391G>A", "5540G>A",  "16291C>T", "13676A>G", "7362G>A",  "11711G>A", "5703G>A" , "15755T>C", "5979G>A", 
                       "16292C>T", "189A>G" ,  "709G>A" ,  "2992G>A"  ,"5112G>A"  ,"2702G>A",  "5814T>C" , "5329T>C" , "8141G>A"  ,"564G>A"   ,"10310G>A",
                       "66G>A"  ,  "7950T>C" , "2571G>A",  "5553T>C"  ,"3591G>A",  "16366C>T" ,"16261C>T" ,"879T>C" ,  "385A>G" ,  "16293A>G" ,"195T>C"  ,
                       "10373G>A" ,"16327C>T" ,"12439T>C", "14384G>A" ,"368A>G" ,  "199T>C"  , "10290G>A", "3380G>A" , "1485G>A" , "249A>G"  , "15729T>C",
                       "12093T>C", "16223C>T", "16166A>G", "294T>C" ,  "2275T>C" , "6943T>C")


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

# Reorder metadata by Celltype
filtered_metadata <- filtered_metadata[order(filtered_metadata$peaks_snn_res.0.5), ]

# Reorder matrix columns accordingly
filtered_mat <- filtered_mat[, rownames(filtered_metadata)]

ordered_indices <- order(filtered_metadata$peaks_snn_res.0.5)  # Sort by Celltype

filtered_metadata <- filtered_metadata[ordered_indices, ]  # Keep annotation in sync

head(filtered_metadata$peaks_snn_res.0.5)
colnames(filtered_metadata)
head(filtered_metadata)




#: Define a function to compute proportion dataframe
get_celltype_prop_df <- function(meta.data, group.label) {
  ab <- table(meta.data$Celltype, meta.data$orig.ident)
  prop <- sweep(x = ab, MARGIN = 2, STATS = colSums(ab), FUN ="/")
  nb.samples <- ncol(ab)
  nb.celltype <- nrow(ab)
  
  celltype.name.vec <- rep(rownames(ab), nb.samples)
  sample.name.vec <- as.vector(sapply(X = colnames(ab), FUN = function(x) rep(x, nb.celltype)))
  prop.cells.vec <- as.vector(prop)
  
  df.plot <- data.frame(
    CellType = celltype.name.vec,
    Sample = sample.name.vec,
    CellProp = prop.cells.vec,
    Group = group.label
  )
  return(df.plot)
}

#Run for both full and subset data
df.full <- get_celltype_prop_df(meta.data, "Full")
df.sub <- get_celltype_prop_df(filtered_metadata, "Subset")

df.combined <- rbind(df.full, df.sub)


# Make a grouped bar plot
df.combined$CellType <- factor(df.combined$CellType, levels = celltype_order)

p1 <- ggplot(data = df.combined, aes(x = Group, y = CellProp, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", color = "gray", linewidth = 0.1) +
  scale_fill_manual(values = colors.Celltype) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  ylab("Proportion of Cells")

colors.Celltype <- c(
  "SMC" = "deepskyblue3",
    "SMC I" = "skyblue2",
    "Fibroblasts" = "darkseagreen3",
  "Macrophages" = "plum2",
  "NK/T Cells" = "lightsalmon",
  "Endothelial Cells" = "brown2")

p1 <- ggplot(df.combined, aes(x = Group, stratum = CellType, alluvium = CellType, y = CellProp, fill = CellType)) +
  geom_flow(width = 0.5) +
  scale_fill_manual(values = c(colors.Celltype)) +
  geom_stratum(width = 0.6, color = "grey80") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Cell Proportion (%)", x = "", fill = "Cell Type") +
  theme(
    #axis.line=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position="none",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank())
p1
ggsave(plot = p1, filename = "PropCells_Allu_Subset_Top50.pdf", width = 4, height = 10, dpi=300)

# Statistical Test

df.combined %>%
  group_by(CellType) %>%
  summarise(p_value = wilcox.test(CellProp ~ Group)$p.value)

