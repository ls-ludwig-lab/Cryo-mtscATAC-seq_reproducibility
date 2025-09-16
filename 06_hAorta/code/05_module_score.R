library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(RColorBrewer)

df <- readRDS(file="preprocessed.rds")

DefaultAssay(df) <- "peaks"
 
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


DefaultAssay(df) <- 'RNA'

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
markers <- FindAllMarkers(df, only.pos = TRUE)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

feature1 <- read_excel(".../Mosquera_etal.xlsx")
feature1

unique_values <- unique(feature1$cluster)
# Create an empty list to store the individual data frames
separated_dfs <- list()

# Iterate over each unique value
for (value in unique_values) {
  # Filter the data frame based on the current value
  separated_dfs[[value]] <- subset(feature1, cluster == value)
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
  cell_marker_gene_list <- separated_dfs[[value]]$gene
  DefaultAssay(df) <- "RNA"
  df <<- AddModuleScore(object = df, features = list(cell_marker_gene_list), name = paste0(value, "_Score"))
}

# Apply the function to each unique value using lapply
lapply(unique_values, add_module_score)
head(df@meta.data)


meta_columns <- c("B.cell_Score1", "Endothelial_Score1", "Fibroblast_Score1", "SMC_Score1", "Mast.cell_Score1", "Neuron_Score1", "Pericyte_Score1", "Macrophage_Score1", "Plasma.cell_Score1", "pDC_Score1", "T.NK_Score1")
# Loop through the metadata columns and create DimPlots
for (meta in meta_columns) {
  # Generate the DimPlot
  p <- FeaturePlot(object = df, features = meta,  pt.size =0.5,  order = TRUE)  & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  p + theme(axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks=element_blank(),
            #legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank(),
            axis.line = element_blank())
  
  print(p)
  # Save the plot
  ggsave(filename = paste0("Feature_", meta, ".svg"), plot = p, width = 6, height = 6, dpi=300)
}

