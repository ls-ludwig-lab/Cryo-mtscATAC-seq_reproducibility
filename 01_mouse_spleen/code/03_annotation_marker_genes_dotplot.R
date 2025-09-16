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

setwd(".../integration/output")

df <- readRDS(file="hm.integrated.rds")

meta_columns <- c("Marginal_zone_B_cell_Score1",
                  "T_cell_Score1", 
                  "NK_cell_Score1", 
                  "Dendritic.cell_S100a4_high_Score1",
                  "Monocyte_Score1",
                  "Neutrophil_S100a8_high_Score1",
                  "Macrophage_Score1",
                  "B_cell_Score1",
                  "Plasmacytoid_dendritic_cell_Siglech_high_Score1",
                  "Erythroid_cell_Score1",
                  "Neutrophil_Ngp_high_Score1")

# Loop through the metadata columns and create DimPlots
for (meta in meta_columns) {
  # Generate the DimPlot
  p <- FeaturePlot(object = df, features = meta,  pt.size =0.5,  order = TRUE)  & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  # Display the plot in the plot panel
  print(p)
  
  # Save the plot
  ggsave(filename = paste0("Feature_", meta, ".png"), plot = p, width = 6, height = 6, dpi=300)
  
}

feature <- as.data.frame(read_excel(".../Spleen_CellMarkers_mouse.xlsx"))
feature

unique_values <- unique(feature$alias)
unique_values

library(dplyr)

get_top_unique_genes <- function(feature, n = 5) {
  feature <- feature %>%
    arrange(p_val) %>%         # sort by p-value globally
    group_by(alias) %>%
    mutate(rank = row_number()) %>%
    ungroup()
  
  selected <- c()
  output <- list()
  
  for (al in unique(feature$alias)) {
    alias_genes <- feature %>%
      filter(alias == al) %>%
      pull(gene)
    
    # take genes not yet selected
    unique_genes <- alias_genes[!alias_genes %in% selected]
    
    # if fewer than n left, just take what’s there
    chosen <- head(unique_genes, n)
    
    selected <- c(selected, chosen)  # update selected list
    output[[al]] <- chosen
  }
  
  output
}

# Example usage
top_unique <- get_top_unique_genes(feature, n = 5)
top_unique


# 0) Choose the assay you actually want to plot genes from
#    (RNA for gene expression, or 'GeneActivity' if using gene activity)
DefaultAssay(df) <- "RNA"   # change if needed

# 1) Make a single unique vector of genes from your list
#    (you currently have `top_unique`, not `top_genes_per_alias`)
features_vec <- unique(unlist(top_unique))

# 2) Match to what exists in the object (row names of the chosen assay)
feats_in_obj <- rownames(GetAssayData(df, assay = DefaultAssay(df), slot = "data"))

# If your object stores Ensembl IDs, map first; otherwise just intersect:
features_found <- intersect(features_vec, feats_in_obj)
features_missing <- setdiff(features_vec, feats_in_obj)

if (length(features_missing) > 0) {
  message("Skipping missing features: ", paste(features_missing, collapse = ", "))
}

# 3) Ensure the plotting vector has no duplicates (prevents factor-level error)
features_found <- unique(features_found)

manual_genes <- c("Cd8a","Cd4")
features_found <- unique(c(features_found, manual_genes))


# 4) Plot
dot_plot <- DotPlot_scCustom(
  seurat_object = df,
  features = features_found,
  group.by = "seurat_clusters"
) +
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(x = NULL, y = "Cluster")

dot_plot
DimPlot(df, group.by = "seurat_clusters", label = TRUE)


Idents(df) <- "seurat_clusters"
# wilcox is the default option for test.use
cluster13 <- FindMarkers(
  object = df,
  ident.1 = 13,
  test.use = 'wilcox',
  min.pct = 0.1
)
cluster13

df <- RenameIdents(df,"0" = "Marginal Zone B Cell",
                   "1" = "Marginal Zone B Cell",
                   "2" = "Erythroid Progenitor Cell",
                   "3" = "CD4 T Cell",
                   "4" = "Marginal Zone B Cell",
                   "5" = "Marginal Zone B Cell",
                   "6" = "Marginal Zone B Cell",
                   "7" = "Erythroid Progenitor Cell",
                   "8" = "B Cell",
                   "9" = "Erythroid Progenitor Cell",
                   "10" = "Macrophage",
                   "11" = "CD8 T Cell",
                   "12" = "CD8 T Cell",
                   "13" = "other",
                   "14" = "CD4 T Cell",
                   "15" = "other",
                   "16" = "Myeloid Cell",
                   "17" = "Dendritic Cell",
                   "18" = "other",
                   "19" = "Marginal Zone B Cell",
                   "20" = "Erythroid Progenitor Cell",
                   "21" = "Dendritic Cell",
                   "22" = "NK Cell",
                   "23" = "NK/ T Cell",
                   "24" = "Marginal Zone B Cell",
                   "25" = "Myeloid Cell",
                   "26" = "Dendritic Cell",
                   "27" = "Marginal Zone B Cell",
                   "28" = "Myeloid Cell",
                   "29" = "Erythroid Progenitor Cell",
                   "30" = "Myeloid Cell")

# Save cluster identities as a new column in meta.data
df$possible_Celltype <- df@active.ident
unique(df@meta.data$possible_Celltype)

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
df$possible_Celltype <- factor(df$possible_Celltype, levels = cell_type_order)


dot_plot <- DotPlot_scCustom(
  seurat_object = df,
  features = features_found,
  group.by = "possible_Celltype"
) +
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(x = NULL, y = "Cluster")

dot_plot
ggsave(filename = "Marker_genes.pdf", plot = dot_plot, width = 14, height = 8)
saveRDS(df, "integrated_annotated.rds")
write.csv(meta.data, "meta.data_annotated.csv")

