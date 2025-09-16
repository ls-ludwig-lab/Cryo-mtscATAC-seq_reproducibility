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



setwd("/data/cephfs-1/work/projects/ludwig-mtscatac-als/ALS/2024_Region/test_output")

df <- readRDS(file="merged.rds")

df <- subset(df,
              subset = mtDNA_depth < 100)

unique(df@meta.data$dataset)

df3 <- subset(df, subset = dataset == "Region_5")
df4 <- subset(df, subset = dataset == "Region_10")
df5 <- subset(df, subset = dataset == "Region_28")
df6 <- subset(df, subset = dataset == "Region_29")


# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(df3, df4, df5, df6),
  anchor.features = rownames(df3),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = df[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
DefaultAssay(integrated) <- "ATAC"
integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
integrated <- FindClusters(object = integrated, resolution = 0.6, verbose = FALSE, algorithm = 3)

p1 <- DimPlot(integrated, group.by="dataset", label.size = 4, label=FALSE) +
  # -- Here you tell ggplot to eliminate all the axis info from your original plot -- #
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        #legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) 
p1

ggsave(plot = p1, filename = "DimPlot_Integrated_dataset.png", width = 7, height = 7, dpi=300)



#proceeded with Harmony integration and used this data rds file for downstream analysis
hm.integrated <- RunHarmony(
  object = df,
  group.by.vars = 'dataset',
  reduction.use = 'lsi',
  nclust = 8,
  assay.use = 'ATAC',
  project.dim = FALSE
)


# re-compute the UMAP using corrected LSI embeddings
hm.integrated <- FindNeighbors(object = hm.integrated, reduction = 'harmony', dims = 2:30)
hm.integrated <- FindClusters(object = hm.integrated, resolution = 1.0, verbose = FALSE, algorithm = 3)

hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony')

p5 <- DimPlot(hm.integrated, group.by = 'dataset', pt.size = 0.1) + ggplot2::ggtitle("Harmony integration")
p5

feature <- as.data.frame(read_excel(".../marker_genes.xlsx"))

unique_values <- unique(feature$Cluster)

# Create an empty list to store the individual data frames
separated_dfs <- list()


gene.activities <- GeneActivity(hm.integrated)

# add to the Seurat object as a new assay
hm.integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)

hm.integrated <- NormalizeData(
  object = hm.integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(hm.integrated$nCount_RNA)
)


# Iterate over each unique value
for (value in unique_values) {
  # Filter the data frame based on the current value
  separated_dfs[[value]] <- subset(feature, Cluster == value)
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
  DefaultAssay(hm.integrated) <- "RNA"
  hm.integrated <<- AddModuleScore(object = hm.integrated, features = list(cell_marker_gene_list), name = paste0(value, "_Score"))
}

# Apply the function to each unique value using lapply
lapply(unique_values, add_module_score)
head(hm.integrated@meta.data)

meta_columns <- c("Oligodendrocytes_Score1",
                  "OPCs_Score1", 
                  "Astrocytes_Score1", 
                  "Neuronal.Cells_Score1",
                  "Microglia_Score1")
# Loop through the metadata columns and create DimPlots
for (meta in meta_columns) {
  # Generate the DimPlot
  p <- FeaturePlot(object = hm.integrated, features = meta,  pt.size =0.5,  order = TRUE)  & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  print(p)
  # Save the plot
  ggsave(filename = paste0("Feature_", meta, ".png"), plot = p, width = 6, height = 6, dpi=300)
  
}

saveRDS(hm.integrated, "hm.integrated.rds")
