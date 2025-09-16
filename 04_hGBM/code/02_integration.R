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


setwd("/../integration/output")

df <- readRDS(file="merged.rds")

unique(df@meta.data$dataset)
# "patient1_primary_single" "patient2_recu"        "patient1_recu"        "patient2_primary"       

#An object of class Seurat 
#183609 features across 30936 samples within 1 assay 
#Active assay: ATAC (183609 features, 183609 variable features)
#2 dimensional reductions calculated: lsi, umap

df@meta.data <- df@meta.data %>%
  mutate(dataset = ifelse(dataset == "patient1_primary_single", "patient1_primary", dataset))
#[1] "patient1_primary" "patient2_recu" "patient1_recu" "patient2_primary"

#subset mtDNA depth to 

df <- subset(df,
              subset = mtDNA_depth < 100)

#An object of class Seurat 
#183609 features across 30431 samples within 1 assay 
#Active assay: ATAC (183609 features, 183609 variable features)
#2 dimensional reductions calculated: lsi, umap

df1 <- subset(df, subset = dataset == "patient1_primary")
#187802 features across 5703 samples within 1 assay 
#Active assay: ATAC (187802 features, 187799 variable features)

df2 <- subset(df, subset = dataset == "patient2_recu")
#187802 features across 7497 samples within 1 assay 
#Active assay: ATAC (187802 features, 187799 variable features)

df3 <- subset(df, subset = dataset == "patient1_recu")
#187802 features across 4020 samples within 1 assay 
#Active assay: ATAC (187802 features, 187799 variable features)

df4 <- subset(df, subset = dataset == "patient2_primary")
#187802 features across 4705 samples within 1 assay 
#Active assay: ATAC (187802 features, 187799 variable features)

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(df1, df2, df3, df4),
  anchor.features = rownames(df1),
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
integrated <- FindClusters(object = integrated, resolution = 0.8, verbose = FALSE, algorithm = 3)

p1 <- DimPlot(integrated, label.size = 4, label=FALSE) +
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

ggsave(plot = p1, filename = "DimPlot_Integrated.png", width = 7, height = 7, dpi=300)


gene.activities <- GeneActivity(integrated)

# add to the Seurat object as a new assay
integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)

integrated <- NormalizeData(
  object = integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated$nCount_RNA)
)


feature <- read_excel("/../output/meta_module_list.xlsx")

unique_values <- unique(feature$Module)
# Create an empty list to store the individual data frames
separated_dfs <- list()

# Iterate over each unique value
for (value in unique_values) {
  # Filter the data frame based on the current value
  separated_dfs[[value]] <- subset(feature, Module == value)
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
  DefaultAssay(integrated) <- "RNA"
  integrated <<- AddModuleScore(object = integrated, features = list(cell_marker_gene_list), name = paste0(value, "_Score"))
}

# Apply the function to each unique value using lapply
lapply(unique_values, add_module_score)
head(integrated@meta.data)
#the code above is doing the code below for all the different modules(seperated_dfs)
#cell <- separated_dfs$G2_M
#cell_marker_gene_list <- cell$Gene
#DefaultAssay(df) <- "RNA"
#df <- AddModuleScore(object = df, features = list(cell_marker_gene_list), name = "G2_M")


meta_columns <- c("MES2_Score1", "MES1_Score1", "AC_Score1", "OPC_Score1", "NPC1_Score1", "NPC2_Score1", "G1_S_Score1", "G2_M_Score1")
# Loop through the metadata columns and create DimPlots
for (meta in meta_columns) {
  # Generate the DimPlot
  p <- FeaturePlot(object = integrated, features = meta,  pt.size =0.5,  order = TRUE)  & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  print(p)
  # Save the plot
  ggsave(filename = paste0("Feature_", meta, ".png"), plot = p, width = 10.1, height = 9, dpi=300)
}

saveRDS(integrated, "integrated.rds")

