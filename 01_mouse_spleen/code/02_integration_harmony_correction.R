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
library(harmony)
library(readxl)

# Define base data directory and set working directory (for portability)
setwd("../../integration/output")
# Load combined Seurat object
df <- readRDS(file="combined_1.rds")

unique(df@meta.data$dataset)

# Subset for datasets/fixation conditions
df3 <- subset(df, subset = dataset == "FF_10min")
df4 <- subset(df, subset = dataset == "FF30_min")
df5 <- subset(df, subset = dataset == "Cells")

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(df3, df4, df5),
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

ggsave(plot = p1, filename = "DimPlot_Integrated_dataset.png", width = 7, height = 7, dpi=300)


saveRDS(integrated, "integrated.rds")

#Harmony Integration

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
hm.integrated <- FindClusters(object = hm.integrated, resolution = 1, verbose = FALSE, algorithm = 3)

hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony')

p5 <- DimPlot(hm.integrated, group.by = 'dataset', pt.size = 0.1) + ggplot2::ggtitle("Harmony integration")
p5

p6 <- DimPlot(hm.integrated, pt.size = 0.1) 
p6

p8 <- p5 + p6
p8

gene.activities <- GeneActivity(hm.integrated)

# add to the Seurat object as a new assay
hm.integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)

hm.integrated <- NormalizeData(
  object = hm.integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(hm.integrated$nCount_RNA)
)
DefaultAssay(hm.integrated) <- "RNA"
marker<- FindAllMarkers(hm.integrated,  min.pct = 0.25)


#Features are manually curated from Ambre Gigulay in the lab
feature <- as.data.frame(read_excel("/../../Spleen_CellMarkers_mouse.xlsx"))
feature

unique_values <- unique(feature$alias)
unique_values

# Create an empty list to store the individual data frames
separated_dfs <- list()

# Iterate over each unique value
for (value in unique_values) {
  # Filter the data frame based on the current value
  separated_dfs[[value]] <- subset(feature, alias == value)
}

# Access individual data frames
for (value in unique_values) {
  cat(paste("Data frame for '", value, "':\n", sep = ""))
  print(separated_dfs[[value]])
  cat("\n")
}

#[1] "Marginal_zone_B_cell"                    
#[2] "T_cell"                                  
#[3] "NK_cell"                                 
#[4] "Dendritic cell_S100a4_high"              
#[5] "Monocyte"                                
#[6] "Neutrophil_S100a8_high"                  
#[7] "Macrophage"                              
#[8] "B_cell"                                  
#[9] "Plasmacytoid_dendritic_cell_Siglech_high"
#[10] "Erythroid_cell"                          
#[11] "Neutrophil_Ngp_high" 


add_module_score <- function(value) {
  cell_marker_gene_list <- separated_dfs[[value]]$gene
  print(paste("Processing:", value))  # Check which alias is being processed
  print(cell_marker_gene_list)  # Check if genes are present
  DefaultAssay(hm.integrated) <- "RNA"
  hm.integrated <<- AddModuleScore(object = hm.integrated, features = list(cell_marker_gene_list), name = paste0(value, "_Score"))
}

#Verify That separated_dfs Contains Gene Lists
print(separated_dfs[[unique_values[1]]])
print(colnames(feature))
print(head(feature))
cell_marker_gene_list <- separated_dfs[[value]]$gene  # Adjust column name if needed
cell_marker_gene_list

for (value in unique_values) {
  if (!is.null(separated_dfs[[value]]$gene) && length(separated_dfs[[value]]$gene) > 0) {
    add_module_score(value)
  } else {
    print(paste("Skipping", value, "due to empty gene list."))
  }
}


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
  p <- FeaturePlot(object = hm.integrated, features = meta,  pt.size =0.5,  order = TRUE)  & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  # Display the plot in the plot panel
  print(p)
  
  # Save the plot
  ggsave(filename = paste0("Feature_", meta, ".png"), plot = p, width = 6, height = 6, dpi=300)
  
}

T_cell <- c("Cd4", "Foxp3", "Cd8a","Ifng", "Ncr1", "Fcgr3")

DefaultAssay(hm.integrated) <- "RNA"
p1 <- FeaturePlot(hm.integrated, 
            features =  T_cell,
            order = TRUE)
p1
ggsave(p1, filename = "T_cell_marker.png", width = 6, height = 6, dpi=300)

saveRDS(hm.integrated, "hm.integrated.rds")
