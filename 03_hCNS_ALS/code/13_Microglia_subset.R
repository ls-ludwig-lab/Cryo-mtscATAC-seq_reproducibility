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

unique(d@meta.data$possible_Celltype)

df <- subset(d, subset = possible_Celltype == "Microglia" | possible_Celltype == "proinflammatory Microglia")

p <- DimPlot(df, group.by = "possible_Celltype")
p

# Count the number of cells per category combination
cell_counts <- df@meta.data %>%
  group_by(dataset, possible_Celltype) %>%
  summarise(count = n(), .groups = "drop")
cell_counts


#  a new UMAP using the integrated embeddings
df <- RunUMAP(df, reduction = "harmony", dims = 2:30)

DefaultAssay(df) <- "ATAC"
df <- FindNeighbors(object = df, reduction = 'harmony', dims = 2:30)
df <- FindClusters(object = df, resolution = 0.4, verbose = FALSE, algorithm = 3)

p1 <- DimPlot(df, group.by = "possible_Celltype", label.size = 4, label=TRUE) +
  # -- Here you tell ggplot to eliminate all the axis info from your original plot -- #
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) 
p1


p2 <- DimPlot(df, group.by = "seurat_clusters", label.size = 4, label=TRUE) +
  # -- Here you tell ggplot to eliminate all the axis info from your original plot -- #
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) 
p2


p3 <- DimPlot(df, group.by = "dataset", label.size = 4, label=TRUE) +
  # -- Here you tell ggplot to eliminate all the axis info from your original plot -- #
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) 
p3


saveRDS(df, "Microglia_subset.rds")



