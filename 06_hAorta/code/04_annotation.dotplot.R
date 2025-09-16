library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(scCustomize)
library(RColorBrewer)
library(readxl)

setwd(".../output")

df <- readRDS(file="preprocessed.rds")

DefaultAssay(df) <- "RNA"


p1 <- DimPlot(df, label.size = 4, label=TRUE) +
  # -- Here you tell ggplot to eliminate all the axis info from your original plot -- #
  theme(aspect.ratio=1/1,
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
p1
ggsave(filename = "UMAP_peaks_res0.5.pdf", plot = p1, width = 5, height = 5)


#manual curated marker of vasculature
# ----------------  Canonical markers ----------------
canonical_markers <- c("PECAM1", "VWF","CLDN5","FLT1","ESAM","ENG","ICAM1","DCN","LUM","COL1A1","COL1A2","FBLN1","FBLN2","CXCL12","MMP2","COL6A2","TIMP1",
                       "ACTA2","TAGLN","MYH11","CNN1","CALD1","LMOD1","MYL9","TPM2","CD68","CD14","C1QA","C1QB","C1QC","LYZ","AIF1","MS4A7","HLA-DRA","APOE",
                       "CD3D","CD3E","CD3G","CD2","CD247","IL7R","CD8A","CD4","NKG7","GZMA","CCL5")

df
dot_plot <- DotPlot_scCustom(
  seurat_object = df,
  features = canonical_markers,
  group.by = "peaks_snn_res.0.5"
) +
  scale_color_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu"))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(x = NULL, y = "Cluster")
dot_plot
ggsave(filename = "canonical_marker.pdf", plot = dot_plot, width = 12, height = 5)


