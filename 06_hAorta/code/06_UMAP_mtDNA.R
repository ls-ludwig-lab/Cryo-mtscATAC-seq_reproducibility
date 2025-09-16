library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(readxl)

setwd(".../output/")


df <- readRDS(file="preprocessed.rds")

meta.data <- df@meta.data

umap_coords <- Embeddings(df, "umap")
umap_df <- as.data.frame(umap_coords)
meta.data <- cbind(meta.data, umap_df)



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

clustering.color <- c("SMC" = "deepskyblue3",
                      "SMC I" = "skyblue2",
                      "Fibroblasts" = "darkseagreen3",
                      "Macrophages" = "plum2",
                      "NK/T Cells" = "lightsalmon",
                      "Endothelial Cells" = "brown2")

p <- ggplot(meta.data, aes(x = UMAP_1, y =UMAP_2, color = Celltype)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = clustering.color) +
  theme(aspect.ratio=1/1,
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        axis.line = element_blank())
p
ggsave(plot = p, filename = "UMAP_Cluster_res0p5.pdf", width = 5, height = 5, dpi=300)


p1 <- ggplot(meta.data, aes(x = Celltype, y = mtDNA_depth, fill = Celltype)) +
  scale_fill_manual(values = clustering.color) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", color = "grey30") +
  theme_classic() +
  ylab("mtDNA depth") +
  ylim(0, 150) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 20),
    axis.line.x = element_line(),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_blank()
  )
p1
ggsave(plot = p1, filename = "VlnPlot_mtDNA.pdf", width = 4, height = 6, dpi = 300)
