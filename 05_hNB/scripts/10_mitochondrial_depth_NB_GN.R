# Load required libraries
library(Seurat)      
library(ggplot2)     
library(tidyverse)   
library(dplyr)       
library(zoo)       
library(gghalves)  


# Load Objects and Metadata
NB <- readRDS("path_to/Harmony_Integrated_Object.rds")
GN <- readRDS("path_to/GN01.rds")
metadata_NB <- read.table("Metadata_NB.tsv", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
metadata_GN <- read.table("Metadata_GN.tsv", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

# ============================
# Mitochondrial DNA Depth — NB
# ============================

# UMAP plot of mtDNA depth
ggplot(metadata_NB, aes(umap_1, umap_2, color = mtDNA_depth)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_continuous(type = "viridis", limits = c(5, 100)) +
  theme_classic() + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(
    aspect.ratio = 1,
    legend.background = element_rect(fill = "white"),
    legend.key.size = unit(0.5, 'cm'),
    legend.text = element_text(size = 10),
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line.x = element_line(linetype = "solid", size = 0.5,
                               arrow = grid::arrow(angle = 20, length = unit(0.1, "inches"), type = "closed")),
    axis.line.y = element_line(linetype = "solid", size = 0.5,
                               arrow = grid::arrow(angle = 20, length = unit(0.1, "inches"), type = "closed")),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title   = element_text(hjust = 0)
  ) +
  guides(x = "axis", y = "axis", colour = guide_colorbar(reverse = FALSE))

# UMAP plot faceted by dataset (timepoint)
ggplot(metadata_NB, aes(umap_1, umap_2, color = mtDNA_depth)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_continuous(type = "viridis", limits = c(5, 100)) +
  theme_classic() + xlab("UMAP 1") + ylab("UMAP 2") +
  facet_wrap(~Dataset) +
  theme(aspect.ratio = 1,
        legend.background = element_rect(fill = "white"),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 10),
        legend.position = "right",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_line(linetype = "solid", size = 0.5,
                                   arrow = grid::arrow(angle = 20, length = unit(0.1, "inches"), type = "closed")),
        axis.line.y = element_line(linetype = "solid", size = 0.5,
                                   arrow = grid::arrow(angle = 20, length = unit(0.1, "inches"), type = "closed")),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title   = element_text(hjust = 0)) +
  guides(x = "axis", y = "axis", colour = guide_colorbar(reverse = FALSE))

# ============================
# mtDNA depth per cell type and dataset (NB)
# ============================
metrics_df <- data.frame(
  Celltype = NB$Manual_Annotation,
  Dataset  = NB$Dataset,
  mtdepth  = NB$mtDNA_depth
)

# Define cell type order and dataset-specific transparency
celltype_order  <- c("Neuroendocrine", "Immune", "Endothelial", "Fibroblast")
celltype_colors <- c("#F1656E", "#b37bb3", "#A8D7E3", "#616887")
alpha_values    <- c("Primary" = 0.4, "Second-look" = 0.9)

metrics_df$Celltype <- factor(metrics_df$Celltype, levels = celltype_order)
metrics_df$Dataset  <- factor(metrics_df$Dataset)

# Half-violin and half-boxplot comparison between datasets
ggplot(metrics_df, aes(x = Celltype, y = mtdepth)) +
  # Left half = dataset 1
  geom_half_violin(
    data = subset(metrics_df, Dataset == levels(metrics_df$Dataset)[1]),
    aes(fill = Celltype, alpha = Dataset),
    side = "l", width = 1.1, position = position_nudge(x = 0), color = NA
  ) +
  geom_half_boxplot(
    data = subset(metrics_df, Dataset == levels(metrics_df$Dataset)[1]),
    side = "l", position = position_nudge(x = 0),
    width = 0.2, outlier.shape = NA, color = "black", fill = "white"
  ) +
  # Right half = dataset 2
  geom_half_violin(
    data = subset(metrics_df, Dataset == levels(metrics_df$Dataset)[2]),
    aes(fill = Celltype, alpha = Dataset),
    side = "r", width = 1.1, position = position_nudge(x = 0.05), color = NA
  ) +
  geom_half_boxplot(
    data = subset(metrics_df, Dataset == levels(metrics_df$Dataset)[2]),
    side = "r", position = position_nudge(x = 0.05),
    width = 0.2, outlier.shape = NA, color = "black", fill = "white"
  ) +
  scale_fill_manual(values = celltype_colors) +
  scale_alpha_manual(values = alpha_values) +
  labs(y = "mtDNA depth", x = "", fill = "Celltype", alpha = "Dataset") +
  coord_cartesian(ylim = c(0, 130)) +
  scale_x_discrete(expand = c(0.05, 0)) +
  theme_classic() +
  theme(
    text = element_text(colour = "black", size = 18),
    axis.line = element_line(color = "black", size = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 0),
    axis.text.y = element_text(colour = "black", size = 18),
    axis.ticks.y = element_line(color = "black", size = 0.8),
    axis.ticks.x = element_line(color = "black", size = 0.8)
  )

# ============================
# Mitochondrial DNA Depth — GN
# ============================
qc_metrics_df <- data.frame(
  Cells   = GN$Annotation_Manual,
  mtdepth = GN$mtDNA_depth
)

qc_metrics_df$Cells <- factor(qc_metrics_df$Cells,
                              levels = c("Endothelial", "B_cells", "Other", "Fibroblast", "Schwann_cells",
                                         "Monocyte_Macrophage", "T_cells")
)

ggplot(qc_metrics_df, aes(x = Cells, y = mtdepth, fill = Cells)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.9, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", fill = "white", alpha = 1) +
  labs(y = "mtDNA depth", x = "Cell Types") +
  ylim(0, 120) +
  theme_classic() +
  theme(
    text = element_text(colour = "black", size = 18),
    axis.line = element_line(color = "black", size = 0.8),
    axis.text.x = element_text(colour = "black", size = 0, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 18),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5)
  ) +
  scale_fill_manual(values = c("#A8D7E3", "#e05b02", "#658375", "#616887", "#704D6C", "#95BBA1", "#F9C071"))

# ============================
# Mitochondrial DNA Coverage
# ============================
pull_coverage <- function(SE, cells, resolution = 5) {
  zoo::rollmean(rowMeans(assays(SE)[['coverage']]), resolution)
}

mito.NB01 <- readRDS("cryo_mtscATAC_NB01.rds")
mito.NB02 <- readRDS("cryo_mtscATAC_NB02.rds")
mito.GN01 <- readRDS("cryo_mtscATAC_GN01.rds")

# Smooth coverage along mitochondrial genome
cov_NB01 <- data.frame(pos = zoo::rollmean(1:16569, 5), mito_gotcha = pull_coverage(mito.NB01))
cov_NB02 <- data.frame(pos = zoo::rollmean(1:16569, 5), mito_gotcha = pull_coverage(mito.NB02))
cov_GN01 <- data.frame(pos = zoo::rollmean(1:16569, 5), mito_gotcha = pull_coverage(mito.GN01))

cov <- cbind(cov_NB01, cov_NB02[, -1], cov_GN01[, -1])
colnames(cov)[2:4] <- c("NB01", "NB02", "GN01")

ggplot(cov, aes(x = pos)) +
  geom_line(aes(y = NB01), color = "#76C3A9", linetype = "solid", size = 0.75) +
  geom_line(aes(y = NB02), color = "#355C4C", linetype = "solid", size = 0.75) +
  geom_line(aes(y = GN01), color = "#b37bb3", linetype = "solid", size = 0.75) +
  labs(x = "Position mitochondrial chromosome", y = "mean mtDNA depth") +
  scale_y_continuous(limits = c(0, 65)) +
  theme_classic() +
  theme(
    text = element_text(colour = "black", size = 18),
    axis.line = element_line(color = "black", size = 0.8),
    axis.text.x = element_text(colour = "black", size = 18),
    axis.text.y = element_text(colour = "black", size = 18),
    axis.ticks.y = element_line(color = "black", size = 0.8),
    axis.ticks.x = element_line(color = "black", size = 0.8)
  )


