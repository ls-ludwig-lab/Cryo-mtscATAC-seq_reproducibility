# Load Libraries
library(Seurat)      
library(ggplot2)     
library(tidyverse)   
library(dplyr)       

# Load Seurat Objects
NB <- readRDS("path_to/Harmony_Integrated_Object.rds")   # Neuroblastoma datasets (NB01 & NB02)
GN <- readRDS("path_to/GN01.rds")                        # Ganglioneuroma dataset

# ============================
# nCount ATAC QC
# ============================
## Extract QC metrics from NB
qc_NB <- data.frame(
  Sample = NB$Sample, 
  nCount_peaks = NB$nCount_peaks
)

## Extract QC metrics from GN
qc_GN <- data.frame(
  Sample = "GN01",            # Label GN01 for plotting
  nCount_peaks = GN$nCount_peaks
)

## Combine NB and GN QC data
qc_combined <- rbind(qc_NB, qc_GN)

## Ensure Sample is a factor with the desired order for plotting
qc_combined$Sample <- factor(qc_combined$Sample, levels = c(unique(seurat$Sample), "GN01"))

## Violin plot for nCount ATAC per sample
ggplot(qc_combined, aes(x = Sample, y = nCount_peaks, fill = Sample)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.9, alpha = 0.3) +
  geom_boxplot(width = 0.4, outlier.shape = NA, color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("#76C3A9", "#355C4C", "#b37bb3")) +
  labs(y = "nCount ATAC", x = "") +
  theme_classic() + ylim(0, 10000) +
  theme(text = element_text(colour = "black", size = 15),
        axis.line = element_line(color = "black", size = 0.8),
        axis.text.x = element_text(colour = "black", size = 15, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.ticks.y = element_line(color = "black", size = 0.8),
        axis.ticks.x = element_line(color = "black", size = 0))

# ============================
# TSS Enrichment QC
# ============================
qc_NB <- data.frame(Sample = NB$Sample, TSS_enrichment = NB$TSS.enrichment)
qc_GN <- data.frame(Sample = "GN01", TSS_enrichment = GN$TSS.enrichment)
qc_combined <- rbind(qc_NB, qc_GN)
qc_combined$Sample <- factor(qc_combined$Sample, levels = c(unique(seurat$Sample), "GN01"))

ggplot(qc_combined, aes(x = Sample, y = TSS_enrichment, fill = Sample)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.9, alpha = 0.3) +
  geom_boxplot(width = 0.4, outlier.shape = NA, color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("#76C3A9", "#355C4C", "#b37bb3")) +
  labs(y = "TSS enrichment", x = "") +       # Note: previously "nCount ATAC", should reflect TSS
  theme_classic() + ylim(0, 10) +
  theme(text = element_text(colour = "black", size = 15),
        axis.line = element_line(color = "black", size = 0.8),
        axis.text.x = element_text(colour = "black", size = 15, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.ticks.y = element_line(color = "black", size = 0.8),
        axis.ticks.x = element_line(color = "black", size = 0))

# ============================
# Mitochondrial Depth QC
# ============================
qc_NB <- data.frame(Sample = NB$Sample, mtdepth = NB$mtDNA_depth)
qc_GN <- data.frame(Sample = "GN01", mtdepth = GN$mtDNA_depth)
qc_combined <- rbind(qc_NB, qc_GN)
qc_combined$Sample <- factor(qc_combined$Sample, levels = c(unique(seurat$Sample), "GN01"))

ggplot(qc_combined, aes(x = Sample, y = mtdepth, fill = Sample)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.9, alpha = 0.3) +
  geom_boxplot(width = 0.4, outlier.shape = NA, color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("#76C3A9", "#355C4C", "#b37bb3")) +
  labs(y = "mtDNA depth", x = "") +
  theme_classic() + ylim(0, 120) +
  theme(text = element_text(colour = "black", size = 15),
        axis.line = element_line(color = "black", size = 0.8),
        axis.text.x = element_text(colour = "black", size = 15, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.ticks.y = element_line(color = "black", size = 0.8),
        axis.ticks.x = element_line(color = "black", size = 0))

# ============================
# QC Comparison by Celltype - Neuroblastoma
# ============================
qc_NB <- data.frame(
  Sample = seurat$Sample,
  Celltype = seurat$Manual_Annotation,
  nCount_peaks = seurat$nCount_peaks,
  TSS_enrichment = seurat$TSS.enrichment
)
qc_NB$Celltype <- factor(qc_NB$Celltype, levels = c("Neuroendocrine", "Endothelial", "Fibroblast", "Immune"))

ggplot(qc_NB, aes(x = Celltype, y = TSS_enrichment, fill = Celltype)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.9, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", fill = "white", alpha = 1) +
  scale_fill_manual(values = c("#F1656E", "#A8D7E3", "#616887", "#b37bb3")) +
  labs(y = "TSS") +
  ylim(0, 10) +
  theme_classic() +
  theme(text = element_text(colour = "black", size = 15),
        axis.line = element_line(color = "black", size = 0.8),
        axis.text.x = element_text(colour = "black", size = 0, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.ticks.y = element_line(color = "black", size = 0.8),
        axis.ticks.x = element_line(color = "black", size = 0))

ggplot(qc_NB, aes(x = Celltype, y = nCount_peaks, fill = Celltype)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.9, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", fill = "white", alpha = 1) +
  scale_fill_manual(values = c("#F1656E", "#A8D7E3", "#616887", "#b37bb3")) +
  labs(y = "nCount ATAC") +
  ylim(0, 10000) +
  theme_classic() +
  theme(text = element_text(colour = "black", size = 15),
        axis.line = element_line(color = "black", size = 0.8),
        axis.text.x = element_text(colour = "black", size = 0, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 15),
        axis.ticks.y = element_line(color = "black", size = 0.8),
        axis.ticks.x = element_line(color = "black", size = 0))

# ============================
# QC Comparison by Celltype - Ganglioneuroma
# ============================
qc_gn <- data.frame(
  Sample = "GN01",
  nCount_peaks = GN$nCount_peaks,
  Celltype = GN$Annotation_Manual,
  TSS_enrichment = GN$TSS.enrichment
)
qc_gn$Celltype <- factor(qc_gn$Celltype, levels = c("T_cells", "B_cells", "Endothelial", "Monocyte_Macrophage", "Other", "Fibroblast", "Schwann_cells"))

ggplot(qc_gn, aes(x = Celltype, y = nCount_peaks, fill = Celltype)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.9, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", fill = "white", alpha = 1) +
  scale_fill_manual(values = c("#F9C071", "#e05b02", "#A8D7E3", "#95BBA1", "#658375", "#616887", "#704D6C")) +
  labs(y = "nCount ATAC") +
  ylim(0, 10000) +
  theme_classic() +
  theme(
    text = element_text(colour = "black", size = 15),
    axis.line = element_line(color = "black", size = 0.8),
    axis.text.x = element_text(colour = "black", size = 0, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 15),
    axis.ticks.y = element_line(color = "black", size = 0.8),
    axis.ticks.x = element_line(color = "black", size = 0)
  )

ggplot(qc_gn, aes(x = Celltype, y = TSS_enrichment, fill = Celltype)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.9, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", fill = "white", alpha = 1) +
  scale_fill_manual(values = c("#F9C071", "#e05b02", "#A8D7E3", "#95BBA1", "#658375", "#616887", "#704D6C")) +
  labs(y = "TSS") +
  ylim(0,10) +
  theme_classic() +
  theme(
    text = element_text(colour = "black", size = 15),
    axis.line = element_line(color = "black", size = 0.8),
    axis.text.x = element_text(colour = "black", size = 0, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 15),
    axis.ticks.y = element_line(color = "black", size = 0.8),
    axis.ticks.x = element_line(color = "black", size = 0)
  )

