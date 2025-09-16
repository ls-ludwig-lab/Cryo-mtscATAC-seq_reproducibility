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
library(BuenColors)
library(rstatix)
library(ComplexHeatmap)
library(ggrepel)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BuenColors)
library(Matrix)
library(GenomicRanges)


setwd(".../output/")
df <- read.csv("TopSMC_clones_meta.csv")
seurat <- readRDS(file="preprocessed.rds")

seurat
DimPlot(seurat)
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

DefaultAssay(seurat) <- "peaks"
# add motif information


seurat <- AddMotifs(
  object = seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm)


main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- as.logical(seqnames(granges(seurat)) %in% main.chroms)
seurat <- seurat[keep.peaks, ]

seurat <- RunChromVAR(
  object = seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(seurat) <- 'chromvar'

seurat@meta.data <- seurat@meta.data %>%
  mutate(Celltype = case_when(
    peaks_snn_res.0.5 == 0 ~ "SMC",
    peaks_snn_res.0.5 == 1 ~ "SMC I",
    peaks_snn_res.0.5 == 2 ~ "Fibroblasts",
    peaks_snn_res.0.5 == 3 ~ "Macrophages",
    peaks_snn_res.0.5 == 4 ~ "NK/T Cells",
    peaks_snn_res.0.5 == 5 ~ "Endothelial Cells",
    TRUE ~ "Unknown"
  ))
Idents(seurat) <- "Celltype"

unique(seurat$Celltype)



saveRDS(seurat, "chromVar_seurat.rds")

