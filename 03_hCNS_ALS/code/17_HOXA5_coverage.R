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
library(tidyr)
library(ggrepel)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)


#### Microglia subset analysis ####
#############################

setwd(".../output")

df <- readRDS(file=".../Microglia_subset.rds")

R_5 <- read.csv("df_non_wt_R5.csv")
R_10 <- read.csv("df_non_wt_R10.csv")
R_28 <- read.csv("df_non_wt_R28.csv")
R_29 <- read.csv("df_non_wt.csv")

combined <- rbind(R_5, R_10, R_28, R_29)

head(df@meta.data)

barcodes <- combined$X

# Create a logical vector indicating if each cell barcode is in the `barcode_vector`
df$X1345G_A_clone <- ifelse(rownames(df@meta.data) %in% barcodes, "clone", "wt")


p <- DimPlot(df, group.by = "X1345G_A_clone")
p

meta.data <- df@meta.data
umap_coords <- Embeddings(df, "umap")
umap_df <- as.data.frame(umap_coords)
meta.data <- cbind(meta.data, umap_df)

color.palette <-c("firebrick","grey" )
p <- ggplot(meta.data, aes(x = UMAP_1, y =UMAP_2, color = X1345G_A_clone)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c(color.palette)) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        axis.line = element_blank())
p


DefaultAssay(df) <- 'ATAC'

Idents(df) <- "X1345G_A_clone"
# wilcox is the default option for test.use

# change back to working with peaks instead of gene activities
DefaultAssay(df) <- 'ATAC'
Idents(df)
# wilcox is the default option for test.use
da_peaks <- FindMarkers(
  object = df,
  ident.1 = "wt",
  ident.2 = "clone",
  test.use = 'wilcox',
  min.pct = 0.1
)

head(da_peaks)

open_wt <- rownames(da_peaks[da_peaks$avg_log2FC > 1.5, ])
open_clone <- rownames(da_peaks[da_peaks$avg_log2FC < -1.5, ])
closest_open_clone <- ClosestFeature(df, regions = open_clone)
closest_open_wt <- ClosestFeature(df, regions = open_wt)

write.csv(closest_open_wt, "closest_open_wt.csv")
write.csv(closest_open_clone, "closest_open_clone.csv")

p1 <- CoveragePlot(
  object = df,
  region = "HOXA5",
  annotation = TRUE,
  peaks = FALSE,
  extend.upstream = 10000,
  extend.downstream = 10000,
)
p1 
ggsave(filename = "HOXA5_coverage.pdf", plot = p1,  width = 10, height = 5)