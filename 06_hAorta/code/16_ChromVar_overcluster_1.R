library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(JASPAR2020)
library(TFBSTools)

df <- readRDS(file="chromVar_seurat.rds")
df
DefaultAssay(df) <- "peaks"
df <- FindClusters(df, resolution = 1.8, algorithm = 3)



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
ggsave(filename = "UMAP_peaks_res1.8.pdf", plot = p1, width = 5, height = 5)

df <- FindClusters(df, resolution = 1.7, algorithm = 3)

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
ggsave(filename = "UMAP_peaks_res1.7.pdf", plot = p1, width = 5, height = 5)


Idents(df) <- "peaks_snn_res.1.8"
DimPlot(df, label = TRUE)
# wilcox is the default option for test.use

da_peaks <- FindMarkers(
  object = df,
  ident.1 = 9,
  ident.2 = 13,
  test.use = 'wilcox',
  min.pct = 0.1
)

head(da_peaks)
open_ident1 <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_ident2 <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])

write.csv(closest_genes_ident1, "open_Cluster9_clonal_fibro.csv")
write.csv(closest_genes_ident2, "open_Cluster13_fibro.csv")

head(closest_genes_ident1)

DefaultAssay(df2) <- 'peaks'




pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

DefaultAssay(df) <- 'chromvar'

differential.activity <- FindMarkers(
  object = df,
  ident.1 = 9,
  ident.2 = 13,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

differential.activity

# look at the activity of Mef2c
p2 <- FeaturePlot(
  object = df,
  features = "MA0083.3",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  order = TRUE
)
p2
ggsave(plot = p2, filename = "MA0083.3_ChromVar.pdf", width = 5, height = 5, dpi=300)

write.csv(differential.activity, "differential.activity.Fibro9vs13.csv")


