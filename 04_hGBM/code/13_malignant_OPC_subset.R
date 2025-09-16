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

setwd("/data/cephfs-1/work/projects/ludwig-ff-mtscatac/mtscATAC_brain/10_2024_resequence_GBM/integration/output")

df <- readRDS(file="integrated_annotated.rds")

head(df@meta.data)

unique(df@meta.data$possible_Celltype_broader)


df <- subset(df,possible_Celltype_broader == "OPC-like" |
               possible_Celltype_broader == "Malignant OPC-like"
             )

p <- DimPlot(df, group.by = "possible_Celltype_broader")
p


# create a new UMAP using the integrated embeddings
df <- RunUMAP(df, reduction = "integrated_lsi", dims = 2:30)
DefaultAssay(df) <- "ATAC"
df <- FindNeighbors(object = df, reduction = 'integrated_lsi', dims = 2:30)
df <- FindClusters(object = df, resolution = 0.2, verbose = FALSE, algorithm = 3)

DimPlot(df)
DimPlot(df, group.by = "possible_Celltype_broader")
DimPlot(df, group.by = "dataset")

head(df@meta.data)


Idents(df)

library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

DefaultAssay(df) <- "ATAC"
# add motif information
df <- AddMotifs(df, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)

Idents(df) <- "possible_Celltype_broader"

da_peaks <- FindMarkers(
  object = df,
  ident.1 = "OPC-like",
  ident.2 = "Malignant OPC-like",
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_ATAC'
)

#filter for chr
filtered_da_peaks <- da_peaks[grepl("^chr", row.names(da_peaks)), ]

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.2 > 0.2, ])


top.da.peak <- rownames(filtered_da_peaks[filtered_da_peaks$p_val < 0.005, ])


# test enrichment
enriched.motifs <- FindMotifs(
  object = df,
  features = top.da.peak
)

p1 <- MotifPlot(
  object = df,
  motifs = head(rownames(enriched.motifs))
)
p1



head(da_peaks)

open_cd_OPC <- rownames(da_peaks[da_peaks$avg_log2FC > 1.5, ])
open_maligOPC <- rownames(da_peaks[da_peaks$avg_log2FC < -1.5, ])

closest_genes_OPC <- ClosestFeature(df, regions = open_cd_OPC)
closest_genes_malig_OPC <- ClosestFeature(df, regions = open_maligOPC)


ggplot(data=da_peaks, aes(x=avg_log2FC, y=-log10(p_val))) + geom_point() 

# gather the footprinting information for sets of motifs
df <- Footprint(
  object = df,
  motif.name = c("SOX2", "OLIG2", "BACH1", "JUN", "SOX10"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = TRUE
)


# plot the footprint data for each group of cells
p2 <- PlotFootprint(df, features = c("SOX2", "OLIG2", "BACH1", "JUN", "SOX10"))
p2
ggsave(plot = p2, filename = "TF.pdf", width = 50, height = 20, dpi = 300, limitsize = FALSE)


"MYT1"


regions_to_plot <- c("MYT1","SOX2", "OLIG2", "BACH1", "JUN", "SOX10")

# Loop through each gene in the list
for (gene in regions_to_plot) {
  # Generate the coverage plot for the current gene
  p <- CoveragePlot(
    object = df,
    region = gene,
    extend.upstream = 1000,
    extend.downstream = 1000
  )
  
  # Save the plot with a unique filename
  ggsave(plot = p, filename = paste0("Coverage_", gene, ".pdf"),
         width = 8, height = 6, dpi = 300, limitsize = FALSE)
}
p
