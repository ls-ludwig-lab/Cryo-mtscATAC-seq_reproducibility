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

output_dir <- "..."

df <- readRDS(".../filtered_0p5pct.rds")
df

#filtering, optional
df <- subset(df, doublet == "singlet")
df <- subset(df, mtDNA_depth > 5)
df


d_mito <- read_tsv(".../f0p5pct_mitoBender_data_meta.tsv", show_col_types = FALSE)
head(d_mito)
meta_df <- dplyr::select(d_mito, cell_id, assign) %>%
  tibble::column_to_rownames("cell_id")

df <- AddMetaData(df, metadata = meta_df)
df <- subset(df, assign == "donor1" | assign == "donor0")
df


DefaultAssay(df) <- "peaks"
df <- RunTFIDF(df)
df <- FindTopFeatures(df, min.cutoff = 'q0')
df <- RunSVD(df)
df <- RunUMAP(object = df, reduction = 'lsi', dims = 2:30)
df <- FindNeighbors(object = df, reduction = 'lsi', dims = 2:30)
df <- FindClusters(object = df, resolution = 0.6, verbose = FALSE, algorithm = 3)
DimPlot(object = df, label = TRUE) + NoLegend()


p1 <- DimPlot(df, label.size = 4, label=TRUE) +
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


saveRDS(df, "subset_mtitoFiltered_both.rds")



