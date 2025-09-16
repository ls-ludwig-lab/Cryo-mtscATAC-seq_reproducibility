library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(BuenColors)
"%ni%" <- Negate("%in%")

df <- readRDS(file=".../rds_files/unfiltered_0p5pct.rds")
d0 <- readRDS(file=".../output/variants_donor_specific.rds")
d0


head(colnames(df))

cell_variant.matrix <- assays(d0)$coverage*assays(d0)$allele_frequency

#depth per barcode
d0@colData@listData[["depth"]]

#mean Depth/Barcode
mean(d0@colData@listData[["depth"]])

df_prim <- data.frame(rowData(d0))


#set parameters: minimum number of reads and cells for supporting mutation event
coverage_min <- 5
depth_min <- 5
mean_min <- 0.001
nreads_min <- 5

#Operator
`%ni%` <- Negate(`%in%`)

df
metadata <- df@meta.data

# Subset SummarizedExperiment based on data frame column names
d0 <- d0[, colnames(d0) %in% colnames(df)]
d0

#annotated heatmap function
#make cell x variant matrix with heteroplasmy values only with mutations supported by at least n number of mutant reads
#I already have my variants of interest in this SE
cell_variant.matrix <- assays(d0)$coverage*assays(d0)$allele_frequency >= nreads_min
cell_variant.matrix <- cell_variant.matrix*assays(d0)$allele_frequency #sqrt(cell_variant.matrix*assays(se)$allele_frequency)
Cell_Variant.dim <- dim(cell_variant.matrix)
Cell_Variant.dim


#subset to mito and atac quality filtered cells
cells <- as.data.frame(colData(d0)) %>% 
  filter(depth > depth_min) %>% rownames(.) %>% 
  intersect(., rownames(metadata))

cell_variant.matrix <- cell_variant.matrix[,cells]
dim(cell_variant.matrix)

#identify cells with at least one informative variant
pos.cells <- colSums(cell_variant.matrix)>0
cells <- names(pos.cells[pos.cells==TRUE])

#subset cell variant matrix to cells of interest
cell_variant.matrix <- cell_variant.matrix[,cells]
cell_variant.dim <- dim(cell_variant.matrix)
cell_variant.dim

#compute dendrogram
dist.mtx = dist(x = t(cell_variant.matrix>0), method = "binary")
hca = hclust(d = dist.mtx, method = "ward.D2")
dend = as.dendrogram(hca)

##make variant annotation
#mean
mean.annotation <- setNames(rowData(d0)$mean, rownames(d0))

#ncells
ncells.annotation <- setNames(rowData(d0)$n_cells_over_5, rownames(d0))
ncells.annotation <- ncells.annotation[cells]

depth.annotation <-  setNames(colData(d0)$depth, colnames(d0))
depth.annotation <- depth.annotation[cells]


donor.annotation <- setNames(metadata$vireo_donor, rownames(metadata))[cells]

nCount_peaks.annotation <- setNames(metadata$nCount_peaks, rownames(metadata))[cells]
doublet.annotation <- setNames(metadata$vireo_doublet_logLikRatio, rownames(metadata))[cells]
nFeature_mito.annotation <- setNames(metadata$nFeature_mito, rownames(metadata))[cells]

# Annotations
donor.annotation <- setNames(metadata$vireo_donor, rownames(metadata))[cells]
nCount_peaks.annotation <- setNames(metadata$nCount_peaks, rownames(metadata))[cells]
doublet.annotation <- setNames(metadata$vireo_doublet_logLikRatio, rownames(metadata))[cells]
nFeature_mito.annotation <- setNames(metadata$nFeature_mito, rownames(metadata))[cells]
amulet.annotation  <- setNames(metadata$doublet, rownames(metadata))[cells]


continuous_palette1 <- jdb_palette("flame_macaw", type = "discrete")
continuous_palette2 <- jdb_palette("blue_cyan", type = "discrete")
continuous_palette3 <- jdb_palette("brewer_violet", type = "discrete")

col_fun_seqdepth <- colorRamp2(c(min(nCount_peaks.annotation), max(nCount_peaks.annotation)), c(continuous_palette1[1], continuous_palette2[length(continuous_palette2)]))
col_fun_doubletq <- colorRamp2(c(min(doublet.annotation), max(doublet.annotation)), c(continuous_palette2[1], continuous_palette1[length(continuous_palette2)]))
col_fun_nFeature <- colorRamp2(c(min(nFeature_mito.annotation), max(nFeature_mito.annotation)), c(continuous_palette3[1], continuous_palette3[length(continuous_palette3)]))

library(viridisLite)

col_fun_seqdepth <- colorRamp2(
  quantile(nCount_peaks.annotation, probs = c(0.05, 0.95), na.rm = TRUE),
  viridis(2)
)

ha <- HeatmapAnnotation(
  depth = anno_barplot(depth.annotation),
  donor = donor.annotation,
  seq.depth = nCount_peaks.annotation, # ✅ Fixed here
  amulet = amulet.annotation,
  
  annotation_name_side = "left",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 8),
  simple_anno_size = unit(6, "mm"),
  
  col = list(
    donor = c("donor0" = "#5FCDD9", "donor1" = "#F00D09", "doublet" = "orange2", "unassigned" = "#ccd5ae"),
    assign = c("donor0" = "#5FCDD9", "donor1" = "#F00D09", "Collision" = "orange1", "Low_coverage" = "#D7E0B8"),
    amulet = c("singlet" = "#9284F0", "multiplet" = "#F0E285"),
    `seq.depth` = col_fun_seqdepth  # ✅ Keep this too
  ),
  
  border = TRUE
)


ht <- Heatmap(
  as.matrix(cell_variant.matrix),
  col = as.character(jdb_palette("solar_rojos", type = "continuous")),
  show_column_names = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  cluster_columns = dend,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  top_annotation = ha,
  border = TRUE
)
ht

pdf(file="ComplexHeatmap_0p5_unfiltered.pdf", width = 10, height = 5)
ht
dev.off()



