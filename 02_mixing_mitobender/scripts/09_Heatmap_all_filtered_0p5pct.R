library(Seurat)
library(readr)
library(dplyr)

# Define paths to input
condition_list <- list(
  donor_specific = list(
    variant_path = ".../variants_donor_specific.rds",
    meta_path = ".../rds_files/f0p5pct_data_meta.tsv"
  ),
  mitoBender = list(
    variant_path = "../mitoBender/variants_donor_specific.rds",
    meta_path = "../rds_files/f0p5pct_mitoBender_data_meta.tsv"
  )
)

# Load Seurat object once
df <- readRDS(".../rds_files/filtered_0p5pct.rds")
dim(df)

add_assign_to_seurat <- function(seurat_obj, tsv_path, colname = "assign") {
  stopifnot(is.character(colname), length(colname) == 1)
  
  tab <- readr::read_tsv(tsv_path, show_col_types = FALSE)
  
  # Try common variants of column names
  cell_cols   <- intersect(c("cell_id","cell","barcode","cellID","CellID"), colnames(tab))
  assign_cols <- intersect(c("assign","assignment","label","call","predicted"), colnames(tab))
  
  if (length(cell_cols) == 0 || length(assign_cols) == 0) {
    stop(sprintf(
      "Couldn't find required columns. Available columns: %s\nLooked for cell in %s; assign in %s",
      paste(colnames(tab), collapse = ", "),
      paste(c("cell_id","cell","barcode","cellID","CellID"), collapse = ", "),
      paste(c("assign","assignment","label","call","predicted"), collapse = ", ")
    ))
  }
  
  tab <- dplyr::select(tab, cell_id = dplyr::all_of(cell_cols[1]),
                       assign  = dplyr::all_of(assign_cols[1])) |>
    dplyr::distinct(cell_id, .keep_all = TRUE)
  
  assign_vec <- setNames(tab$assign, tab$cell_id)[colnames(seurat_obj)]
  seurat_obj[[colname]] <- assign_vec
  seurat_obj
}

df <- add_assign_to_seurat(
  seurat_obj = df,
  tsv_path   = condition_list$donor_specific$meta_path,
  colname    = "assign_donor_specific"
)

# (optional) add mitoBender too
df <- add_assign_to_seurat(
  seurat_obj = df,
  tsv_path   = condition_list$mitoBender$meta_path,
  colname    = "assign_mitoBender"
)

# sanity check
table(df$assign_mitoBender, useNA = "ifany")


df

df <- subset(
  df,
  subset = doublet == "singlet" &
    mtDNA_depth > 5 &
    assign_donor_specific %in% c("donor0","donor1")
  # If you prefer the mitoBender column instead, swap to assign_mitoBender
)
df
# Refresh metadata object you pass to your heatmap code
metadata <- df@meta.data
head(metadata)


d0 <- readRDS(file="...t/data/mitoBender/variants_donor_specific.rds")
d0


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
assign.annotation <- setNames(metadata$assign_donor_specific, rownames(metadata))[cells]
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

#Donor0 and Donor1 are assigned randomly by vireo. 0,1% and 0,5% are not assigned the same. Here Donor0 (0,1%) == Donor 1 (0,5%)

ha <- HeatmapAnnotation(
  depth = anno_barplot(depth.annotation),
  donor = donor.annotation,
  assign = assign.annotation,
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

pdf(file="ComplexHeatmap_0p5_mitobender_all_filtered.pdf", width = 10, height = 5)
ht
dev.off()




