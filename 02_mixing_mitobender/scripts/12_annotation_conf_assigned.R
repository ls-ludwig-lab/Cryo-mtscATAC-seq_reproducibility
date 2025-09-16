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
# Define color schemes for plotting cell types and dataset regions
output_dir <- "..."

df <- readRDS(file.path(output_dir,file="subset_mitoFiltered_both.rds"))
DefaultAssay(df) <- "peaks"

DimPlot(df)


gene.activities <- GeneActivity(df)
# add the gene activity matrix to the Seurat object as a new assay and normalize it

df[['RNA']] <- CreateAssayObject(counts = gene.activities)

df <- NormalizeData(
  object = df,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(df$nCount_RNA)
)
Idents(df)

DefaultAssay(df) <- "RNA"
marker<- FindAllMarkers(df,  min.pct = 0.25)

# Identify unique values
unique_values <- unique(marker$cluster)

# Create an empty list to store the individual data frames
separated_dfs <- list()

# Iterate over each unique value
for (value in unique_values) {
  # Filter the data frame based on the current value
  separated_dfs[[value]] <- subset(marker, cluster == value)
}

# Access individual data frames
for (value in unique_values) {
  cat(paste("Data frame for '", value, "':\n", sep = ""))
  print(separated_dfs[[value]])
  cat("\n")
}

# Define the path where CSV files will be saved
output_directory <- "..."

# Ensure the output directory exists
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# Save individual data frames as CSV files
for (value in unique_values) {
  cat(paste("Data frame for '", value, "':\n", sep = ""))
  print(separated_dfs[[value]])
  cat("\n")
  
  # Define the file name for the CSV
  file_name <- file.path(output_directory, paste("data_frame_mitoFiltered_", value, ".csv", sep = ""))
  
  # Save the data frame as a CSV
  write.csv(separated_dfs[[value]], file_name, row.names = FALSE)
  cat(paste("Data frame for '", value, "' saved as '", file_name, "'.\n", sep = ""))
}

DimPlot(df)



marker_select <- list(
  "excitatory neurons" = c('ENC1', 'SLC17A7', 'TAFA2', 'NTNG1', 'PCP4', 'FSTL4'),
  "Astrocytes" = c('SLC1A2', 'GFAP', 'AQP4'),
  "OPCs" = c('PTPRZ1', 'VCAN', 'PCDH15'),
  "Oligodendrocytes" = c('TMTC2', 'ANLN', 'LRP2', 'PLP1', 'ST18', 'BCAS1'),
  "Microglia" = c('SPP1', 'APBB1IP', 'CD74'),
  "Endothelial Cells" = c('FLT1', 'MECOM', 'COBLL1'),
  "inhibitory neurons" = c('INPP4B', 'EML5', 'GRM1', 'KCNC2', 'KIAA1217', 'NXPH1'),
  "GABA inhibitory neurons" = c('RGS12', 'CXCL14', 'GAD1', 'FGF13', 'TMEM132D', 'KCNC2'),
  "Glia Cells" = c('ADAMTSL1', 'LINC01727', 'ALDH1A1'),
  "Fibroblasts" = c('CEMIP', 'DCN', 'CEMIP'),
  "Ependymal" = c('CFAP54', 'DNAH9', 'CFAP299')
)

# Flatten marker list
marker_df <- stack(marker_select)
marker_df
# Remove duplicated genes (keep first occurrence)
marker_df <- marker_df[!duplicated(marker_df$values), ]
marker_df
# Preserve original order of cell types
marker_df$ind <- factor(marker_df$ind, levels = names(marker_select))

# Group genes per cell type
marker_list_clean <- marker_df$values
marker_list_clean# Run DotPlot with scCustomize
p1 <- DotPlot( df,
  features = marker_list_clean
)
p1
ggsave(plot = p1, "DotPlot_markerGenes.pdf", width = 20, height = 10)


DimPlot(df, label = TRUE)

FeaturePlot(df, c('PTPRZ1', 'VCAN', 'PCDH15'))

genes <- c(
  "SOX4", "CDK4", "TUBB2B", "SOX11", "MDK", "MEST", "DLL3", "SOX2", "DCTN2", "CCND2", 
  "C1orf61", "VIM", "FABP7", "NPPA", "MT3", "ETV1", "MEG3", "EGFR", "PTPRZ1", "PTN", 
  "CLU", "MT1X", "DTX3", "SPARCL1", "CST3", "BCAN", "GPR17", "SCRG1", "TNR", "VCAN", 
  "FXYD6", "GPM6A", "FYN"
)

FeaturePlot(df, c('CD68'), order = TRUE)


df <- RenameIdents(df,"0" = "Oligodendrocytes",
                   "1" = "Glia Cells",
                   "2" = "Oligodendrocytes I",
                   "3" = "Projection Nerons (cortical glutameric lineage)",
                   "4" = "Microglia",
                   "5" = "Astrocytes",
                   "6" = "Vascular Cells",
                   "7" = "Ependymal Cells",
                   "8" = "Excitatory Neurons (pyramidal, glutamatergic)",
                   "9" = "Astrocyte-like neural stem cells",
                   "10" = "OPCs",
                   "11" = "Astrocytes I",
                   "12" = "Lymphocytes",
                   "13" = "non neuronal mesoderm")


Idents(df)
# Save cluster identities as a new column in meta.data
df$possible_Celltype <- df@active.ident

meta.data <- df@meta.data

umap_coords <- Embeddings(df, "umap")
umap_df <- as.data.frame(umap_coords)
meta.data <- cbind(meta.data, umap_df)

DimPlot(df, group.by = "assign", label = TRUE)

write.csv(meta.data, "meta.data.csv")
saveRDS(df, file = file.path(output_dir, "mitoFilter_df_annotated.rds"))

