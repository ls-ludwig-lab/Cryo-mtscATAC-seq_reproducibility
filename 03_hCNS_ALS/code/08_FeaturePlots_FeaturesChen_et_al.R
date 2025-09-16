library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(readxl)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(scCustomize)

setwd("/data/cephfs-1/work/projects/ludwig-mtscatac-als/ALS/2024_BB28/test_output")

df <- readRDS(file="integrated_annotated.rds")
p <- DimPlot(df, group.by = "dataset")
p

DefaultAssay(df) <- "RNA"

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

desired_order <- c(
  "Microglia",
  "proinflammatory Microglia",
  "Oligodendrocytes",
  "Oligodendrocytes I",
  "OPCs",
  "OPCs I",
  "OPCs II",
  "Astrocytes",
  "Astrocytes I",
  "Glia Cells",
  "Glia Cells I",
  "Glia Cells II",
  "Glia Cells III",
  "excitatory neurons",
  "GABA inhibitory neurons",
  "inhibitory neurons",
  "Endothelial Cells",
  "CD8 T-Cells"
)

# Make sure Idents are set
Idents(df) <- "possible_Celltype"

# Reorder factor levels
df$possible_Celltype <- factor(df$possible_Celltype, levels = desired_order)



p1 <- DotPlot_scCustom(
  seurat_object = df,
  features = marker_list_clean,
  cluster.idents = TRUE,
  x_lab_rotate = TRUE
)
p1
ggsave(plot = p1, filename = "DotPlot_Chen_etal.pdf", width = 20, height = 10, dpi = 300, limitsize = FALSE)

