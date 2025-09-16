library(data.table)
library(dplyr)
library(stringr)
library(Signac)
library(SummarizedExperiment)
library(ggplot2)
library(ComplexHeatmap)
library(tibble)
library(tidyr)
library(circlize)


setwd("/data/cephfs-1/work/projects/ludwig-mtscatac-als/ALS/2024_Region")

#load the variants called using the 00_variants_calling.R output rds
#data will be available upon request

area1 <- readRDS(file="output/Region_5/output/primary.rds")
area2 <- readRDS(file="output/Region_10/output/primary.rds")
area3 <- readRDS(file="output/Region_28/output/primary.rds")
area4 <- readRDS(file="output/Region_29/output/primary.rds")


vars_area1 <- data.frame(rowData(area1)) %>%
  filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 5) %>%
  pull(variant)

vars_area2 <- data.frame(rowData(area2)) %>%
  filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 5) %>%
  pull(variant)

vars_area3 <- data.frame(rowData(area3)) %>%
  filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 5) %>%
  pull(variant)

vars_area4 <- data.frame(rowData(area4)) %>%
  filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 5) %>%
  pull(variant)


vars_combined <- unique(c(vars_area1, vars_area2, vars_area3, vars_area4))

area1_filtered <- area1[vars_combined, ]
area2_filtered <- area2[vars_combined, ]
area3_filtered <- area3[vars_combined, ]
area4_filtered <- area4[vars_combined, ]

allele_freq_df <- data.frame(
  variant = rownames(area1_filtered),
  area1 = rowData(area1_filtered)$mean,
  area2 = rowData(area2_filtered)$mean,
  area3 = rowData(area3_filtered)$mean,
  area4 = rowData(area4_filtered)$mean
)


allele_freq_df <- allele_freq_df %>%
  mutate(
    log10_FC_area2_vs_area1 = log10(area2 / area1),
    log10_FC_area3_vs_area1 = log10(area3 / area1),
    log10_FC_area4_vs_area1 = log10(area4 / area1)
  )

# Assuming your dataframe is called `allele_freq_df`
allele_freq_df <- allele_freq_df[allele_freq_df$variant != "310T>C", ]
allele_freq_df <- allele_freq_df[allele_freq_df$variant != "189A>G", ]


library(ComplexHeatmap)

# Data preparation (similar to above)
heatmap_data <- allele_freq_df %>%
  dplyr::select(variant, area1, area2, area3, area4) %>%
  pivot_longer(cols = starts_with("area"), names_to = "region", values_to = "heteroplasmy") %>%
  pivot_wider(names_from = region, values_from = heteroplasmy) %>%
  column_to_rownames("variant")



pdf(file="ComplexHeatmap_mito.pdf", width = 10, height = 60)
# Create a heatmap
p1 <- Heatmap(
  as.matrix(heatmap_data),
  name = "Heteroplasmy",
  col = colorRamp2(
    c(0, 0.001, 0.005, 0.01, 0.015),
    brewer.pal(n = 5, name = "YlGnBu")
  ),
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  heatmap_legend_param = list(title = "Heteroplasmy %"),
  column_title = "Brain Regions",
  row_title = "Variants",
  border = TRUE)
p1
dev.off()


