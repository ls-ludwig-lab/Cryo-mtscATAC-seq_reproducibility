library(dplyr)
library(ggplot2)
library(viridis)
library(pheatmap) 
library(SummarizedExperiment)
library(ComplexHeatmap)
library(tidyr)
library(circlize)
library(viridis)    # for viridis palette


setwd("/.../Figures")

process_donor <- function(output_dir, primary_path, subset_path, donor_id) {
  
  # Load variants and subset cells
  primary <- readRDS(primary_path)
  df <- readRDS(subset_path)
  primary <- primary[, colnames(df)]
  
  # Filter variants
  misc_df <- data.frame(rowData(primary))
  vars <- misc_df %>%
    dplyr::filter(
      n_cells_conf_detected >= 5,
      strand_correlation > 0.65,
      log10(vmr) > -2,
      mean_coverage >= 5
    ) %>%
    pull(variant)
  
  # Subset to filtered variants
  prim <- primary[vars,]
  rowData(prim)$called_in_primary <- rownames(prim) %in% vars
  
  # Prepare mean heteroplasmy data
  mean_df <- data.frame(
    donor = donor_id,
    mean_heteroplasmy = rowData(prim)$mean * 100,  # convert to %
    variant = rownames(prim)
  )
  
  return(mean_df)
}

# ---- Run for both donors ----

donor1_df <- process_donor(
  output_dir = ".../Donor1_0p5_single/",
  primary_path = ".../variants_donor1_mitoBender.rds",
  subset_path = ".../donor1_subset_mtitoFiltered.rds",
  donor_id = "Donor1"
)
donor1_df


donor0_df <- process_donor(
  output_dir = ".../Donor0_0p5pct_single/",
  primary_path = ".../variants_donor0_mitoBender.rds",
  subset_path = ".../donor0_subset_mtitoFiltered.rds",
  donor_id = "Donor0"
)

donor0_df
# Combine datasets
combined_df <- bind_rows(donor1_df, donor0_df)
head(combined_df)


# Pivot: donors = rows, variants = columns
wide <- combined_df %>%
  dplyr::select(variant, donor, mean_heteroplasmy) %>%
  tidyr::pivot_wider(names_from = variant, values_from = mean_heteroplasmy)

mat <- as.matrix(wide %>% dplyr::select(-donor))
rownames(mat) <- wide$donor
mode(mat) <- "numeric"

# Optional: cap values
mat <- pmin(mat, 100)

# Order rows (donors) by average heteroplasmy
row_order <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)

# Order columns (variants) by maximum heteroplasmy
col_order <- order(apply(mat, 2, max, na.rm = TRUE), decreasing = TRUE)

mat <- mat[row_order, col_order, drop = FALSE]

# Viridis color mapping
rng <- range(mat, na.rm = TRUE)
col_fun <- circlize::colorRamp2(seq(rng[1], rng[2], length.out = 9),
                                viridis::viridis(9))

# Heatmap
p1 <- Heatmap(
  mat,
  name = "Mean heteroplasmy (%)",
  col = col_fun,
  na_col = "grey90",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = rownames(mat),
  column_order = colnames(mat),
  column_names_side = "top"
)
p1
pdf(file="ComplexHeatmap_pseubobulk_heteroplasmy_mitofiltered.pdf", width = 10, height = 3)
p1
dev.off()


# Pivot: variants = rows, donors = columns
wide <- combined_df %>%
  dplyr::select(variant, donor, mean_heteroplasmy) %>%
  tidyr::pivot_wider(names_from = donor, values_from = mean_heteroplasmy)

mat <- as.matrix(wide %>% dplyr::select(-variant))
rownames(mat) <- wide$variant
mode(mat) <- "numeric"

# Optional: cap values
mat <- pmin(mat, 100)

# Order rows (variants) by maximum heteroplasmy
row_order <- order(apply(mat, 1, max, na.rm = TRUE), decreasing = TRUE)

# Order columns (donors) by average heteroplasmy
col_order <- order(colMeans(mat, na.rm = TRUE), decreasing = TRUE)

mat <- mat[row_order, col_order, drop = FALSE]

# Viridis color mapping
rng <- range(mat, na.rm = TRUE)
col_fun <- circlize::colorRamp2(seq(rng[1], rng[2], length.out = 9),
                                viridis::viridis(9))

# Heatmap: donors on columns, variants on rows
p1 <- Heatmap(
  mat,
  name = "Mean heteroplasmy (%)",
  col = col_fun,
  na_col = "grey90",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = rownames(mat),
  column_order = colnames(mat),
  column_names_side = "top"
)
p1
# Save to PDF
pdf(file = "ComplexHeatmap_pseudobulk_heteroplasmy_mitofiltered.pdf", width = 3, height = 6)
p1
dev.off()

