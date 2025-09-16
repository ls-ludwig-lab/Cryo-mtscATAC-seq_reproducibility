# load libraries
library(Signac)
library(Seurat)
library(ggplot2)
library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(data.table)
library(tibble)
library(tidyr)

set.seed(257)
"%ni%" <- Negate("%in%")

# Load Objects
NB <- readRDS("path_to/Harmony_Integrated_Object.rds")
GN <- readRDS("path_to/GN01.rds")

# ============================================
# Filtering Variant Calling Output (ncof >=1)
# ============================================
## Define your dataset names
datasets <- c("NB01", "NB02", "GN01")

# Define base input and output paths
input_base_path <- "..."
output_base_path <- ".../Variant_Calling_Objects"

# Loop over each dataset
for (dataset in datasets) {
  # Construct input file path
  input_rds_path <- file.path(input_base_path, paste0("Variant_Calling_", dataset, ".rds"))
  
  # Construct output file paths
  output_csv_path <- file.path(output_base_path, dataset, "output", "ncof1_only_variants.csv")
  output_rds_path <- file.path(output_base_path, dataset, "output", "ncof1_only_variants.rds")
  
  # Load data
  SE <- readRDS(file = input_rds_path)
  
  # Create a data.frame from rowData
  df_prim <- data.frame(rowData(SE))
  
  # Filter: n_cells_conf_detected > 1 and pull ‘variant’
  df_prim1 <- df_prim %>%
    dplyr::filter(n_cells_conf_detected >= 1) %>%
    dplyr::filter(mean_coverage >= 5) %>%
    dplyr::filter(vmr > 0.01) %>%
    dplyr::filter(strand_correlation > 0.65) %>%
    pull(variant)
  
  # Write filtered variants to CSV
  write.csv(df_prim1, file = output_csv_path, row.names = FALSE)
  
  # Subset SE to keep only selected variants
  df_prim2 <- SE[which(df_prim$variant %in% df_prim1), ]
  
  # Save the subsetted object
  saveRDS(df_prim2, file = output_rds_path)
}


# ============================================
# Filtering Variant Calling Output (ncof >=5)
# ============================================
## Define your dataset names
datasets <- c("NB01", "NB02", "GN01")

# Define base input and output paths
input_base_path <- "..."
output_base_path <- ".../Variant_Calling_Objects"

# Loop over each dataset
for (dataset in datasets) {
  # Construct input file path
  input_rds_path <- file.path(input_base_path, paste0("Variant_Calling_", dataset, ".rds"))
  
  # Construct output file paths
  output_csv_path <- file.path(output_base_path, dataset, "output", "ncof5_only_variants.csv")
  output_rds_path <- file.path(output_base_path, dataset, "output", "ncof5_only_variants.rds")
  
  # Load data
  SE <- readRDS(file = input_rds_path)
  
  # Create a data.frame from rowData
  df_prim <- data.frame(rowData(SE))
  
  # Filter: n_cells_conf_detected > 1 and pull ‘variant’
  df_prim1 <- df_prim %>%
    dplyr::filter(n_cells_conf_detected >= 5) %>%
    dplyr::filter(mean_coverage >= 5) %>%
    dplyr::filter(vmr > 0.01) %>%
    dplyr::filter(strand_correlation > 0.65) %>%
    pull(variant)
  
  # Write filtered variants to CSV
  write.csv(df_prim1, file = output_csv_path, row.names = FALSE)
  
  # Subset SE to keep only selected variants
  df_prim2 <- SE[which(df_prim$variant %in% df_prim1), ]
  
  # Save the subsetted object
  saveRDS(df_prim2, file = output_rds_path)
}

# ============================================
# Plot VMR - Strand Concordance Plot
# ============================================
# Define your sample names
to_subset <- c("NB01", "NB02", "GN01")

# Loop over each sample
for (sample in to_subset) {
  
  # Load variant file
  variant_file <- file.path(
    ".../Variant_Calling_Objects",
    sample, "output", "ncof1_only_variants.rds"
  )
  
  variants <- readRDS(variant_file)
  
  # Extract metadata
  misc_df <- data.frame(rowData(variants))
  
  # Create the strand concordance vs log10 VMR plot
  p <- ggplot(misc_df, aes(x = strand_correlation, y = log10(vmr),
                           color = log10(vmr) > -2 & strand_correlation > 0.65)) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c("black", "firebrick")) +
    labs(color = "HQ", x = "Strand concordance", y = "log VMR") +
    theme_classic() +
    geom_vline(xintercept = 0.6, linetype = 2) +
    geom_hline(yintercept = -2, linetype = 2) +
    theme(
      legend.position = "none",
      text = element_text(colour = "black", size = 15),
      axis.line = element_line(color = "black", size = 0.8),
      axis.text.x = element_text(colour = "black", size = 15, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 15),
      axis.ticks.y = element_line(color = "black", size = 0.8),
      axis.ticks.x = element_line(color = "black", size = 0)
    )
  
  # Save the plot as a PDF
  output_plot <- file.path(
    ".../Variant_Calling_Objects",
    sample, "output", paste0("StrandConcordance_VMR_", sample, ".pdf")
  )
  
  ggsave(plot = p, filename = output_plot, width = 7, height = 4.5, dpi = 300)
  
  message("Finished processing: ", sample)
}

# ============================================
# Mitochondrial DNA Mutational Distribution
# ============================================
## Import ncof5 filtered variant calling output
NB01_mgatk_ncof5 <- readRDS(".../Variant_Calling_Objects/NB01/output/ncof5_only_variants.rds")
NB02_mgatk_ncof5 <- readRDS(".../Variant_Calling_Objects/B2R2/output/ncof5_only_variants.rds")
GN01_mgatk_ncof5 <- readRDS(".../Variant_Calling_Objects/B150/output/ncof5_only_variants.rds")

NB01_mgatk_ncof5 <- data.frame(rowData(NB01_mgatk_ncof5))
NB02_mgatk_ncof5 <- data.frame(rowData(NB02_mgatk_ncof5))
GN01_mgatk_ncof5 <- data.frame(rowData(GN01_mgatk_ncof5))


## 1. Gene annotation function
addgenelabel <- function(bp) {
  breaks <- c(0, 576, 648, 1602, 1671, 3230, 3305, 3307, 4263, 4332, 4401, 4402, 4470,
              5512, 5580, 5587, 5656, 5657, 5730, 5826, 5892, 5904, 7446, 7515, 7518, 7586,
              8270, 8295, 8365, 8366, 8573, 9208, 9991, 10059, 10405, 10470, 10767, 12138,
              12207, 12266, 12337, 14149, 14674, 14743, 14747, 15888, 15954, 15956, 16024, 17000)
  
  labels <- c("Control-Region", "tRNA", "rRNA", "tRNA", "rRNA", "tRNA", "Non-Coding", "ND1", "tRNA", "tRNA",
              "Non-Coding", "tRNA", "ND2", "tRNA", "Non-Coding", "tRNA", "Non-Coding", "tRNA", "tRNA", "tRNA",
              "Non-Coding", "CO1", "tRNA", "Non-Coding", "tRNA", "CO2", "Non-Coding", "tRNA", "Non-Coding",
              "ATP8", "ATP6", "CO3", "tRNA", "ND3", "tRNA", "ND4L", "ND4", "tRNA", "tRNA", "tRNA", "ND5", "ND6",
              "tRNA", "Non-Coding", "CYB", "tRNA", "Non-Coding", "tRNA", "Control-Region")
  
  cut(bp, breaks = breaks, labels = labels, include.lowest = TRUE, right = FALSE)
}

# 2. Gene colors
gene_colors <- c(
  "tRNA" = "#934258",
  "rRNA" = "#F1656E",
  "ND1" = "#E6550D",
  "ND2" = "#FD8D3C",
  "ND3" = "#FDAE6B",
  "ND4" = "#FEE6CE",
  "ND5" = "#31A354",
  "ND6" = "#74C476",
  "CO1" = "#3182BD",
  "CO2" = "#6BAED6",
  "CO3" = "#9ECAE1",
  "ATP6" = "#756BB1",
  "ATP8" = "#BCBDDC",
  "CYB" = "#E7BA52",
  "Non-Coding" = "#CCCCCC",
  "Control-Region" = "#969696"
)

# 3. Gene annotation base
gene_points <- data.frame(pos = 1:16569)
gene_points$color_code <- addgenelabel(gene_points$pos)
gene_points$newAF <- 0.00002

y_base <- 0.00002
gene_ring <- gene_points %>%
  group_by(color_code) %>%
  summarise(start = min(pos), end = max(pos), .groups = "drop") %>%
  mutate(
    ymin = y_base - 0.00001,
    ymax = y_base + 0.000005
  )

gene_labels <- gene_ring %>%
  mutate(mid = (start + end) / 2,
         label = as.character(color_code))

tick_labels <- seq(0, 16000, by = 4000)
tick_df <- data.frame(position = tick_labels, mean = 0.00002)

ref_vals <- c(0.00001, 0.0001, 0.001, 0.07)
ref_df <- expand.grid(position = seq(0, 16569, length.out = 100),
                      mean = ref_vals)

# 4. Sample info
df_list <- list(GN01_mgatk_ncof5, NB01_mgatk_ncof5, NB02_mgatk_ncof5)  
sample_names <- c("GN01", "NB01", "NB02")
sample_colors <- c("#756BB1", "#76C3A9", "#355C4C")
names(sample_colors) <- sample_names

# 5. Loop over samples to create and save plots
for (i in seq_along(df_list)) {
  
  df <- df_list[[i]]
  df$sample <- sample_names[i]  # ensure the sample column matches
  
  df$gene <- addgenelabel(df$position)
  
  # Density curve for variant distribution
  density_data <- density(df$position, bw = 200, from = 0, to = 16569, n = 512)
  
  max_y <- 0.07
  min_y <- 0.000003
  scale_factor <- (max_y - min_y) * 0.2
  y_offset <- max_y + scale_factor * 0.5
  
  density_df <- data.frame(
    position = density_data$x,
    density = density_data$y
  )
  density_df$y <- y_offset + density_df$density * scale_factor / max(density_df$density)
  
  # Plot
  p <- ggplot(df, aes(x = position, y = mean, color = sample)) +
    
    geom_point(size = 3) +
    
    geom_point(data = gene_points[seq(1, nrow(gene_points), length.out = 200),],
               aes(x = pos, y = newAF),
               size = 0.6, color = "gray80", inherit.aes = FALSE) +
    
    geom_rect(data = gene_ring,
              aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = color_code),
              inherit.aes = FALSE,
              color = "grey40", alpha = 0.7) +
    
    geom_text(data = gene_labels,
              aes(x = mid, y = 0.00002, label = label),
              size = 0,
              angle = (gene_labels$mid / 16569) * 360,
              hjust = 0, vjust = -0.5,
              inherit.aes = FALSE) +
    
    geom_text(data = tick_df,
              aes(x = position, y = mean, label = position),
              color = "black", size = 0,
              angle = (tick_df$position / 16569) * 360,
              hjust = 0, vjust = -0.5,
              inherit.aes = FALSE) +
    
    geom_line(data = ref_df,
              aes(x = position, y = mean, group = mean),
              color = "gray60", size = 0.75, inherit.aes = FALSE) +
    annotate("text", x = 0, y = ref_vals, label = ref_vals,
             size = 5, hjust = -0.1, vjust = 0.5) +
    
    geom_ribbon(data = density_df,
                aes(x = position, ymin = y_offset, ymax = y),
                fill = "black", alpha = 0.3, inherit.aes = FALSE) +
    
    geom_point(size = 2) +
    
    scale_color_manual(values = sample_colors, name = "Sample") +
    scale_fill_manual(values = gene_colors, name = "Gene") +
    coord_polar(direction = 1) +
    scale_y_log10(breaks = ref_vals, limits = c(min_y, y_offset + scale_factor)) +
    labs(x = "", y = "Heteroplasmy (mean)", title = sample_names[i]) +
    theme_minimal(base_size = 6) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.position = "right"
    )
  
  # Save PDF
  pdf_path <- paste0(".../Variant_Calling_Objects/Mutation_Distribution_", sample_names[i], ".pdf")
  pdf(pdf_path, width = 8, height = 8)
  print(p)
  dev.off()
}




