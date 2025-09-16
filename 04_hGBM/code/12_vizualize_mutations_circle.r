library(SummarizedExperiment)
library(BuenColors)
library(reshape2)
library(stringr)
library(ggplot2)
library(cowplot)

"%ni%" <- Negate("%in%")
setwd(".../output") 

output_dir <- ("/.../output")

# Load your samples
samples <- list(
  patient1_primary = readRDS("../ncof1_variants.rds"),
  patient1_recu = readRDS("../ncof1_variants.rds"),
  patient2_primary = readRDS("../ncof1_variants.rds"),
  patient2_recu = readRDS("../ncof1_variants.rds")
)

# Gene annotation function
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

# Define gene color map
gene_colors <- c(
  "Control-Region" = "lightblue4", "tRNA" = "black", "rRNA" = "mediumaquamarine", 
  "Non-Coding" = "sienna4", "ND1" = "magenta", "ND2" = "mediumblue", 
  "CO1" = "olivedrab", "CO2" = "orange2", "ATP8" = "orchid4", "ATP6" = "red3",
  "CO3" = "royalblue2", "ND3" = "palegreen4", "ND4L" = "grey0", 
  "ND4" = "pink4", "ND5" = "yellow4", "ND6" = "steelblue4", "CYB" = "tan", "black"
)

dataset.color <- c("N20_1201" ="#5B79BE",
                   "N20_4400" = "#8C233F", 
                   "N20_5563" = "#A7BAF2",
                   "N21_1209" = "#c1121f")

# Precompute static elements
gene_points <- data.frame(pos = 1:16569)
gene_points$color_code <- addgenelabel(gene_points$pos)
gene_points$newAF <- 0.00002

y_base <- 0.00002
gene_ring <- gene_points %>%
  group_by(color_code) %>%
  summarise(start = min(pos), end = max(pos), .groups = "drop") %>%
  mutate(ymin = y_base - 0.000003, ymax = y_base + 0.000003)

ref_vals <- c(0.00001, 0.0001, 0.001, 0.01)
ref_df <- expand.grid(position = seq(0, 16569, length.out = 100), mean = ref_vals)

max_y <- 0.01
min_y <- 0.000003
scale_factor <- (max_y - min_y) * 0.2
y_offset <- max_y + scale_factor * 0.5

# Main plotting loop
for (sample_name in names(samples)) {
  se <- samples[[sample_name]]
  df <- data.frame(position = rowData(se)$position,
                   nucleotide = rowData(se)$nucleotide,
                   mean = rowData(se)$mean)
  df$gene <- addgenelabel(df$position)
  
  # Density curve
  density_data <- density(df$position, bw = 200, from = 0, to = 16569, n = 512)
  density_df <- data.frame(position = density_data$x, density = density_data$y)
  density_df$y <- y_offset + density_df$density * scale_factor / max(density_df$density)
  
  # Generate plot
  p <- ggplot(df, aes(x = position, y = mean, colour = sample_name)) +
    geom_point(size = 0.8) +
    scale_color_manual(values = c(dataset.color)) +
    geom_point(data = gene_points[seq(1, nrow(gene_points), length.out = 200),],
               aes(x = pos, y = newAF),
               size = 0.6, color = "gray80", inherit.aes = FALSE) +
    geom_rect(data = gene_ring,
              aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = color_code),
              inherit.aes = FALSE,
              color = "grey40", alpha = 0.7) +
    geom_line(data = ref_df,
              aes(x = position, y = mean, group = mean),
              color = "gray60", size = 0.2, inherit.aes = FALSE) +
    annotate("text", x = 0, y = ref_vals, label = ref_vals,
             size = 4, hjust = -0.1, vjust = 0.5) +
    geom_ribbon(data = density_df,
                aes(x = position, ymin = y_offset, ymax = y),
                fill = "black", alpha = 0.3, inherit.aes = FALSE) +
    scale_fill_manual(values = gene_colors, name = "Gene") +
    coord_polar(direction = 1) +
    scale_y_log10(breaks = ref_vals, limits = c(min_y, y_offset + scale_factor)) +
    labs(x = "", y = "Heteroplasmy (mean)") +
    theme_minimal(base_size = 6) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.position = "none"
    )
  
  # Save plot
  ggsave2(
    p,
    filename = file.path(output_dir, paste0("circos_", sample_name, ".pdf")),
    height = 5, width = 5
  )
}

