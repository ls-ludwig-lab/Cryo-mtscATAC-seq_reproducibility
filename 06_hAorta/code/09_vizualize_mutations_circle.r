# Load libraries
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(cowplot)

# Define output directory
output_dir <- ".../output/"


# 1. Gene annotation function
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
  "Control-Region" = "lightblue4", "tRNA" = "black", "rRNA" = "mediumaquamarine", 
  "Non-Coding" = "sienna4", "ND1" = "magenta", "ND2" = "mediumblue", 
  "CO1" = "olivedrab", "CO2" = "orange2", "ATP8" = "orchid4", "ATP6" = "red3",
  "CO3" = "royalblue2", "ND3" = "palegreen4", "ND4L" = "grey0", 
  "ND4" = "pink4", "ND5" = "yellow4", "ND6" = "steelblue4", "CYB" ="tan", "black"
)


# 3. Input dataframes and labeling

"%ni%" <- Negate("%in%")
setwd(".../output/")
# Import per-sample calls
sample1 <-  as.data.frame(rowData(readRDS("variants_ncof1_only.rds")))

df <- readRDS(file="preprocessed.rds")
sample1 <- sample1[,colnames(df)]

sample1

# 4. Combine data
combined_df <- sample1
combined_df$gene <- addgenelabel(combined_df$position)
# 5. Gene annotation base
gene_points <- data.frame(pos = 1:16569)
gene_points$color_code <- addgenelabel(gene_points$pos)
gene_points$newAF <- 0.00002
# 6. Align gene ring over grey ring
y_base <- 0.00002
gene_ring <- gene_points %>%
  group_by(color_code) %>%
  summarise(start = min(pos), end = max(pos), .groups = "drop") %>%
  mutate(
    ymin = y_base - 0.000003,
    ymax = y_base + 0.000003
  )
# 7. Gene label positions
gene_labels <- gene_ring %>%
  mutate(mid = (start + end) / 2,
         label = as.character(color_code))
# 8. Tick labels
tick_labels <- seq(0, 16000, by = 4000)
tick_df <- data.frame(position = tick_labels, mean = 0.00002)
# 9. Reference rings
ref_vals <- c(0.00001, 0.0001, 0.001, 0.01)
ref_df <- expand.grid(position = seq(0, 16569, length.out = 100),
                      mean = ref_vals)
# 10. Sample colors



# 11. Density curve for variant distribution
density_data <- density(combined_df$position, bw = 200, from = 0, to = 16569, n = 512)
max_y <- 0.01  # outer radial limit
min_y <- 0.000003  # inner radial limit
scale_factor <- (max_y - min_y) * 0.2  # 20% of radial range
y_offset <- max_y + scale_factor * 0.5  # center above plot
density_df <- data.frame(
  position = density_data$x,
  density = density_data$y
)
density_df$y <- y_offset + density_df$density * scale_factor / max(density_df$density)




# 12. Plot
head(combined_df)
p1 <- ggplot(combined_df, aes(x = position, y = mean)) +
  # Heteroplasmy points
  geom_point(size = 0.5) +
  # Grey circular base
  geom_point(data = gene_points[seq(1, nrow(gene_points), length.out = 200),],
             aes(x = pos, y = newAF),
             size = 0.6, color = "gray80", inherit.aes = FALSE) +
  # Gene blocks
  geom_rect(data = gene_ring,
            aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = color_code),
            inherit.aes = FALSE,
            color = "grey40", alpha = 0.7) +
  # Gene labels (currently hidden)
  geom_text(data = gene_labels,
            aes(x = mid, y = 0.00002, label = label),
            size = 0,
            angle = (gene_labels$mid / 16569) * 360,
            hjust = 0, vjust = -0.5,
            inherit.aes = FALSE) +
  # Tick labels (currently hidden)
  geom_text(data = tick_df,
            aes(x = position, y = mean, label = position),
            color = "black", size = 0,
            angle = (tick_df$position / 16569) * 360,
            hjust = 0, vjust = -0.5,
            inherit.aes = FALSE) +
  # Reference rings
  geom_line(data = ref_df,
            aes(x = position, y = mean, group = mean),
            color = "gray60", size = 0.2, inherit.aes = FALSE) +
  annotate("text", x = 0, y = ref_vals, label = ref_vals,
           size = 4, hjust = -0.1, vjust = 0.5) +
  # Variant density curve
  geom_ribbon(data = density_df,
              aes(x = position, ymin = y_offset, ymax = y),
              fill = "black", alpha = 0.3, inherit.aes = FALSE) +
  # Style
  scale_fill_manual(values = gene_colors, name = "Gene") +
  # D-loop highlight
  #geom_rect(
    #data = data.frame(xmin = c(16024, 1), xmax = c(16569, 576), ymin = min_y, ymax = y_offset + scale_factor),
    #aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    #fill = "lightcyan", alpha = 0.5, inherit.aes = FALSE) + 

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
# 13. Show plot
p1
cowplot::ggsave2(p1, filename = ".../output/circos_ncof1.pdf", 
                 height = 5, width = 5)

