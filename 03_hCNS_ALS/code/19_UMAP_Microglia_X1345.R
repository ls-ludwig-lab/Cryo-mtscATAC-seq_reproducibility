library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readxl)
setwd(".../test_output")


R_5 <- read.csv("Microglia_Clones_R5_meta_1.csv")
R_10 <- read.csv("Microglia_Clones_R10_meta_1.csv")
R_28 <- read.csv("Microglia_Clones_R28_meta_1.csv")
R_29 <- read.csv("Microglia_Clones_R29_meta_1.csv")

combined <- rbind(R_5, R_10, R_28, R_29)
dim(combined)

variant_colors <- c(
  "wt" = "grey90",              
  "1345G>A" = "turquoise1"     
)


# sanity check: do all Clone labels have a color?
setdiff(unique(R_10$Clone), names(variant_colors))
unique(combined$X1345G_A)

combined$.fg <- combined$X1345G_A == "X1345G_A"
combined_ord <- combined[order(combined$.fg), ]


# if this prints anything, add those names to variant_colors
p1 <- ggplot(combined_ord, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = X1345G_A), size = 1) +
  scale_color_manual(values = variant_colors) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank()
  )
p1
ggsave(plot = p1, filename = "X1345_clone.pdf", width = 2, height = 2, dpi=300)