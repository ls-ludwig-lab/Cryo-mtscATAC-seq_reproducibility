library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readxl)


setwd(".../test_output")

#choose the Region of interest
# The workflow shown here is representative for all Regions 
# (Region_5, Region_10, Region_28, Region_29).
#

R_5 <- read.csv("Microglia_Clones_R5_meta_1.csv")
R_10 <- read.csv("Microglia_Clones_R10_meta_1.csv")
R_28 <- read.csv("Microglia_Clones_R28_meta_1.csv")
R_29 <- read.csv("Microglia_Clones_R29_meta_1.csv")



variant_colors <- c(
  "wt" = "grey90",   
  "X8920G_A" = "orangered2",                    # not defined, assigned black
  "X1900A_G" = "orangered4",              # 1900A>G
  "X2075T_C" = "orangered3",              # 2075T>C
  "X1345G_A" = "grey90",              # 1345G>A (for the UMAP here also grey as I dont show it in here)
  "X15346G_A" = "darkorange3",            # 15346G>A
  "X4247T_C" = "darkorange2",                   # not defined
  "X14046A_G" = "deepskyblue4",                  # not defined
  "X16131T_C" = "deepskyblue3",                  # not defined
  "X9445G_A" = "dodgerblue3",                   # not defined
  "X14305G_A" = "dodgerblue4",            # 14305G>A
  "X14046A_G;X9445G_A" = "dodgerblue",       # new color
  "X7561T_C" = "hotpink2",                   # not defined
  "X3705G_A" = "hotpink3",                   # not defined
  "X6872A_G" = "maroon4",                 # 6872A>G
  "X15135G_A;X7561T_C" = "maroon3",       # new color
  "X2778T_C" = "maroon",                  # 2778T>C
  "X12869G_A" = "purple",                 # 12869G>A
  "X3705G_A;X7278T_C" = "hotpink4",        # new color
  "X12083T_C" = "maroon2",                # 12083T>C
  "X15135G_A;X1999A_G;X7561T_C;X8612T_C" = "palevioletred4", # new color
  "X7278T_C" = "purple3",                 # 7278T>C
  "X1661A_G" = "palevioletred3",          # 1661A>G
  "X1999A_G;X7561T_C" = "palevioletred2"         # new color
)


# sanity check: do all Clone labels have a color?
setdiff(unique(R_10$Clone), names(variant_colors))
# if this prints anything, add those names to variant_colors
p1 <- ggplot(R_29, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Clone), size = 1) +
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
ggsave(plot = p1, filename = "Region_R_29_clones.pdf", width = 2, height = 2, dpi=300)