# Load required libraries
library(Seurat)      
library(ggplot2)     
library(ggh4x)       


# Load Object and Metadata
seurat <- readRDS("path_to/GN01.rds")
metadata <- read.table("Metadata_GN.tsv", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

# ============================
# UMAP Plot
# ============================
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

# Custom color palette
celltype_colors <- c(
  "#e05b02", "#A8D7E3", "#616887", "#95BBA1",
  "#658375", "#704D6C", "#F9C071"
)

ggplot(metadata, aes(UMAP_1, UMAP_2, color = Annotation_Manual)) +
  geom_point(size = 0.2, alpha = 1) +
  scale_color_manual(values = celltype_colors) +
  theme_classic() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(
    aspect.ratio = 1,
    legend.background = element_rect(fill = "white"),
    legend.key.size = unit(0.5, 'cm'),
    legend.text = element_text(size = 10),
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks   = element_blank(),
    axis.line.x  = element_line(linetype = "solid", size = 0.5,
                                arrow = grid::arrow(angle = 20, length = unit(0.075, "inches"), type = "closed")),
    axis.line.y  = element_line(linetype = "solid", size = 0.5,
                                arrow = grid::arrow(angle = 20, length = unit(0.075, "inches"), type = "closed")),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title   = element_text(hjust = 0)
  ) +
  guides(x = axis, y = axis, colour = guide_legend(override.aes = list(size = 4)))






