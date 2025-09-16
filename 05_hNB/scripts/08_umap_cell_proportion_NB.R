# Load libraries
library(Seurat)      
library(ggplot2)     
library(dplyr)       
library(ggalluvial)  

# Load Seurat Object and Metadata
seurat <- readRDS("path_to/Harmony_Integrated_Object.rds")  # Load Seurat object
metadata <- read.table("Metadata_NB.tsv", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

# Add Sample Information to Seurat Object
seurat@meta.data <- seurat@meta.data %>%
  mutate(Dataset = case_when(
    Sample == "NB01" ~ "Primary",       # Label NB01 as "Primary"
    Sample == "NB02" ~ "Second-look",   # Label NB02 as "Second-look"
    TRUE ~ NA_character_                 # All others as NA
  ))

# ============================
# UMAP Plot
# ============================
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),  # Lower truncation for axis
  trunc_upper = unit(3, "cm")    # Upper truncation for axis
)

ggplot(metadata, aes(umap_1, umap_2, color = Manual_Annotation)) +
  geom_point(size = 0.2, alpha = 1) + # Plot each cell
  theme_classic() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  scale_color_manual(values = c("#A8D7E3", "#616887", "#b37bb3", "#F1656E")) + # Color by cell type
  theme(
    aspect.ratio = 1/1,
    # legend.position = "none",        # Optional: hide legend
    legend.background = element_rect(fill = "white"),
    legend.key.size = unit(0.5, 'cm'),
    legend.text = element_text(size=10),
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(arrow = arrow(angle = 15, length = unit(.5, "cm"), type = "closed")),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14), 
    axis.title = element_text(hjust = 0),
    axis.line.x = element_line(linetype = "solid", size = 0.5, arrow = grid::arrow(angle=20, length=unit(0.075, "inches"), type = "closed")),
    axis.line.y = element_line(linetype = "solid", size = 0.5, arrow = grid::arrow(angle=20, length=unit(0.075, "inches"), type = "closed"))
  ) +
  guides(x = axis, y = axis) +
  guides(colour = guide_legend(override.aes = list(size=4)))  # Make legend points larger

# ============================
# Alluvial Plot: Temporal Shifts in Cell Type Composition
# ============================

# Ensure factor levels are consistent
seurat$Manual_Annotation <- factor(seurat$Manual_Annotation,
                                   levels = c("Endothelial", "Fibroblast", "Immune", "Neuroendocrine"))

# Create abundance table: counts per cell type per dataset
abundance <- table(seurat$Manual_Annotation, seurat$Dataset)

# Convert counts to proportions
propotion <- sweep(x = abundance, MARGIN = 2, STATS = colSums(abundance), FUN = "/")
prop.cells.vec <- as.vector(propotion)

# Prepare vectors for plotting
celltype.name.vec <- rep(rownames(abundance), 2)  # Repeat cell types for two datasets
sample.name.vec <- sapply(X = colnames(abundance), FUN = function(x) rep(x, 4))  # Repeat dataset names
sample.name.vec <- as.vector(sample.name.vec)

# Combine into data frame for plotting
df.plot <- data.frame(
  CellType = celltype.name.vec,
  Sample = sample.name.vec,
  CellProp = prop.cells.vec
)

df.plot$CellType <- factor(df.plot$CellType, levels = c("Endothelial", "Fibroblast", "Immune", "Neuroendocrine"))
df.plot$AlluviaID <- df.plot$CellType  # ID to track cell types across samples

# Plot alluvial diagram
ggplot(df.plot,
       aes(x = Sample, stratum = CellType, alluvium = AlluviaID,
           y = CellProp, fill = CellType, label = CellType)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", color = "darkgray") +  # Flow between samples
  geom_stratum(width = 0.5, color = "black") +  # Strata for each cell type
  theme_classic() +
  scale_fill_manual(values = c("#A8D7E3", "#616887", "#b37bb3", "#F1656E")) +
  scale_x_discrete() +
  labs(y = "Proportion", x = "Condition") +
  theme(
    text = element_text(colour = "black", size = 18), 
    axis.line = element_line(color = "black", size = 1),
    axis.text.x = element_text(colour = "black", size = 18),
    axis.text.y = element_text(colour = "black", size = 18),
    axis.ticks = element_line(color = "black", size = 1)
  )


