library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(readxl)
library(RColorBrewer)
library(dplyr)
library(tidyr)

setwd(".../integration/output")

df <- readRDS(file="integrated_annotated.rds")

DimPlot(df, label = TRUE)


df$possible_Celltype <- df@active.ident

umap_coords <- Embeddings(df, "umap")
umap_df <- as.data.frame(umap_coords)
df@meta.data <- cbind(df@meta.data, umap_df)
meta.data <- df@meta.data
head(meta.data)

DimPlot(df, group.by = "possible_Celltype")

unique(df$possible_Celltype)

#[1] Dendritic Cell            Myeloid Cell              CD8 T Cell                Marginal Zone B Cell     
#[5] Erythroid Progenitor Cell CD4 T Cell                other                     Macrophage               
#[9] B Cell                    NK/ T Cell                NK Cell                  

cell_type_colors <- c(              # Blue
  "Myeloid Cell" = "lightpink2",      # Orange
  "CD8 T Cell" = "dodgerblue3",                  # Green
  "Monocyte/ Macrophage" = "plum2",        # Red
  "Marginal Zone B Cell" = "darkseagreen3",        # Purple
  "Erythroid Progenitor Cell" = "coral2",                     # Pink
  "other" = "wheat2",                  # Gray
  "Macrophage" = "plum3",                  # Olive
  "unresolved II" = "navajowhite",               # Cyan
  "B Cell" = "darkseagreen2",                      # Light Blue
  "unresolved I" = "wheat3",                # Peach
  "NK/ T Cell" = "dodgerblue",                  # Light Green
  "NK Cell" = "cornflowerblue",                     # Light Red
  "CD4 T Cell" = "deepskyblue3",                  # Lavender
  "Dendritic Cell" = "lightpink", # Light Pink
  "Monocyte" = "plum"                     # Light Brown
)



cell_type_order <- c(
  # Myeloid lineage
  "Myeloid Cell",
  "Monocyte",                      
  "Monocyte/ Macrophage",         
  "Macrophage",                   
  "Dendritic Cell",               
  "Plasmacytoid Dendritic Cell",  
  "Neutrophil / Monocytes",       
  
  # Erythroid lineage
  "Erythroid Progenitor Cell",               
  
  # Lymphoid lineage
  "NK Cell",                      
  "NK/ T Cell",                   
  "CD8 T Cell",                   
  "CD4 T Cell",                   
  
  # B cell lineage
  "B Cell",                       
  "Marginal Zone B Cell",        
  
  # Unresolved clusters
  "other"
)

# Set the factor levels in the biologically-informed order
meta.data$possible_Celltype <- factor(meta.data$possible_Celltype, levels = cell_type_order)


p <- ggplot(meta.data, aes(x = UMAP_1, y =UMAP_2, color = possible_Celltype)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = cell_type_colors) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        axis.line = element_blank())
p

ggsave(plot = p, filename = "UAMP_Celltype.pdf", width = 6, height = 6, dpi=300)



# barplot % cellular compoostion in each Regin

ab <- table(df$possible_Celltype, df$dataset)
ab
nb.samples = ncol(ab)
nb.celltype = nrow(ab)

prop <- sweep(x = ab, MARGIN = 2, STATS = colSums(ab), FUN ="/")
prop
prop.cells.vec <- as.vector(prop)
head(prop.cells.vec)
length(prop.cells.vec)

celltype.name.vec <- rep(rownames(ab), nb.samples)


head(celltype.name.vec)
length(celltype.name.vec)


sample.name.vec = c() #Initialization
for(sample.name in colnames(ab)){
  sample.name.vec = c(sample.name.vec, rep(sample.name, nb.celltype))
}
length(sample.name.vec)
length(sample.name.vec)

truc = sapply(X = 1:33, FUN = sqrt)
truc2 = sapply(X = 1:33, FUN = function(x) sqrt(x))

sample.name.vec = sapply(X = colnames(ab), FUN = function(x) rep(x, nb.celltype))
sample.name.vec = as.vector(sample.name.vec)
head(sample.name.vec)
df.plot = data.frame(CellType = celltype.name.vec, Sample = sample.name.vec, CellProp = prop.cells.vec)
head(df.plot)
ggplot(data = df.plot, aes(x = Sample, y = CellProp, fill = CellType)) +
  geom_bar(stat = "identity")


df.plot$CellType <- factor(df.plot$CellType, levels = cell_type_order)

bar <- ggplot(data = df.plot, aes(x = Sample, y = CellProp, fill = CellType)) +
  geom_bar(stat = "identity", color = "gray", linewidth = 0.1) +
  scale_fill_manual(values =cell_type_colors ) +
  theme(
    #axis.line=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position="none",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank())
bar

ggsave(plot = bar, filename = "PropCells.pdf", width = 5, height = 5, dpi=300)


