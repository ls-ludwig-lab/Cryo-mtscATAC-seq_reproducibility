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

setwd("...")

df <- readRDS("integrated_annotated_1.rds")

df <- read.csv("meta.data.csv")
unique(df$possible_Celltype)
unique(df$dataset)
#[1] "Region_5"  "Region_10" "Region_28" "Region_29"


colors.Celltype <- c("deepskyblue3" = "Oligodendrocytes",
                     "deepskyblue4" = "Oligodendrocytes I",
                     "#F294AD" = "Microglia",
                     "#BF3459" = "proinflammatory Microglia",
                     "#73026B" = "CD8 T-Cells",
                     "#400135" = "Endothelial Cells",
                     "skyblue" = "Glia Cells",
                     "skyblue1" = "Glia Cells I",
                     "skyblue2" = "Glia Cells II",
                     "skyblue3" = "Glia Cells III",
                     "goldenrod2" = "exitatory neurons",
                     "goldenrod1" = "GABAergic interneurons",
                     "orange1" = "excitatory neurons I",
                     "darkseagreen2" = "Astrocytes",
                     "darkseagreen3" = "Astrocytes I",
                     "mediumpurple1" = "OPCs",
                     "mediumpurple2" = "OPCs I",
                     "mediumpurple3" = "OPCs II")

colors.Celltype <- setNames(names(colors.Celltype), colors.Celltype)

p <- ggplot(df, aes(x = UMAP_1, y =UMAP_2, color = possible_Celltype)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = c(colors.Celltype)) +
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

ggsave(plot = p, filename = "UAMP_Celltype_nolegend.pdf", width = 4, height = 4, dpi=300)



# barplot % cellular compoostion in each Regin

ab <- table(df$possible_Celltype, df$Region)
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

truc = sapply(X = 1:72, FUN = sqrt)
truc2 = sapply(X = 1:72, FUN = function(x) sqrt(x))

sample.name.vec = sapply(X = colnames(ab), FUN = function(x) rep(x, nb.celltype))
sample.name.vec = as.vector(sample.name.vec)
head(sample.name.vec)
df.plot = data.frame(CellType = celltype.name.vec, Sample = sample.name.vec, CellProp = prop.cells.vec)
head(df.plot)
ggplot(data = df.plot, aes(x = Sample, y = CellProp, fill = CellType)) +
  geom_bar(stat = "identity")


df.plot$Sample <- factor(df.plot$Sample, levels = c("Region_5", "Region_10", "Region_28", "Region_29"))

celltype_order <- c("Microglia", "proinflammatory Microglia", "CD8 T-Cells", "Endothelial Cells",
                    "Astrocytes", "Astrocytes I",
                    "Glia Cells", "Glia Cells I", "Glia Cells II", "Glia Cells III",
                    "OPCs", "OPCs I", "OPCs II",
                    "Oligodendrocytes", "Oligodendrocytes I",
                    "exitatory neurons", "GABAergic interneurons", "excitatory neurons I")
df.plot$CellType <- factor(df.plot$CellType, levels = celltype_order)

bar <- ggplot(data = df.plot, aes(x = Sample, y = CellProp, fill = CellType)) +
  geom_bar(stat = "identity", color = "gray", linewidth = 0.1) +
  scale_fill_manual(values = colors.Celltype) +
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

ggsave(plot = bar, filename = "PropCells_Region.pdf", width = 3, height = 5, dpi=300)



p <- ggplot(meta.data, aes(x = UMAP_1, y =UMAP_2, color = mtDNA_depth)) +
  geom_point(size = 0.1) +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_color_continuous(type = "viridis",limits = c(5, 100)) +
  theme(aspect.ratio=1/1,
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.1))
p
ggsave(plot = p, filename = "UMAP_mtDNA100.pdf", width = 4, height = 4, dpi=300)


