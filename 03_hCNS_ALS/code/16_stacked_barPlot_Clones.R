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



# barplot % cellular compoostion in each Regin


ab <- table(combined$Clone, combined$Region)
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

truc = sapply(X = 1:88, FUN = sqrt)
truc2 = sapply(X = 1:88, FUN = function(x) sqrt(x))

sample.name.vec = sapply(X = colnames(ab), FUN = function(x) rep(x, nb.celltype))
sample.name.vec = as.vector(sample.name.vec)
head(sample.name.vec)
df.plot = data.frame(CellType = celltype.name.vec, Sample = sample.name.vec, CellProp = prop.cells.vec)
head(df.plot)
ggplot(data = df.plot, aes(x = Sample, y = CellProp, fill = CellType)) +
  geom_bar(stat = "identity")


df.plot$Sample <- factor(df.plot$Sample, levels = c("Region_5", "Region_10", "Region_28", "Region_29"))

variant_colors <- c(
  "wt" = "grey90",   
  "X8920G_A" = "orangered2",                    # not defined, assigned black
  "X1900A_G" = "orangered4",              # 1900A>G
  "X2075T_C" = "orangered3",              # 2075T>C
  "X1345G_A" = "turquoise1",              # 1345G>A
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

head(df.plot)
# Make sure "wt" is first in the CellType factor

variant_levels <- c(
  "X1345G_A",
  setdiff(unique(df.plot$CellType), c("X1345G_A", "wt")),
  "wt"
)

# Set the factor levels
df.plot$CellType <- factor(df.plot$CellType, levels = variant_levels)

bar <- ggplot(data = df.plot, aes(x = Sample, y = CellProp, fill = CellType)) +
  geom_bar(stat = "identity", color = "gray", linewidth = 0.1) +
  scale_fill_manual(values = variant_colors) +
  ylim(0, 0.2) +
  theme_bw() +
  theme(
    #axis.line=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position="none",
    panel.background=element_blank(),
    panel.border=element_blank(),
    #panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank())
bar
ggsave(plot = bar, filename = "PropClones_Region_clones.pdf", width = 3, height = 5, dpi=300)