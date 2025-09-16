# Define color schemes for plotting cell types and dataset regions

output_dir <- "/.../Donor1_0p5_single/"

df <- readRDS(file.path(output_dir,file="mitoFilter_df_annotated.rds"))


meta.data <- df@meta.data

umap_coords <- Embeddings(df, "umap")
umap_df <- as.data.frame(umap_coords)
meta.data <- cbind(meta.data, umap_df)
unique(meta.data$possible_Celltype)
celltype_colors <- c(
  "Oligodendrocytes"    = "deepskyblue3",
  "Oligodendrocytes I" = "deepskyblue2",
  "OPCs"                = "mediumpurple1",
  "Microglia"           = "#F294AD",
  "Lymphocytes" = "plum1",
  "Astrocytes"          = "darkseagreen2",
  "Astrocytes I"          = "darkseagreen3",
  "Astrocyte-like neural stem cells"          = "darkseagreen",
  "Glia Cells"         = "skyblue",
  "Vascular Cells"      = "#400135",  
  "Ependymal Cells"     = "#6B2747",  
  "non neuronal mesoderm" = "#6B2671",
  "Projection Nerons (cortical glutameric lineage)"             = "orange1",
  "Excitatory Neurons (pyramidal, glutamatergic)"             = "orange2"
)



p <- ggplot(meta.data, aes(x = UMAP_1, y =UMAP_2, color = possible_Celltype)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = celltype_colors) +
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
        axis.line =element_blank())
p

ggsave(
  filename = file.path(output_dir, "annotated_UMAP_mitoFilter_all_cells.pdf"),
  plot = p,
  width = 7,
  height = 5,
  dpi = 300)


# barplot % cellular compoostion in each Regin

ab <- table(df$possible_Celltype, df$assign)
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

truc = sapply(X = 1:28, FUN = sqrt)
truc2 = sapply(X = 1:28, FUN = function(x) sqrt(x))

sample.name.vec = sapply(X = colnames(ab), FUN = function(x) rep(x, nb.celltype))
sample.name.vec = as.vector(sample.name.vec)
head(sample.name.vec)
df.plot = data.frame(CellType = celltype.name.vec, Sample = sample.name.vec, CellProp = prop.cells.vec)
head(df.plot)
ggplot(data = df.plot, aes(x = Sample, y = CellProp, fill = CellType)) +
  geom_bar(stat = "identity")

celltype_order <- c("Oligodendrocytes",
                    "Oligodendrocytes I",
                    "Glia Cells" ,
                    "OPCs",
                    "Astrocytes",
                    "Astrocytes I",
                    "Astrocyte-like neural stem cells",
                    "Projection Nerons (cortical glutameric lineage)",
                    "Excitatory Neurons (pyramidal, glutamatergic)",
                    "Microglia",
                    "Lymphocytes",
                    "Vascular Cells",  
                    "Ependymal Cells",  
                    "non neuronal mesoderm")
df.plot$CellType <- factor(df.plot$CellType, levels = celltype_order)


assign_order <- c("donor1", "donor0")

bar <- ggplot(data = df.plot, aes(x = Sample, y = CellProp, fill = CellType)) +
  geom_bar(stat = "identity", color = "gray", linewidth = 0.1) +
  scale_fill_manual(values = celltype_colors) +
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

ggsave(plot = bar, filename = "PropCells_assign.pdf", width = 3, height = 5, dpi=300)

