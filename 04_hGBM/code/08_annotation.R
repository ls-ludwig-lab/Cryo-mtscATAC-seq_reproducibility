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


setwd("../output")

df <- readRDS(file="integrated.rds")

head(df@meta.data)

Idents(df) <- "ATAC_snn_res.0.8"

DimPlot(df, label = TRUE)
df

FeaturePlot(df, "GRIN1", order = TRUE)

df <- RenameIdents(df,"0" = "malignant_other",
                   "1" = "malignant_other",
                   "2" = "Microglia",
                   "3" = "AC-MEs-like",
                   "4" = "OPC-like",
                   "5" = "Oligodendrocytes",
                   "6" = "DLX-GAD-high Neuronal Cells",
                   "7" = "malignant_other",
                   "8" = "Neuronal Cells",
                   "9" = "AC-like",
                   "10" = "Neuronal Cells",
                   "11" = "AC-like",
                   "12" = "Neuronal Cells",
                   "13" = "Endothelial Cells",
                   "14" = "Microglia I",
                   "15" = "NPC1-OPC-like",
                   "16" = "malignant_other",
                   "17" = "OPC-like",
                   "18" = "Pericytes",
                   "19" = "AC-like",
                   "20" = "AC-like",
                   "21" = "T Cells")

# Save cluster identities as a new column in meta.data
df$possible_Celltype_broader <- df@active.ident

umap_coords <- Embeddings(df, "umap")
umap_df <- as.data.frame(umap_coords)
df@meta.data <- cbind(df@meta.data, umap_df)
meta.data <- df@meta.data
meta.data

colnames(df@meta.data)

saveRDS(df, "integrated_annotated_1.rds")
write.csv(meta.data, "meta.data_annotated_1.csv")

unique(meta.data$possible_Celltype_broader)

Celltype.color <- c("malignant_other" = "#829BD4",
                    "Microglia" = "#D9488B",
                    "NPC1-OPC-like" = "lightslateblue",
                    "OPC-like" = "dodgerblue", 
                    "Oligodendrocytes" =  "dodgerblue4",
                    "Neuronal Cells" = "orange",
                    "AC-MEs-like" = "#A7BAF2",
                    "DLX-GAD-high Neuronal Cells" = "sienna2",
                    "OPC-like" = "dodgerblue3",
                    "AC-like" = "paleturquoise3",
                    "Endothelial Cells" = "#400135",
                    "Microglia I" = "#F294AD",
                    "Pericytes" = "#73026B",
                    "T Cells" ="#8C233F")

meta.data <- df@meta.data
p <- ggplot(meta.data, aes(x = UMAP_1, y =UMAP_2, color = possible_Celltype_broader)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = Celltype.color) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        #legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        axis.line = element_blank())
p

ggsave(plot = p, filename = "UAMP_annotated.pdf", width = 7, height = 5, dpi=300)

