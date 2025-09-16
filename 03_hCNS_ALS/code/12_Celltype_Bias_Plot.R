library(Seurat)
library(reshape2)
library(ggplot2)
library(cowplot)

setwd(".../")


# The workflow shown here is representative for all Regions 
# csv files are provided with this github

'%ni%' <- Negate('%in%')

# Define color schemes for plotting cell types and dataset regions
colors.possible_Celltype <- c(
  "Oligodendrocytes" = "deepskyblue3",
  "Oligodendrocytes I" = "deepskyblue4",
  "Microglia" = "#F294AD",
  "proinflammatory Microglia" = "#BF3459",
  "CD8 T-Cells" = "#73026B",
  "Endothelial Cells" = "#400135",
  "Glia Cells" = "skyblue",
  "Glia Cells I" = "skyblue1",
  "Glia Cells II" = "skyblue2",
  "Glia Cells III" = "skyblue3",
  "inhibitory neurons" = "goldenrod2",
  "GABA inhibitory neurons" = "goldenrod1",
  "excitatory neurons" = "orange1",
  "Astrocytes" = "darkseagreen2",
  "Astrocytes I" = "darkseagreen3",
  "OPCs" = "mediumpurple1",
  "OPCs I" = "mediumpurple2",
  "OPCs II" = "mediumpurple3"
)


# Define output directories for CSV and plot files
dir_csv <- ".../test_output/"
dir_plots <- ".../test_output/"


celltype_mito_combined <- read.csv(("Region_29_possible_Celltype_bias.csv"), row.names = 1) #combination of output of summarized experiment with celltype annotation, rownames = cellbc, colnames = mtDNA + celltype annotation

var_artifactX <- c('X301A>C', 'X302A>C', 'X309C>T','X310T>C', 'X316G>C', 'X3109T>C')

celltype_mito_combined <- celltype_mito_combined %>% 
  dplyr::filter(mut %ni% var_artifactX)

head(celltype_mito_combined)

p <- ggplot(celltype_mito_combined, aes(
  x = -log10(kruskal_pvalue_adj), 
  y = Heteroplasmy_celltype * 100,
  col = possible_Celltype, 
  label = mut
)) +
  geom_point(aes(size = n_cells_celltype, alpha = 0.5)) +
  scale_y_log10(breaks = c(0.01, 0.05, 0.1, 0.2, 0.5) * 100) +
  scale_color_manual(values = colors.possible_Celltype) +
  theme_bw() + 
  theme(panel.border = element_blank()) +
  xlim(0, 40) +   # set x-axis range
  theme(
    axis.title.x = element_text(color = "black", size = 5, face = "bold"),
    axis.title.y = element_text(color = "black", size = 5, face = "bold"),
    axis.text = element_text(size = 5),
    legend.position = "none"
  )

p
ggsave2(p, filename ="Celltype_bias_Region_29.pdf", width = 3, height = 3, limitsize = FALSE)
