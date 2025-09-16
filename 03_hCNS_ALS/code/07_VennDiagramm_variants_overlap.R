library(Signac)
library(Seurat)
library(SummarizedExperiment)
library(ggplot2)
library(patchwork)
library(readxl)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(VennDiagram)

setwd("...")



#load the called variants 

area1 <- readRDS(file="Region_5/output/primary.rds")
area2 <- readRDS(file="Region_10/output/primary.rds")
area3 <- readRDS(file="Region_28/output/primary.rds")
area4 <- readRDS(file="Region_29/output/primary.rds")


#heteroplasmic

vars_area1 <- data.frame(rowData(area1)) %>%
  filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 5) %>%
  pull(variant)

vars_area2 <- data.frame(rowData(area2)) %>%
  filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 5) %>%
  pull(variant)

vars_area3 <- data.frame(rowData(area3)) %>%
  filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 5) %>%
  pull(variant)

vars_area4 <- data.frame(rowData(area4)) %>%
  filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 5) %>%
  pull(variant)

sets <- list(Region_5 = vars_area1,
             Region_10 = vars_area2,
             Region_28 = vars_area3, 
             Region_29 = vars_area4
)


library(VennDiagram)

# Create the Venn diagram as a grid object (don't save to file automatically)
venn_plot <- venn.diagram(
  x = sets,
  filename = NULL,  # this avoids automatic saving
  lwd = 2,
  lty = 'solid',
  fill = NA,
  col = c("firebrick", "orange2", "dodgerblue3", "plum"),
  cex = 0.6,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.5, 
  cat.fontfamily = "sans",
  cat.default.pos = "outer"
)
venn_plot
# Save manually to PDF
pdf("VennHeteroplsmicc.pdf", width = 8, height = 8)
grid.draw(venn_plot)
dev.off()

