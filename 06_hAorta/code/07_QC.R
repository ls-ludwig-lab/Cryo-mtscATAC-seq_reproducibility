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


setwd(".../output/")


df <- readRDS(file="preprocessed.rds")

meta.data <- df@meta.data
library(ggplot2)

# List of y variables
y_vars <- c("nCount_peaks", "nFeature_peaks", "mtDNA_depth", "TSS.enrichment")  # Replace with actual variable names

meta.data <- meta.data %>%
  mutate(Celltype = case_when(
    peaks_snn_res.0.5 == "0" ~ "SMC",
    peaks_snn_res.0.5 == "1" ~ "SMC I",
    peaks_snn_res.0.5 == "2" ~ "Fibroblasts",
    peaks_snn_res.0.5 == "3" ~ "Macrophages",
    peaks_snn_res.0.5 == "4" ~ "NK/T Cells",
    peaks_snn_res.0.5 == "5" ~ "Endothelial Cells",
    TRUE ~ "Unknown"
  ))
head(meta.data)


clustering.color <- c("SMC" = "deepskyblue3",
                      "SMC I" = "skyblue2",
                      "Fibroblasts" = "darkseagreen3",
                      "Macrophages" = "plum2",
                      "NK/T Cells" = "lightsalmon",
                      "Endothelial Cells" = "brown2")


# Loop through y_vars to create and save plots
for (y_var in y_vars) {
  # Create the plot for the current y_var
  p <- ggplot(meta.data, aes_string(x = "Celltype", y = y_var, fill = "Celltype")) +
    geom_boxplot(trim = TRUE) +
    scale_fill_manual(values = clustering.color) +
    theme_classic() +
    scale_y_log10() +  # Apply logarithmic scale to y-axis
    theme(
      axis.title.y = element_text(size = 20),
      axis.line.x = element_line(),
      axis.text.y = element_text(size = 15),
      axis.title.x = element_blank(),
      legend.position = "none"
    )
  print(p)
  # Save the plot
  ggsave(paste0("boxplot_", y_var, ".pdf"), plot = p, width = 5, height = 5)
}

