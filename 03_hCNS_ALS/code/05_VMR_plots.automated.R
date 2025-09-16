library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(readxl)
library(RColorBrewer)
library(SummarizedExperiment)



setwd("/data/cephfs-1/work/projects/ludwig-mtscatac-als/ALS/2024_Region/test_output")

#Here load the variants, which were subsetted by only filtering to ncof1 (this rds file will only be used to generate this plot (procudes similar plot to signac vignette))

df <- readRDS(file="integrated_annotated.rds")

to_subset <- c("Region_5", "Region_10", "Region_28", "Region_29")

for (sample in to_subset) {
  message("Processing: ", sample)
  
  # Subset the Seurat object
  df_sub <- subset(df, subset = dataset == sample)
  
  # Modify cell names to remove the prefix
  new_cell_names <- gsub(pattern = paste0("^", sample, "_"), replacement = "", x = Cells(df_sub))
  df_sub <- RenameCells(object = df_sub, new.names = new_cell_names)
  
  # Load corresponding variants file (corrected path)
  variant_file <- paste0(".../output/", sample, "/output/ncof1_only_variants.rds")
  
  if (!file.exists(variant_file)) {
    warning("Variant file not found for ", sample, ": ", variant_file)
    next  # Skip this sample if the file does not exist
  }
  
  variants <- readRDS(variant_file)
  
  # Subset variants to match cells in df_sub
  variants <- variants[, colnames(df_sub)]
  misc_df <- data.frame(rowData(variants))
  
  # Add allele frequency assay
  assay_data <- assay(variants, "allele_frequency")
  df_sub[["alleles"]] <- CreateAssayObject(counts = assay_data)
  DefaultAssay(df_sub) <- "alleles"
  
  # Generate plot
  p <- ggplot(misc_df, aes(x = strand_correlation, y = log10(vmr), color = log10(vmr) > -2 & strand_correlation > 0.65)) +
    geom_point(size = 0.2) + 
    scale_color_manual(values = c("black", "firebrick")) +
    labs(color = "HQ", x = "Strand concordance", y = "log VMR") +
    theme_classic() +
    geom_vline(xintercept = 0.6, linetype = 2) +
    geom_hline(yintercept = -2, linetype = 2) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(-5, 0)) + 
    theme(legend.position = "none")
  
  print(p)
  
  # Save plot with correct output path
  output_plot <- paste0("/.../output/", sample, "/output/VariantPlot_ncof1.pdf")
  ggsave(plot = p, filename = output_plot, width = 5, height = 3, dpi = 300)
  
  message("Finished processing: ", sample)
}
