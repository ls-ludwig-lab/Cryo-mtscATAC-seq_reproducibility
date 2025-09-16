library(Seurat)
library(ggplot2)
library(SummarizedExperiment)
datasets <- c("patient1_primary", "patient1_recu", "patient2_primary","patient2_recu")  
unique_idents <- c("patient1_primary_single_", "patient2_recu_single_", "patient1_recu_single_", "patient2_primary_single_")  

# Define dataset info
dataset_name <- c("patient1_primary", "patient2_recu", "patient1_recu", "patient2_primary")  
datasets_cell <-c("patient1_primary_single_", "patient2_recu_single_", "patient1_recu_single_", "patient2_primary_single_")
dataset_direction <- c("patient1_primary/output/", 
                       "patient2_recu/output/", 
                       "patient1_recu/output/", 
                       "patient2_recu/output/")

# Set working directory to root
setwd("...")

# Load integrated Seurat object
df <- readRDS("integration/output/integrated_annotated.rds")

# Define detection thresholds
thresholds <- c("_ncof1")

for (i in seq_along(dataset_name)) {
  current_name <- dataset_name[i]
  current_prefix <- datasets_cell[i]
  current_dir <- dataset_direction[i]
  
  message(paste("Processing dataset:", current_name))
  
  # Subset the integrated Seurat object — NO renaming!
  df_subset <- subset(df, subset = dataset == current_name)
  
  # Loop over thresholds
  for (threshold in thresholds) {
    rds_path <- paste0(current_dir,current_name, thresholds, "_only_variants", ".rds")
    
    if (file.exists(rds_path)) {
      variants <- readRDS(rds_path)
      
      # Add prefix to colnames of the variants matrix to match integrated Seurat object
      colnames(variants) <- paste0(current_prefix, colnames(variants))
      
      # Only keep cells that are present in df_subset
      common_cells <- intersect(colnames(df_subset), colnames(variants))
      variants <- variants[, common_cells]
      
      # Defensive check
      if (length(common_cells) == 0) {
        message("No matching cells after prefixing. Skipping.")
        next
      }
      
      misc_df <- data.frame(rowData(variants))
      assay_data <- assay(variants, "allele_frequency")
      
      # Subset the assay data too
      assay_data <- assay_data[, common_cells]
      
      # Add assay to Seurat object
      df_subset[["alleles"]] <- CreateAssayObject(counts = assay_data)
      DefaultAssay(df_subset) <- "alleles"
      
      # Plot
      misc_df <- misc_df[is.finite(log10(misc_df$vmr)) & !is.na(misc_df$strand_correlation), ]
      p <- ggplot(misc_df, aes(x = strand_correlation, y = log10(vmr), 
                               color = log10(vmr) > -2 & strand_correlation > 0.65)) +
        geom_point(size = 0.4) + 
        scale_color_manual(values = c("black", "firebrick")) +
        labs(color = "HQ", x = "Strand concordance", y = "log VMR") +
        geom_vline(xintercept = 0.6, linetype = 2) +
        geom_hline(yintercept = -2, linetype = 2) +
        theme_minimal() +
        theme(legend.position = "none",axis.line = element_line(color = "black", linewidth = 0.2),
              axis.ticks = element_line(color = "black", linewidth = 0.2))
      
      # Save
      output_path <- paste0(current_dir, "VariantPlot_", current_name, "_ncof_1only_", ".pdf")
      ggsave(plot = p, filename = output_path, width = 7, height = 5, dpi = 300)
      
    } else {
      message(paste("File not found:", rds_path))
    }
  }
}
