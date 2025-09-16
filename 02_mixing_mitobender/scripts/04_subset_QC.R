library(ggplot2)
library(Seurat)

# Define input files and corresponding output labels
datasets <- list(file = "unfiltered_0p5pct.rds", label = "0p5pct")

# Set output directory
setwd("...")

# Define donor colors
donor.color <- c(
  "donor0" = "#5FCDD9", "donor1" = "#F00D09", "doublet" = "orange2", "unassigned" = "#ccd5ae")
orig.ident <- c(
  "orig.ident" = "#ccd5ae")

# Define variables to plot
y_vars <- c("nCount_peaks", "nFeature_peaks", "mtDNA_depth", "TSS.enrichment")

# Loop over each dataset
for (dataset in datasets) {
  # Read the input file
  df <- readRDS(paste0(".../rds_files/", dataset$file))
  
  # Apply filtering
  df1 <- subset(df, subset = nCount_peaks > 1000 & TSS.enrichment > 1.5 & mtDNA_depth >= 5)
  
  # Extract metadata
  meta.data <- df1@meta.data
  
  # Generate and save plots
  for (y_var in y_vars) {
    p <- ggplot(meta.data, aes_string(x = "orig.ident", y = y_var, fill = "orig.ident")) +
      geom_boxplot(trim = TRUE) +
      scale_fill_manual(values = donor.color) +
      theme_classic() +
      scale_y_log10() +
      theme(
        axis.title.y = element_text(size = 20),
        axis.line.x = element_line(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank()
      )
    print(p)
    ggsave(paste0("boxplot_", dataset$label, "_", y_var, ".pdf"), plot = p, width = 8, height = 6)
  }
  
  # Save the filtered object
  saveRDS(df1, paste0("filtered_", dataset$label, ".rds"))
}
