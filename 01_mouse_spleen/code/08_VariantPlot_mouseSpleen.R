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

#setwd to the folder where you have stored the mgatk output files (for mutation calling)
#to reprpoduce the figure from Lareau et al. Nature Biotcechnology 2021 cells were filtered only on ncof1 (see )
setwd(".../FF_10min/output/")
setwd(".../FF_30min/output/")
setwd(".../Cells/output/")

#Loop through the different variants data frames 

# Define the different n_cells_conf_detected thresholds and corresponding file names
thresholds <- c(1)

# Loop over each threshold
for (threshold in thresholds) {
  # Construct the file path for the RDS file based on the threshold
  rds_path <- paste0("ncof", thresholds, "_only_variants", ".rds")
  
  # Check if the file exists before attempting to read it
  if (file.exists(rds_path)) {
    # Read the variants RDS file
    variants <- readRDS(rds_path)
    
    # Subset to HQ cells that exist in the df
    misc_df <- data.frame(rowData(variants))
    
    # Extract the assay data, assuming "allele_frequency" is the assay name
    assay_data <- assay(variants, "allele_frequency")
    
    # Create the plot
    p <- ggplot(misc_df, aes(x = strand_correlation, y = log10(vmr), color = log10(vmr) > -2 & strand_correlation > 0.65)) +
      geom_point(size = 0.4) + scale_color_manual(values = c("black", "firebrick")) +
      ylim(-5,0) +
      xlim(-1,1)+
      labs(color = "HQ", x = "Strand concordance", y = "log VMR") +
      theme(legend.position = "bottom",
            plot.background = element_blank(),        # Remove background color for the entire plot
            panel.background = element_blank(),       # Remove background color for the plot area
            panel.grid = element_blank(),             # Remove grid lines
            panel.border = element_blank())+ 
      geom_vline(xintercept = 0.6, linetype = 2) +
      geom_hline(yintercept = -2, linetype = 2) + 
      theme(legend.position = "none",axis.line = element_line(color = "black", linewidth = 0.2),
            axis.ticks = element_line(color = "black", linewidth = 0.2))
    print(p)
    # Save the plot with a filename based on the threshold
    plot_filename <- paste0("VariantPlot_ncof", threshold, ".pdf")
    ggsave(plot = p, filename = plot_filename, width = 7, height = 5, dpi = 300)
    
  } else {
    message(paste("File not found:", rds_path))
  }
}

