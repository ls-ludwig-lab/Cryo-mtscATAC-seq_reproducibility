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

setwd(".../output/")


#Loop through the different variants data frames 

df<- readRDS(file="preprocessed.rds")

# Read the variants RDS file
variants <- readRDS("variants_ncof1_only.rds")
    
# Subset to HQ cells that exist in the df
variants <- variants[, colnames(df)]
misc_df <- data.frame(rowData(variants))
    
# Extract the assay data, assuming "allele_frequency" is the assay name
assay_data <- assay(variants, "allele_frequency")
    
# Add the new assay to the Seurat object
df[["alleles"]] <- CreateAssayObject(counts = assay_data)
DefaultAssay(df) <- "alleles"
    
# Create the plot
p <- ggplot(misc_df, aes(x = strand_correlation, y = log10(vmr), color = log10(vmr) > -2 & strand_correlation > 0.65)) +
      geom_point(size = 0.4) + scale_color_manual(values = c("black", "firebrick")) +
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
    ggsave(plot = p, filename = "VariantPlot.pdf", width = 7, height = 5, dpi = 300)
