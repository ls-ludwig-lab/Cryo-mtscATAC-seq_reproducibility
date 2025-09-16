library(Seurat)
library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(ggplot2)


setwd(".../called_variants/")

d <- readRDS(file=".../output/integrated_annotated.rds")

  # Subset the Seurat object for the current dataset
df <- subset(d, subset = Donor == "patient1_primary")

  # Modify the cell names to remove the prefix
new_cell_names <- gsub(pattern = "^patient1_primary_single_", replacement = "", x = Cells(df))
df <- RenameCells(object = df, new.names = new_cell_names)
  
  # Assuming variantsR is a SummarizedExperiment object
variantsR <- readRDS(".../called_variants/primary.rds")
  
# Adding the variants as assay data to the Seurat object
variantsR <- variantsR[, colnames(df)]
assay_data <- assay(variantsR, "allele_frequency")
  
  # Add the new assay to the Seurat object
df[["alleles"]] <- CreateAssayObject(counts = assay_data)
DefaultAssay(df) <- "alleles"
  

# Loop through each variant to create the plots
to_subset <- c("12865A>G", "12008G>A")
  
# Extract the data and convert it into a data frame
variants_df <- as.data.frame(assay(variantsR, "allele_frequency"))
  
# Now you can filter the variants data frame
variants1 <- variants_df %>% 
    rownames_to_column(var = "value") %>%  # Ensure the rownames (variant names) are included in the data
    dplyr::filter(value %in% to_subset)  # Apply the filter to get the variants of interest
  
for (i in 1:nrow(variants1)) {
    column <- variants1$column[i]
    value <- variants1$value[i]
    
    # Extract the allele frequency values
    allele_values <- GetAssayData(df, assay = "alleles", slot = "data")[value, ]
    
    # Get UMAP coordinates and prepare for plotting
    umap_coords <- Embeddings(df, "umap")
    umap_df <- as.data.frame(umap_coords)
    umap_df$allele_values <- allele_values
    
    # Order the UMAP dataframe by allele values
    umap_df <- umap_df[order(umap_df$allele_values), ]
    df@meta.data <- df@meta.data[rownames(umap_df), ]
    
    # Create and save the plot
    p2 <- ggplot(data = umap_df, aes(x = UMAP_1, y = UMAP_2, color = allele_values)) +
      geom_point(size = 1) +
      scale_color_gradientn(colors = c("grey", "darkred"), limits = c(0, 1)) +
      theme_void() +
      labs(title = paste("Variant:", value))
    p2
    # Save the plot with a dynamic filename
    p2_filename <- paste0("patient1_primary_", value, ".pdf")
    ggsave(plot = p2, filename = p2_filename, width = 6, height = 5, dpi = 300, limitsize = FALSE)
  }
