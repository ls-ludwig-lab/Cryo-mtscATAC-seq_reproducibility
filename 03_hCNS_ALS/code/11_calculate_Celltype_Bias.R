library(Seurat)
library(SummarizedExperiment)  
library(Matrix)                
library(ggplot2)
library(cowplot)               
library(reshape2)              
library(dplyr)                 
library(tibble)



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

region.colors <- c("Region_5" = "firebrick",
                   "Region_10" = "orange2",
                   "Region_28" = "dodgerblue3",
                   "Region_29" = "plum")


# Define the dataset names and corresponding variant RDS files
datasets <- c("Region_5", "Region_10", "Region_28", "Region_29")
variants_list <- c("R5", "R10", "R28", "R29")

# Define output directories for CSV and plot files
dir_csv <- "..."
dir_plots <- "..."

variant_dir <- ".../output"
variants_list <- file.path(variant_dir, c("Region_5/output/variants_ncof5.rds",
                                          "Region_10/output/variants_ncof5.rds",
                                          "Region_28/output/variants_ncof5.rds",
                                          "Region_29/output/variants_ncof5.rds"))




# Define a wrapper function to process each dataset
process_dataset <- function(dataset_name, variants_list) {
  
  # Load the integrated annotated Seurat object and subset for the current dataset
  d <- readRDS("integrated_annotated.rds")
  df <- subset(d, subset = dataset == dataset_name)
  df <- RenameCells(df, new.names = gsub(paste0("^", dataset_name, "_"), "", Cells(df)))
  
  # Load variant RDS file and align with the Seurat object
  variants <- readRDS(variants_list)
  variants <- variants[, colnames(df)]
  
  # Add allele frequency assay to the Seurat object
  assay_data <- assay(variants, "allele_frequency")
  df[["alleles"]] <- CreateAssayObject(counts = assay_data)
  
  # Prepare metadata and allele matrix for downstream calculations
  metadata_df <- df@meta.data[, "possible_Celltype", drop = FALSE]
  alleles_counts <- df@assays[["alleles"]]@counts
  alleles_counts_df <- as.data.frame(t(alleles_counts[, rownames(metadata_df)]))
  colnames(alleles_counts_df) <- ifelse(grepl("^[0-9]", colnames(alleles_counts_df)),
                                        paste0("X", colnames(alleles_counts_df)),
                                        colnames(alleles_counts_df))
  
  # Merge metadata and alleles into one dataframe
  mtDNA_vars_celltype <- cbind(metadata_df, alleles_counts_df)
  mtDNA_vars_celltype$patient <- dataset_name
  
  # Compute Kruskal-Wallis p-values for each variant vs. cell type presence
  type <- "possible_Celltype"
  bias <- data.frame(possible_Celltype = NA, mut = NA, bias = NA)
  
  for(possible_Celltype in unique(mtDNA_vars_celltype$possible_Celltype)) {
    y <- as.numeric(mtDNA_vars_celltype$possible_Celltype == possible_Celltype)
    for(mut in grep("^X", colnames(mtDNA_vars_celltype), value = TRUE)) {
      tmp <- data.frame(possible_Celltype, mut,
                        bias = kruskal.test(mtDNA_vars_celltype[[mut]] ~ y)$p.value)
      bias <- rbind(tmp, bias)
    }
  }
  
  # Adjust p-values using Benjamini-Hochberg correction
  bias$kruskal_pvalue_adj <- p.adjust(bias$bias, method = "BH")
  bias <- na.omit(bias)
  
  # Calculate mean heteroplasmy per mutation per cell type
  celltypes <- unique(mtDNA_vars_celltype$possible_Celltype)
  mutation_cols <- grep("^X", colnames(mtDNA_vars_celltype), value = TRUE)
  
  heteroplasmy_matrix <- do.call(cbind, lapply(celltypes, function(ct) {
    df_ct <- mtDNA_vars_celltype[mtDNA_vars_celltype$possible_Celltype == ct, mutation_cols, drop = FALSE]
    colMeans(df_ct)
  }))
  colnames(heteroplasmy_matrix) <- celltypes
  
  heteroplasmy_celltype <- melt(heteroplasmy_matrix, value.name = "Heteroplasmy_celltype")
  heteroplasmy_celltype$possible_Celltype <- heteroplasmy_celltype$Var2
  heteroplasmy_celltype$mut <- heteroplasmy_celltype$Var1
  heteroplasmy_celltype$Var1 <- NULL
  heteroplasmy_celltype$Var2 <- NULL
  
  # Count number of cells with heteroplasmy > 0.2 per mutation and cell type
  n_cells_matrix <- do.call(cbind, lapply(celltypes, function(ct) {
    df_ct <- mtDNA_vars_celltype[mtDNA_vars_celltype$possible_Celltype == ct, mutation_cols, drop = FALSE]
    colSums(df_ct > 0.2)
  }))
  colnames(n_cells_matrix) <- celltypes
  
  n_cells_per_celltype <- melt(n_cells_matrix, value.name = "n_cells_celltype")
  n_cells_per_celltype$possible_Celltype <- n_cells_per_celltype$Var2
  n_cells_per_celltype$mut <- n_cells_per_celltype$Var1
  n_cells_per_celltype$Var1 <- NULL
  n_cells_per_celltype$Var2 <- NULL
  
  # Merge p-values, heteroplasmy and cell count data into one summary dataframe
  complete_bias_df <- merge(bias, merge(heteroplasmy_celltype, n_cells_per_celltype, by=c("possible_Celltype","mut")), by=c("possible_Celltype","mut"))
  
  # Export the summary statistics as CSV for this dataset
  write.csv(complete_bias_df, paste0(dir_csv, dataset_name, "_possible_Celltype_bias.csv"))
  
  print(paste0("Finished processing ", dataset_name))
}


# Run the wrapper function for all datasets using a loop
for (i in seq_along(datasets)) {
  process_dataset(datasets[i], variants_list[i])
}

