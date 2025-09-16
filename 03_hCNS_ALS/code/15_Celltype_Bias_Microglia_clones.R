library(dplyr)
library(readr)

# Data availability:
# Raw sequencing data (FASTQ files) have been deposited at the 
# European Genome-phenome Archive (EGA) under controlled access 
# [accession number to be provided upon acceptance].
# Processed data required to reproduce this analysis, including 
# fragment files and mgatk outputs for all Regions, are available 
# at the Gene Expression Omnibus (GEO) [accession number to be provided].



setwd("...")

'%ni%' <- Negate('%in%')


# List of CSV files
datasets <- c("Region_5_possible_Celltype_bias.csv",
              "Region_10_possible_Celltype_bias.csv",
              "Region_28_possible_Celltype_bias.csv",
              "Region_29_possible_Celltype_bias.csv")

# Known artifacts to filter out
var_artifactX <- c('X301A>C', 'X302A>C', 'X309C>T', 'X310T>C', 'X316G>C', 'X3109T>C', "X189A>G")

# Initialize an empty list to collect filtered results
all_variants_list <- list()

for (file in datasets) {
  # Read data
  df <- read.csv(file, row.names = 1)
  
  # Filter out known artifacts
  df <- df %>% 
    dplyr::filter(!(mut %in% var_artifactX))
  
  # Filter for Microglia with bias ≤ 0.05 and at least 4 cells
  subset_df <- df %>%
    dplyr::filter(grepl("Microglia", possible_Celltype),
                  bias <= 0.05,
                  n_cells_celltype > 3)
  
  # Add dataset name to trace the origin
  subset_df$dataset <- file
  
  # Store in list
  all_variants_list[[file]] <- subset_df
}

# Combine all filtered data into one dataframe
combined_variants <- bind_rows(all_variants_list)

# Optional: Get just the unique variants
variants <- unique(combined_variants$mut)
variants
# Output
head(combined_variants)
variants

write.csv(combined_variants, "celltype_mito_combined.csv")