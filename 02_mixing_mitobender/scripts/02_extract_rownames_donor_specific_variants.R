library(data.table)
library(dplyr)
library(stringr)
library(Signac)
library(SummarizedExperiment)

#subsetting to donor specific variants

datasets <- c("mitoBender")

# Loop over each dataset
for (dataset in datasets) {
  
  # Construct file paths
  input_rds_path <- paste0(".../data/", dataset,
                            "/primary.rds")
  output_csv_path <- paste0(".../data/", dataset,
                             "/variants_donor_specific.csv")
  output_rds_path <- paste0("/.../data/", dataset,
                            "/variants_donor_specific.rds")
  
  primary
  # Load data
  primary <- readRDS(file = input_rds_path)
  
  # Create a data.frame from rowData
  df_prim <- data.frame(rowData(primary))
  head(df_prim)
  # Filter: n_cells_conf_detected > 1 and pull 'variant'
  df_prim1 <- df_prim %>%
    filter(n_cells_conf_detected > ncol(primary) *0.1 & n_cells_conf_detected < ncol(primary)*0.9 &
             strand_correlation > 0.65 &
             vmr > 0.1 &
             mean > 0.1 & mean < 0.9) %>%
    pull(variant)
  df_prim1
  
  # Write filtered variants to CSV
  write.csv(df_prim1, file = output_csv_path, row.names = FALSE)
  
  # Subset primary to keep only selected variants
  df_prim2 <- subset(primary, variant %in% df_prim1)
  
  # Save the subsetted object
  saveRDS(df_prim2, file = output_rds_path)
  
  # (Optional) print progress
  message("Finished processing ", dataset)
}
