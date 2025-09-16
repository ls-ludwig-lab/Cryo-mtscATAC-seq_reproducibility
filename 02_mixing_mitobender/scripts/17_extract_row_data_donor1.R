library(data.table)
library(dplyr)
library(stringr)
library(Signac)
library(SummarizedExperiment)

input_dir <- (".../Donor1_0p5_single/")
output_dir <- (".../Donor1_0p5_single")

primary <- readRDS(file.path(input_dir, "variants_donor1.rds"))

df_prim <- data.frame(rowData(primary))


# Define the different thresholds for n_cells_conf_detected
thresholds <- c(1, 3, 5)

# Loop over each threshold
for (threshold in thresholds) {
  # Filter the data based on the current threshold
  df_prim_filtered <- df_prim %>%
    dplyr::filter(n_cells_conf_detected >= threshold & 
                  strand_correlation > 0.65 & 
                  vmr > 0.01 & 
                  mean_coverage >= 5) %>%
    dplyr::pull(variant)
  
  # Define file paths with the current threshold included in the filename
  variants_csv_path <- paste0(file.path(output_dir), threshold, "_variants.csv")
  misc_df_csv_path <- paste0(file.path(output_dir), threshold, ".csv")
  rds_path <- paste0(file.path(output_dir), threshold, ".rds")
  
  # Write the filtered variants to a CSV file
  write.csv(df_prim_filtered, variants_csv_path)
  
  # Subset the primary dataset based on the filtered variants
  df_prim_subset <- subset(primary, variant %in% df_prim_filtered)
  
  # Convert the subsetted data to a data frame and save as CSV
  misc_df <- data.frame(rowData(df_prim_subset))
  write.csv(misc_df, misc_df_csv_path, row.names = TRUE)
  
  # Save the subsetted data as an RDS file
  saveRDS(df_prim_subset, rds_path)
}


