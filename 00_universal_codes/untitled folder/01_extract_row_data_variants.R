library(data.table)
library(dplyr)
library(stringr)
library(Signac)
library(SummarizedExperiment)


#Code is based on Lareau et al. Nature Biotechnology 2021 
#see https://github.com/caleblareau/mtscATACpaper_reproducibility


#Set working directory to the folder where you stored the mgatk called variants (from 00_variant_calling.R)
setwd("/../../")


#read in variants
primary <- readRDS(file="variants.rds")
primary

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
  variants_csv_path <- paste0("/../output/ncof", threshold, "_variants.csv")
  misc_df_csv_path <- paste0("/../misc_df_ncof", threshold, ".csv")
  rds_path <- paste0("/../output/variants_ncof", threshold, ".rds")
  
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


#additional for homoplasmic variants

df_prim1 <- df_prim %>%
    dplyr::filter(n_cells_conf_detected >= 5 & 
                  strand_correlation > 0.65 & 
                  vmr < 0.01 & 
                  mean_coverage >= 5) %>%
    dplyr::pull(variant)

write.csv(df_prim1, "/../homoplasmic_variants.csv")
df_prim2 <- subset(primary, variant %in% df_prim1)

misc_df <- data.frame(rowData(df_prim2))
write.csv(misc_df, "/../misc_df_homoplasmic.csv", row.names = TRUE)

saveRDS(df_prim2, "/../variants_homoplasmic.rds")


#additionally for ncof 1 but not filtered on any other QC matrix (will be used for VMR plot only!)

df_prim1 <- df_prim %>% 
  dplyr::filter(n_cells_conf_detected >1) %>% 
  dplyr::pull(variant) 
df_prim1
write.csv(df_prim1, "/.../output/ncof1_only_variants.csv")
df_prim2 <- subset(primary, variant %in% df_prim1)

misc_df <- data.frame(rowData(df_prim2))
write.csv(misc_df, "/.../output/misc_df_ncof1_only.csv", row.names = TRUE)

saveRDS(df_prim2, "/.../output/variants_ncof1_only.rds")
