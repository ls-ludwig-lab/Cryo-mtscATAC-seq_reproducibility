library(dplyr)
library(SummarizedExperiment)

# Set input and output directories
output_dir <- ".../output_plots"
# 1. Set input directory first
input_dir <- ".../output_rds/"

# 2. List files ending in _1.rds (ncof1 variants)
rds_files_1 <- list.files(input_dir, pattern = "_1\\.rds$", full.names = TRUE)

# 3. Read the files
rds_list <- lapply(rds_files_1, readRDS)

# 4. Name the list by the base filename
names(rds_list) <- basename(rds_files_1)

summary_df <- data.frame(
  file = names(rds_list),
  n_rows = sapply(rds_list, function(x) nrow(rowData(x)))
)

print(summary_df)

# add meta data 

#Age

summary_df <- summary_df %>%
  mutate(
    project = case_when(
      grepl("^ALS", file) ~ "ALS",
      grepl("^GBM", file) ~ "GBM",
      grepl("^GN01|^NB", file) ~ "pediatric neoplasm",
      TRUE ~ "Unknown"
    ),
    sample = case_when(
      grepl("ALS_PFC_1", file) ~ "ALS_PFC",
      grepl("ALS_1MNC_1", file) ~ "ALS_1MNC",
      grepl("ALS_MO_1", file) ~ "ALS_MO",
      grepl("ALS_CSC_1", file) ~ "ALS_CSC",
      grepl("GBM_Patient1_primary_1", file) ~ "Patient1_primary",
      grepl("GBM_Patient_1_recurrence_1", file) ~ "Patient1_recurrence",
      grepl("GBM_Patient_2_primary_1", file) ~ "Patient2_primary",
      grepl("GBM_Patient_2_recurrence_1", file) ~ "Patient2_recurrence",
      grepl("GN01_1", file) ~ "GN01",
      grepl("NB01_1", file) ~ "NB01",
      grepl("NB02_1", file) ~ "NB02",
      TRUE ~ NA_character_
    ),
    tissue = case_when(
      project == "ALS" ~ "human brain",
      project == "GBM" ~ "human brain cancer",
      project == "pediatric neoplasm" ~ "pediatric neoplasm",
      TRUE ~ "Unknown"
    ),
    Age = case_when(
      project == "ALS" ~ "80 years",
      sample == "Patient1_primary" ~ "71 years",
      sample == "Patient1_recurrence" ~ "72 years",
      sample == "Patient2_primary" ~ "61 years",
      sample == "Patient2_recurrence" ~ "61 years",
      sample == "GN01" ~ "6,7 years",
      sample == "NB01" ~ "1,8 years",
      sample == "NB02" ~ "2,5 years",
    ),
    Biology = case_when(
      sample == "ALS_PFC" ~ "prefrontal cortex",
      sample == "ALS_1MNC" ~ "1 MN cortex",
      sample == "ALS_MO" ~ "Medulla Oblongata",
      sample == "ALS_CSC" ~ "Cervical Spinal Cord",
      sample == "Patient1_primary" ~ "primary tumor",
      sample == "Patient1_recurrence" ~ "recurrence tumor",
      sample == "Patient2_primary" ~ "primary tumor",
      sample == "Patient2_recurrence" ~ "recurrence tumor",
      sample == "GN01" ~ "Ganglioneuroma bening tumor",
      sample == "NB01" ~ "Neuroblastoma first look",
      sample == "NB02" ~ "Neuroblastoma second look",
    ),
    Patient = case_when(
      project == "ALS" ~ "Patient 2",
      sample == "Patient1_primary" ~ "Patient 3",
      sample == "Patient1_recurrence" ~ "Patient 3",
      sample == "Patient2_primary" ~ "Patient 4",
      sample == "Patient2_recurrence" ~ "Patient 4",
      sample == "GN01" ~ "Patient 5",
      sample == "NB01" ~ "Patient 6",
      sample == "NB02" ~ "Patient 6",
    ))


# Create a summary of sequencing depth per sample
depth_summary <- lapply(names(rds_list), function(f) {
  se <- rds_list[[f]]
  depths <- colData(se)$depth
  data.frame(
    file = f,
    total_depth = sum(depths, na.rm = TRUE),
    mean_depth = mean(depths, na.rm = TRUE),
    median_depth = median(depths, na.rm = TRUE)
  )
}) %>%
  bind_rows()

summary_df <- summary_df %>%
  left_join(depth_summary, by = "file") %>%
  mutate(
    n_rows_per_100k_depth = n_rows / (total_depth / 1e5),
    n_rows_per_mean_depth = n_rows / mean_depth
  )

summary_df

#Union calling per Patient:


# 1. Extract mapping: which file belongs to which patient
patient_map <- summary_df %>%
  dplyr::select(file, Patient)

# 2. For each patient, collect union of variant names across files
patient_variants <- lapply(unique(patient_map$Patient), function(pat) {
  # get files for this patient
  files <- patient_map %>% filter(Patient == pat) %>% pull(file)
  
  # get variant names from all files for this patient
  all_variants <- unlist(lapply(files, function(f) rownames(rds_list[[f]])))
  
  # return unique variants and patient label
  data.frame(
    Patient = pat,
    n_unique_variants = length(unique(all_variants))
  )
}) %>%
  bind_rows()


patient_variant_sets <- lapply(unique(patient_map$Patient), function(pat) {
  files <- patient_map %>% filter(Patient == pat) %>% pull(file)
  variants <- unlist(lapply(files, function(f) rownames(rds_list[[f]])))
  unique(variants)
})

names(patient_variant_sets) <- unique(patient_map$Patient)
patient_variants


# 1. Create patient-to-file mapping
patient_map <- summary_df %>% dplyr::select(file, Patient)

# 2. Calculate unique variant counts per patient
patient_union_counts <- lapply(unique(patient_map$Patient), function(pat) {
  files <- patient_map %>% filter(Patient == pat) %>% pull(file)
  all_variants <- unlist(lapply(files, function(f) rownames(rds_list[[f]])))
  data.frame(
    Patient = pat,
    n_unique_variants_union = length(unique(all_variants))
  )
}) %>%
  bind_rows()

# 3. Merge union count back into summary_df
summary_df <- summary_df %>%
  left_join(patient_union_counts, by = "Patient")


depth_per_patient <- summary_df %>%
  group_by(Patient) %>%
  summarise(total_depth_patient = sum(total_depth, na.rm = TRUE), .groups = "drop")

union_df <- summary_df %>%
  group_by(Patient, Age, project) %>%
  summarise(n_union = unique(n_unique_variants_union), .groups = "drop") %>%
  left_join(depth_per_patient, by = "Patient") %>%
  mutate(
    union_variants_per_100k_depth = n_union / (total_depth_patient / 1e5)
  )

# Merge union variant counts with depth per patient
union_df <- patient_union_counts %>%
  left_join(depth_per_patient, by = "Patient") %>%
  left_join(summary_df %>% dplyr::select(Patient, Age, project) %>% distinct(), by = "Patient") %>%
  mutate(
    union_variants_per_100k_depth = n_unique_variants_union / (total_depth_patient / 1e5)
  ) %>%
  arrange(project, Patient)


# View result
print(summary_df)

write.table(summary_df, "summary_variants_per_age.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


