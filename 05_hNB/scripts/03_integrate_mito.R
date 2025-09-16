# Load required libraries
library(Signac)                
library(Seurat)                
library(ggplot2)               
library(tidyverse)             
library(dplyr)                 
library(devtools)              
library(SummarizedExperiment)  

set.seed(257)                  


# Load Seurat Objects
data_dir <- paste0("...")      # Placeholder: folder containing Seurat RDS files

# List all RDS files in the directory
seurat_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)

# Loop through files and load each object into environment
for (file_path in seurat_files) {
  file_name <- basename(file_path)
  # Extract object identifier from filename (placeholder)
  identifier <- sub("...", "...", file_name)  
  seurat_obj <- readRDS(file_path)
  
  # Assign object to environment
  assign(identifier, seurat_obj, envir = .GlobalEnv)
}


# =========================
# Define functions
# =========================

# --- Function to add mitochondrial genotyping (MGATK) information to a Seurat object ---
Add_MGATK <- function(seurat, dir.data.sample, sample.ID, MinCellVar, MinStrandCor, MinVMR, minMitoDepth) {
  # Load MGATK output for the sample
  mgatk <- ReadMGATK(paste0(dir.data.sample, "/mgatk/"))  # Read MGATK counts and reference
  
  # Create Seurat assay from MGATK counts
  mgatk.assay <- CreateAssayObject(counts = mgatk$counts)
  
  # Subset assay to include only cells present in Seurat object
  mgatk.assay <- subset(mgatk.assay, cells = Cells(seurat))
  
  # Add the assay to Seurat object under name "mito"
  seurat[["mito"]] <- mgatk.assay
  
  # Add mtDNA depth as metadata
  seurat <- AddMetaData(seurat, mgatk$depth, col.name = "mtDNA_depth")
  
  # Identify variable mitochondrial variants
  variable.sites <- IdentifyVariants(seurat, assay = "mito", refallele = mgatk$refallele)
  
  # Save all identified variants
  write.table(variable.sites, paste0(dir.output, "/AllVariants_", sample.ID, ".txt"))
  
  # Filter high-confidence variants based on thresholds
  high.conf <- subset(variable.sites, 
                      n_cells_conf_detected >= MinCellVar &
                        strand_correlation >= MinStrandCor &
                        vmr > MinVMR)
  
  # Save high-confidence variants
  write.table(high.conf, paste0(dir.output, "/HighConfVariants_", sample.ID, ".txt"))
  
  return(seurat)
}


# --- Function to calculate heteroplasmy for a sample based on variants ---
Calculate_Heteroplasmy <- function(sample.ID, variants, dir.output){
  # Load previously saved Seurat object for the sample
  seurat <- readRDS(paste0(dir.output, "/SeuratMito_", sample.ID, ".rds"))
  
  # Compute allele frequencies for specified variants
  seurat <- AlleleFreq(object = seurat, variants = variants, assay = "mito")
  
  # Assign updated Seurat object back to environment (optional)
  assign(sample.ID, seurat, envir = .GlobalEnv)
  
  # Save updated Seurat object
  saveRDS(seurat, paste0(dir.output, "/SeuratMito_", sample.ID, ".rds")) 
}


# =========================
# Sequential processing of all Seurat objects
# =========================

# Placeholder vector: names of Seurat objects in the environment
seurat_names <- c("...")  

# Load filtering thresholds for each sample
thresholds <- read.table("Filtering_Thresholds.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Placeholder: directories containing MGATK output per sample
directories <- paste0("...")  

# Output directory for processed objects and variant tables
dir.output <- "path_to_output_dir"  # Replace with actual output path

# Loop over all Seurat objects and add mitochondrial genotyping
for (name in seurat_names) {
  sample.ID <- name  # Name of Seurat object is used as sample ID
  
  # Retrieve thresholds for this sample
  sample.row <- thresholds %>% dplyr::filter(Sample_ID == sample.ID)
  
  # Retrieve directory containing MGATK output
  dir.row <- directories %>% dplyr::filter(Sample_ID == sample.ID)
  dir.data.sample <- dir.row$Directory
  
  # Get Seurat object from environment
  seurat.obj <- get(name)
  
  # Add MGATK information
  seurat.obj <- Add_MGATK(
    seurat = seurat.obj,
    dir.data.sample = dir.data.sample,
    sample.ID = sample.ID,
    MinCellVar = sample.row$MinCellVar,
    MinStrandCor = sample.row$MinStrandCor,
    MinVMR = sample.row$MinVMR,
    minMitoDepth = sample.row$minMitoDepth
  )
  
  # Assign updated object back to environment
  assign(name, seurat.obj, envir = .GlobalEnv)
  
  # Save processed Seurat object
  saveRDS(seurat.obj, paste0(dir.output, "/SeuratMito_", sample.ID, ".rds"))
  
  cat("Processed and saved Seurat object for sample:", sample.ID, "\n")
}


# =========================
# Union calling by patient
# =========================
## NB01 and NB02 

# Identify unique patients from thresholds table
patients <- unique(thresholds$Patient)

# Sequential processing per patient
for (patient in patients) {
  
  # Get all samples belonging to this patient
  patient_samples <- thresholds %>% dplyr::filter(Patient == patient) %>% pull(Sample_ID)
  
  # Collect variant files for the patient
  file.variants <- list.files(paste0(dir.output, "..."), full.names = TRUE, pattern = "HighConf")
  
  # Keep only files corresponding to the patient's samples
  file.variants <- grep(paste(patient_samples, collapse = "|"), file.variants, value = TRUE)
  
  # Extract unique variants across all samples
  variants <- unique(unlist(sapply(file.variants, function(x) read.table(x, header = TRUE)$variant)))
  
  # Compute heteroplasmy for each sample in the patient group
  for (sample.ID in patient_samples) {
    Calculate_Heteroplasmy(sample.ID = sample.ID, variants = variants, dir.output = dir.output)
  }
}

