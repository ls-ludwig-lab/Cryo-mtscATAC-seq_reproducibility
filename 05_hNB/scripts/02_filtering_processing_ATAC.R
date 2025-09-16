# load libraries
library(Signac)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(devtools)
library(SeuratDisk)
library(SeuratData)
library(patchwork)
library(SummarizedExperiment)
library(irlba)
library(cowplot)
library(scales)
library(hdf5r)
library(GenomeInfoDb)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(RSpectra)
library(SeuratObject)
library(BiocGenerics)
library(biovizBase)
library(foreach)
library(doParallel)
library(readxl)
library(Matrix)

set.seed(257)

# =========================
# Load Seurat Objects
# =========================
data_dir <- paste0("...")        # Placeholder: folder containing Seurat objects

# List all Seurat RDS files starting with 'Seurat_'
seurat_files <- list.files(data_dir, pattern = "^Seurat_", full.names = TRUE)

# Loop over files and load each object into environment
for (file_path in seurat_files) {
  file_name <- basename(file_path)
  
  # Extract identifier: remove 'Seurat_' prefix
  identifier <- sub("^Seurat_", "", file_name)
  
  # Load Seurat object
  seurat_obj <- readRDS(file_path)
  
  # Assign object to environment with identifier name
  assign(identifier, seurat_obj, envir = .GlobalEnv)
}


# =========================
# Define filtering function
# =========================
Filter_ATAC <- function(seurat, minPeakFrag, maxPeakFrag, minPercPeak, maxNucleoSig, minTSS) {
  # Subset Seurat object based on QC metrics
  seurat <- subset(
    seurat,
    subset = pct_reads_in_peaks >= minPercPeak &
      peak_region_fragments >= minPeakFrag &
      peak_region_fragments < maxPeakFrag &
      TSS.enrichment >= minTSS &
      nucleosome_signal < maxNucleoSig
  )
  return(seurat)
}


# =========================
# Load threshold parameters
# =========================
thresholds <- read.table(
  "Filtering_Thresholds.tsv", 
  sep = "\t", 
  header = TRUE, 
  stringsAsFactors = FALSE
)


# =========================
# Identify all Seurat objects to filter
# =========================
seurat_names <- "..."            # Placeholder: list of Seurat object names in environment
output_dir <- paste0("...")      # Placeholder: directory to save filtered objects


# =========================
# Apply filtering sequentially
# =========================
for (seurat_name in seurat_names) {
  
  # Retrieve Seurat object
  seurat_obj <- get(seurat_name)
  
  # Retrieve sample-specific threshold parameters
  params <- thresholds[thresholds$Sample_ID == seurat_name, ]
  
  if (nrow(params) > 0) {
    
    # Apply QC filtering
    filtered_obj <- Filter_ATAC(
      seurat = seurat_obj,
      minPeakFrag = params$minPeakFrag,
      maxPeakFrag = params$maxPeakFrag,
      minPercPeak = params$minPercPeak,
      maxNucleoSig = params$maxNucleoSig,
      minTSS = params$minTSS
    )
    
    # Save filtered Seurat object
    saveRDS(
      filtered_obj, 
      file = file.path(output_dir, paste0(seurat_name, "_filtered.rds"))
    )
    
    # Update object in environment
    assign(seurat_name, filtered_obj, envir = .GlobalEnv)
    
    cat("Filtered and saved:", seurat_name, "\n")
    
  } else {
    cat("No matching thresholds found for:", seurat_name, "\n")
  }
}


# =========================
# Normalization and linear dimensional reduction
# =========================
dir.output <- "path_to_output_directory"   # Placeholder: directory for output
selected_resolutions <- c(0.5)

filtered_object_names <- "..."  # Placeholder: list of filtered Seurat object names

for (name in filtered_object_names) {
  
  # Retrieve filtered Seurat object
  seurat <- get(name)
  
  # ATAC preprocessing
  seurat <- RunTFIDF(seurat)
  seurat <- FindTopFeatures(seurat, min.cutoff = 'q0')
  seurat <- RunSVD(seurat)
  
  # Set default assay to peaks for downstream analysis
  DefaultAssay(seurat) <- "peaks"
  
  # Dimensionality reduction and clustering
  seurat <- RunUMAP(seurat, reduction = 'lsi', dims = 2:30, reduction.name = "UMAP")
  seurat <- FindNeighbors(seurat, reduction = 'lsi', dims = 2:30)
  seurat <- FindClusters(seurat, resolution = selected_resolutions, verbose = FALSE, algorithm = 3)
  
  # Update object in environment
  assign(name, seurat)
}


# =========================
# Create Gene Activity Assay
# =========================
object_names <- "..."  # Placeholder: list of Seurat objects to process

for (seurat_name in object_names) {
  
  # Retrieve Seurat object
  seurat <- get(seurat_name)
  
  # Compute gene activity matrix
  gene.activities <- GeneActivity(seurat)
  
  # Add gene activity matrix as new RNA assay
  seurat[['RNA']] <- CreateAssayObject(counts = gene.activities)
  
  # Normalize RNA assay
  seurat <- NormalizeData(
    seurat,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(seurat$nCount_RNA)
  )
  
  # Save updated Seurat object
  assign(seurat_name, seurat)
  saveRDS(seurat, file = file.path(output_dir, paste0(seurat_name, "_filtered.rds")))
}