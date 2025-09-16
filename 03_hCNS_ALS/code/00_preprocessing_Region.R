# packages 
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)


# Define the list of regions you want to process
regions <- c("Region_5","Region_10", "Region_28", "Region_29")  # Add more regions as needed

# Loop through each region
for (region in regions) {
  
  # Define directories for each region
  input_dir <- paste0("/../data/", region)
  output_dir <- paste0("/../test_output/", region)
  mgatk_dir <- file.path(input_dir, "mgatk/final")
  amulet_dir <- file.path(input_dir, "amulet/MultipletBarcodes_01.txt")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load counts and metadata from cellranger-atac
  counts <- Read10X_h5(filename = file.path(input_dir, "filtered_peak_bc_matrix.h5"))
  
  metadata <- read.csv(
    file = file.path(input_dir, "singlecell.csv"),
    header = TRUE,
    row.names = 1
  )
  
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = file.path(input_dir, "fragments.tsv.gz"),
    min.cells = 5,
    min.features = 50
  )
  
  df <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )
  
  # Extract gene annotations from EnsDb
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- 'UCSC'
  Annotation(df) <- annotations
  
  # Compute nucleosome signal and TSS enrichment
  df <- NucleosomeSignal(object = df)
  df <- TSSEnrichment(object = df, fast = FALSE)
  
  df$pct_reads_in_peaks <- df$peak_region_fragments / df$passed_filters * 100
  df$blacklist_ratio <- df$blacklist_region_fragments / df$peak_region_fragments
  df$high.tss <- ifelse(df$TSS.enrichment > 2, 'High', 'Low')
  
  # Remove low-quality cells
  df1 <- subset(df,
                subset = nCount_peaks > 1000 &
                  TSS.enrichment > 1.5 & 
                  nucleosome_signal < 1.5
  )
  
  # Load mgatk output
  mito.data <- ReadMGATK(dir = mgatk_dir)
  mito <- CreateAssayObject(counts = mito.data$counts)
  mito <- subset(mito, cells = colnames(df1))
  df1[["mito"]] <- mito
  df1 <- AddMetaData(df1, metadata = mito.data$depth, col.name = "mtDNA_depth")
  
  df1 <- subset(df1, mtDNA_depth >= 5)
  
  # Adding Amulet /Doublet Information
  multi <- read.table(amulet_dir, sep = "\t", header = FALSE)
  
  cell_ids <- as.data.frame(WhichCells(df1))
  
  colnames(multi) <- c("cell_barcodes")
  colnames(cell_ids) <- c("cell_barcodes")
  
  # Identify common characters
  common_chars <- intersect(multi$cell_barcodes, cell_ids$cell_barcodes)
  
  cell_ids$doublet <- ifelse(cell_ids$cell_barcodes %in% common_chars, "multiplet", "singlet")
  rownames(cell_ids) <- cell_ids$cell_barcodes
  cell_ids <- cell_ids[, -1]
  
  df1 <- AddMetaData(df1, metadata = cell_ids, col.name = "doublets")
  
  df1 <- subset(df1, subset = doublets == "singlet")
  
  # Save the preprocessed RDS file
  saveRDS(df1, file.path(output_dir, paste0("preprocessed_", region, ".rds")))
  
}
