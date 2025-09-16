# Required libraries
library(Signac)
library(Seurat)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(data.table)

# Define dataset info
datasets <- list(
    name = "0p5pct",
    base_path = ".../.../"
)

# Loop over datasets
for (ds in datasets) {
  
  name <- ds$name
  data_path <- file.path(ds$base_path, "data")
  output_path <- file.path(ds$base_path, "output")
  
  # Set working directory
  setwd(data_path)
  
  # Read data
  counts <- Read10X_h5("filtered_peak_bc_matrix.h5")
  metadata <- read.csv("singlecell.csv", header = TRUE, row.names = 1)
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = "fragments.tsv.gz",
    min.cells = 5,
    min.features = 50
  )
  df <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = metadata)
  
  # Annotation
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- 'UCSC'
  Annotation(df) <- annotations
  
  # QC metrics
  df <- NucleosomeSignal(df)
  df <- TSSEnrichment(df, fast = FALSE)
  df$pct_reads_in_peaks <- df$peak_region_fragments / df$passed_filters * 100
  df$blacklist_ratio <- df$blacklist_region_fragments / df$peak_region_fragments
  df$high.tss <- ifelse(df$TSS.enrichment > 2, 'High', 'Low')
  df$nucleosome_group <- ifelse(df$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  
  # Plots
  setwd(output_path)
  FragHisto <- FragmentHistogram(df, group.by = 'nucleosome_group')
  ggsave(FragHisto, filename = paste0("FragHisto_", name, ".png"), width = 5, height = 5, dpi=300)
  
  Vln <- VlnPlot(df, c("TSS.enrichment", "nCount_peaks", "nucleosome_signal", "pct_reads_in_peaks"), pt.size = 0, ncol = 2, cols = "#72ACA7")
  ggsave(Vln, filename = paste0("Vln_", name, ".png"), width = 8, height = 8, dpi=300)
  
  DensSca <- DensityScatter(df, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
  ggsave(DensSca, filename = paste0("DensSca_", name, ".png"), width = 20, height = 20, dpi=300)
  
  # Doublets via AMULET
  amulet_path <- file.path(data_path, "amulet", "MultipletBarcodes_01.txt")
  if (file.exists(amulet_path)) {
    multi <- read.table(amulet_path, sep = "\t", header = FALSE)
    cell_ids <- data.frame(cell_barcodes = colnames(df))
    colnames(multi) <- "cell_barcodes"
    cell_ids$doublet <- ifelse(cell_ids$cell_barcodes %in% multi$cell_barcodes, "multiplet", "singlet")
    rownames(cell_ids) <- cell_ids$cell_barcodes
    df <- AddMetaData(df, metadata = cell_ids["doublet"])
  }
  
  # Add mito assay
  mito_dir <- file.path(data_path, "mgatk", "final")
  if (dir.exists(mito_dir)) {
    mito.data <- ReadMGATK(mito_dir)
    mito <- CreateAssayObject(counts = mito.data$counts)
    mito <- subset(mito, cells = colnames(df))
    df[["mito"]] <- mito
    df <- AddMetaData(df, metadata = mito.data$depth, col.name = "mtDNA_depth")
    
    features <- c("TSS.enrichment", "nCount_peaks", "nucleosome_signal", "pct_reads_in_peaks", "mtDNA_depth")
    p <- VlnPlot(df, features, pt.size = 0, ncol = 3, cols = "#72ACA7")
    ggsave(p, filename = paste0("Vln_afterQC_mtDNA_", name, ".png"), width = 20, height = 20, dpi=300)
  }
  
  # Add Vireo donor info
  vireo_path <- file.path(data_path, "vireo", "donor_ids.tsv")
  if (file.exists(vireo_path)) {
    vireo <- read.csv(vireo_path, sep='\t', header=TRUE, row.names=1)
    vireo_filtered <- vireo[rownames(vireo) %in% colnames(df), ]
    df <- AddMetaData(df, metadata = vireo_filtered$donor_id, col.name = "vireo_donor")
    df <- AddMetaData(df, metadata = vireo_filtered$doublet_logLikRatio, col.name = "vireo_doublet_logLikRatio")
  }
  
  # Save object
  saveRDS(df, file = file.path("/.../rds_files", paste0("unfiltered_", name, ".rds")))
}
