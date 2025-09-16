# Load libraries
library(Signac)
library(Seurat)
library(SummarizedExperiment)
library(irlba)
library(GenomeInfoDb)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(RSpectra)
library(Rsamtools)
library(SeuratObject)
library(BiocGenerics)
library(biovizBase)


# ============================================
# Functions
# ============================================

# Function to find common peaks across samples
FindCommonPeaks_scATACseq <- function(dir.data.samples, MaxPeakWidth, MinPeakWidth, dir.output){
  
  # Read peaks for each sample
  peak.list <- lapply(dir.data.samples, function(sample_dir) {
    read.table(file.path(sample_dir, "filtered_peak_bc_matrix/peaks.bed"),
               col.names = c("chr", "start", "end"))
  })
  
  # Convert to GRanges
  gr.list <- lapply(peak.list, makeGRangesFromDataFrame)
  
  # Combine all peaks
  combined.peaks <- Reduce(c, gr.list)
  
  # Merge overlapping peaks
  combined.peaks <- Signac::reduce(combined.peaks)
  
  # Filter by peak width
  peak.width <- width(combined.peaks)
  combined.peaks <- combined.peaks[peak.width < MaxPeakWidth & peak.width > MinPeakWidth]
  
  # Save BED file
  write.table(as.data.frame(combined.peaks),
              file.path(dir.output, "CommonSetOfPeaks.bed"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Function to create a Seurat object for one ATAC sample
CreateATACobject <- function(sample.ID, dir.data.sample, MinCounts, dir.output, gannotation, name.ID, integration=TRUE){
  
  # Select peak set
  if(integration){
    peaks <- read.table(file.path(dir.output, "CommonSetOfPeaks.bed"))
  } else {
    peaks <- read.table(file.path(dir.data.sample, "filtered_peak_bc_matrix/peaks.bed"))
  }
  
  peaks <- makeGRangesFromDataFrame(peaks)
  
  # Load metadata
  meta <- read.table(file.path(dir.data.sample, "singlecell.csv"),
                     sep = ",", header = TRUE, stringsAsFactors = FALSE, row.names = 1)[-1,]
  
  # Filter cells by minimum counts
  meta <- meta[meta$passed_filters > MinCounts, ]
  
  # Keep only barcodes present in the matrix
  cells.to.keep <- read.table(file.path(dir.data.sample, "filtered_peak_bc_matrix/barcodes.tsv"))$V1
  cells.to.keep <- intersect(cells.to.keep, rownames(meta))
  meta <- meta[cells.to.keep, ]
  
  # Create fragment object
  frag <- CreateFragmentObject(path = file.path(dir.data.sample, "fragments.tsv.gz"),
                               cells = rownames(meta))
  
  # Build feature matrix
  counts <- FeatureMatrix(fragments = frag, features = peaks, cells = rownames(meta))
  
  # Create Seurat object with ChromatinAssay
  chrom_assay <- CreateChromatinAssay(counts = counts, fragments = frag)
  seurat_obj <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = meta)
  
  # Add sample ID
  seurat_obj <- AddMetaData(seurat_obj, metadata = sample.ID, col.name = name.ID)
  
  # Remove doublets if AMULET output exists
  doublets_file <- file.path(dir.data.sample, "amulet/MultipletBarcodes_01.txt")
  if(file.exists(doublets_file)){
    doublets <- read.table(doublets_file)$V1
    seurat_obj <- subset(seurat_obj, cells = setdiff(Cells(seurat_obj), doublets))
  }
  
  # Add gene annotation
  Annotation(seurat_obj) <- gannotation
  
  # Calculate quality metrics
  seurat_obj <- NucleosomeSignal(seurat_obj)
  seurat_obj <- TSSEnrichment(seurat_obj, fast = FALSE)
  seurat_obj <- AddMetaData(seurat_obj,
                            metadata = seurat_obj[["peak_region_fragments"]] / seurat_obj[["passed_filters"]] * 100,
                            col.name = "pct_reads_in_peaks")
  
  # Save Seurat object
  saveRDS(seurat_obj, file.path(dir.output, "Objects", paste0("Seurat_", sample.ID, ".rds")))
  
  return(seurat_obj)
}

# ============================================
# Parameters
# ============================================

samples <- read.table("path_to_samples.tsv", header = TRUE, stringsAsFactors = FALSE) # Columns: Alias, Directory
dir.output <- "path_to_output_dir"
integration <- TRUE
MinCounts <- 5
MaxPeakWidth <- 10000
MinPeakWidth <- 20
name.ID <- "Sample"

# ============================================
# Run common peak finding
# ============================================

if(integration){
  FindCommonPeaks_scATACseq(dir.data.samples = samples$Directory,
                            MaxPeakWidth = MaxPeakWidth,
                            MinPeakWidth = MinPeakWidth,
                            dir.output = dir.output)
}

# ============================================
# Gene annotation
# ============================================

gene_annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
gene_annotation <- renameSeqlevels(gene_annotation, mapSeqlevels(seqlevels(gene_annotation), "UCSC"))
genome(gene_annotation) <- "hg38"

# ============================================
# Create Seurat objects for each sample
# ============================================

seurat_objects <- list()

for(i in 1:nrow(samples)){
  seurat_obj <- CreateATACobject(sample.ID = samples$Alias[i],
                                 dir.data.sample = samples$Directory[i],
                                 MinCounts = MinCounts,
                                 dir.output = dir.output,
                                 gannotation = gene_annotation,
                                 name.ID = name.ID,
                                 integration = integration)
  seurat_objects[[samples$Alias[i]]] <- seurat_obj
}