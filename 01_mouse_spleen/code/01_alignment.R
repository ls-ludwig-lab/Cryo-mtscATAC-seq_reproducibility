library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Mmusculus.v79)

plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM


# read in peak sets
FF_10min <- read.table(
  file = ".../FF_10min/data/CR-ATAC_mSpleen_FF_10min_mtMask/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)


FF_30min <- read.table(
  file = ".../FF_30min/data/CR-ATAC_mSpleen_FF_30min_mtMask/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)


Cells <- read.table(
  file = ".../fresh_isolated/data/peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
gr.FF_10min <- makeGRangesFromDataFrame(FF_10min)


gr.FF_30min <- makeGRangesFromDataFrame(FF_30min)


gr.Cells <- makeGRangesFromDataFrame(Cells)


# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.FF_10min, gr.FF_30min, gr.Cells))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
#head(peakwidths)
#[1]  886  961 1090  989  917  906
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks


# Create Fragment objects
#quantify our combined set of peaks we’ll need to create a Fragment object for each experiment. The Fragment class is a specialized class defined in Signac to hold all the information related to a single fragment file.

# load metadata
md.FF_10min <- read.table(
  file = ".../FF_10min/data/CR-ATAC_mSpleen_FF_10min_mtMask/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.FF_30min <- read.table(
  file = ".../FF_30min/data/CR-ATAC_mSpleen_FF_30min_mtMask/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]


md.Cells <- read.table(
  file = ".../fresh_isolated/data/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]


# perform an initial filtering of low count cells
md.FF_10min <- md.FF_10min[md.FF_10min$passed_filters > 200, ]
md.FF_30min <- md.FF_30min[md.FF_30min$passed_filters > 200, ]
md.Cells <- md.Cells[md.Cells$passed_filters > 200, ]

# create fragment objects
frags.FF_10min <- CreateFragmentObject(
  path = ".../FF_10min/data/CR-ATAC_mSpleen_FF_10min_mtMask/outs/fragments.tsv.gz",
  cells = rownames(md.FF_10min)
)

frags.FF_30min<- CreateFragmentObject(
  path = ".../FF_30min/data/CR-ATAC_mSpleen_FF_30min_mtMask/outs/fragments.tsv.gz",
  cells = rownames(md.FF_30min)
)

frags.Cells <- CreateFragmentObject(
  path = ".../fresh_isolated/data/fragments.tsv.gz",
  cells = rownames(md.Cells)
)


# Quantify peaks in each dataset
FF_10min.counts <- FeatureMatrix(
  fragments = frags.FF_10min,
  features = combined.peaks,
  cells = rownames(md.FF_10min)
)

FF_30min.counts <- FeatureMatrix(
  fragments = frags.FF_30min,
  features = combined.peaks,
  cells = rownames(md.FF_30min)
)

Cells.counts <- FeatureMatrix(
  fragments = frags.Cells,
  features = combined.peaks,
  cells = rownames(md.Cells)
)

# Create the objects

FF_10min_assay <- CreateChromatinAssay(FF_10min.counts, fragments = frags.FF_10min)
FF_10min <- CreateSeuratObject(FF_10min_assay, assay = "ATAC", meta.data=md.FF_10min)

FF_10min


FF_30min_assay <- CreateChromatinAssay(FF_30min.counts, fragments = frags.FF_30min)
FF_30min <- CreateSeuratObject(FF_30min_assay, assay = "ATAC", meta.data=md.FF_30min)
FF_30min


Cells_assay <- CreateChromatinAssay(Cells.counts, fragments = frags.Cells)
Cells <- CreateSeuratObject(Cells_assay, assay = "ATAC", meta.data=md.Cells)

Cells


FF_10min$dataset <- 'FF_10min'
FF_10min <- AddMetaData(object = FF_10min, metadata = "fresh frozen", col.name = "Specimen")
FF_10min <- AddMetaData(object = FF_10min, metadata = "10 min", col.name = "Fixation")



#meta.data from single processed data
df <- readRDS(file=".../FF_10min/output/preprocessed.rds")
meta.data <- df@meta.data

cell_ids <- rownames(meta.data)
FF_10min <- subset(FF_10min, cells = cell_ids)
FF_10min <- AddMetaData(object = FF_10min, metadata = meta.data$mtDNA_depth, col.name = "mtDNA_depth")
FF_10min <- AddMetaData(object = FF_10min, metadata = meta.data$doublets, col.name = "doublets")
FF_10min

#FF30_min
FF_30min$dataset <- 'FF30_min'

df <- readRDS(file=".../FF_30min/output/preprocessed.rds")
meta.data <- df@meta.data
cell_ids <- rownames(meta.data)
FF_30min <- subset(FF_30min, cells = cell_ids)

FF_30min <- AddMetaData(object = FF_30min, metadata = meta.data$mtDNA_depth, col.name = "mtDNA_depth")
FF_30min <- AddMetaData(object = FF_30min, metadata = meta.data$doublets, col.name = "doublets")
FF_30min

#Cells
Cells$dataset <- 'Cells'

df <- readRDS(file=".../fresh_isolated/output/preprocessed.rds")
meta.data <- df@meta.data

cell_ids <- rownames(meta.data)
Cells <- subset(Cells, cells = cell_ids)

# Check if row names match
all(rownames(meta.data) %in% colnames(Cells))

# Subset meta.data to only include cells present in Cells
meta.filtered <- meta.data[colnames(Cells), , drop = FALSE]

Cells <- AddMetaData(object = Cells, metadata = meta.filtered$mtDNA_depth, col.name = "mtDNA_depth")
Cells <- AddMetaData(object = Cells, metadata = meta.filtered$doublets, col.name = "doublets")

Cells


# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = FF_10min,
  y = list(FF_30min, Cells),
  add.cell.ids = c("FF_10min", "FF_30min", "Cells")
)

combined

granges(combined)
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(combined) <- annotations

# compute nucleosome signal score per cell
combined <- NucleosomeSignal(object = combined)

# compute TSS enrichment score per cell
combined <- TSSEnrichment(object = combined, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments

combined$high.tss <- ifelse(combined$TSS.enrichment > 2, 'High', 'Low')
TSS <- TSSPlot(combined, group.by = 'high.tss') + NoLegend()
TSS

combined$nucleosome_group <- ifelse(combined$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragHisto <- FragmentHistogram(object = combined, group.by = 'nucleosome_group')
FragHisto


# Augment QC metrics that were computed by cellranger-atac
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
combined$pct_reads_in_DNase <- combined$DNase_sensitive_region_fragments / combined$passed_filters * 100
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments

Vln <- VlnPlot(combined, c("TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks", "peak_region_fragments"), pt.size = 0, ncol = 2, cols = "#72ACA7") 
Vln

DensSca <- DensityScatter(combined, x = 'peak_region_fragments', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
DensSca

combined


Vln <- VlnPlot(combined, c("TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks", "peak_region_fragments", "mtDNA_depth"), pt.size = 0, ncol = 2, cols = "#72ACA7") 
Vln
ggsave(plot = Vln, filename = ".../integration/output/Vln_afterQC.png", width = 20, height = 20, dpi=300)


DefaultAssay(combined) <- "ATAC"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 10)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = "lsi", dims = 2:50)
combined <- FindNeighbors(combined, reduction = "lsi", dims = 2:50)
combined <- FindClusters(combined, resolution = 0.6, algorithm = 3)


Dim <- DimPlot(combined, group.by="dataset", raster = FALSE) 
Dim

ggsave(plot = Dim, filename = ".../integration/output/Dim_afterQC_1.png", width = 7.75, height = 4.55, dpi=300)
saveRDS(combined, ".../integration/output/combined_1.rds")
