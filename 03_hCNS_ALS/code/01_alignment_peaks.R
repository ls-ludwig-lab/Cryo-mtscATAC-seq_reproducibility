library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)



plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

# Define base path
dataset_base_path <- "/../data/"


Region_5 <- read.table(
  file = file.path(dataset_base_path,"Region_5/peaks.bed"),
  col.names = c("chr", "start", "end")
)

Region_10 <- read.table(
  file = file.path(dataset_base_path,"Region_10/peaks.bed"),
  col.names = c("chr", "start", "end")
)

Region_28 <- read.table(
  file = file.path(dataset_base_path,"Region_28/peaks.bed"),
  col.names = c("chr", "start", "end")
)
Region_29 <- read.table(
  file = file.path(dataset_base_path,"Region_29/peaks.bed"),
  col.names = c("chr", "start", "end")
)


gr.Region_5 <- makeGRangesFromDataFrame(Region_5)
gr.Region_10 <- makeGRangesFromDataFrame(Region_10)
gr.Region_28 <- makeGRangesFromDataFrame(Region_28)
gr.Region_29 <- makeGRangesFromDataFrame(Region_29)


# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.Region_5, gr.Region_10, gr.Region_28, gr.Region_29))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks


md.Region_5 <- read.table(
  file = file.path(dataset_base_path,"Region_5/singlecell.csv"),
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row


md.Region_10 <- read.table(
  file = file.path(dataset_base_path,"Region_10/singlecell.csv"),
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row


md.Region_28 <- read.table(
  file = file.path(dataset_base_path,"Region_28/singlecell.csv"),
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row


md.Region_29 <- read.table(
  file = file.path(dataset_base_path,"Region_29/singlecell.csv"),
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.Region_5 <- md.Region_5[md.Region_5$passed_filters > 500, ]
md.Region_10 <- md.Region_10[md.Region_10$passed_filters > 500, ]
md.Region_28 <- md.Region_28[md.Region_28$passed_filters > 500, ]
md.Region_29 <- md.Region_29[md.Region_29$passed_filters > 500, ]


# create fragment objects

frags.Region_5 <- CreateFragmentObject(
  path = file.path(dataset_base_path,"Region_5/fragments.tsv.gz"),
  cells = rownames(md.Region_5)
)
frags.Region_10 <- CreateFragmentObject(
  path = file.path(dataset_base_path,"Region_10/fragments.tsv.gz"),
  cells = rownames(md.Region_10)
)

frags.Region_28 <- CreateFragmentObject(
  path = file.path(dataset_base_path,"Region_28/fragments.tsv.gz"),
  cells = rownames(md.Region_28 )
)
frags.Region_29 <- CreateFragmentObject(
  path = file.path(dataset_base_path,"Region_29/fragments.tsv.gz"),
  cells = rownames(md.Region_29)
)


Region_5.counts <- FeatureMatrix(
  fragments = frags.Region_5,
  features = combined.peaks,
  cells = rownames(md.Region_5)
)

Region_10.counts <- FeatureMatrix(
  fragments = frags.Region_10,
  features = combined.peaks,
  cells = rownames(md.Region_10)
)

Region_28.counts <- FeatureMatrix(
  fragments = frags.Region_28,
  features = combined.peaks,
  cells = rownames(md.Region_28)
)
Region_29.counts <- FeatureMatrix(
  fragments = frags.Region_29,
  features = combined.peaks,
  cells = rownames(md.Region_29)
)


Region_5_assay <- CreateChromatinAssay(Region_5.counts, fragments = frags.Region_5)
Region_5 <- CreateSeuratObject(Region_5_assay, assay = "ATAC", meta.data=md.Region_5)
Region_5

Region_10_assay <- CreateChromatinAssay(Region_10.counts, fragments = frags.Region_10)
Region_10 <- CreateSeuratObject(Region_10_assay, assay = "ATAC", meta.data=md.Region_10)
Region_10

Region_28_assay <- CreateChromatinAssay(Region_28.counts, fragments = frags.Region_28)
Region_28 <- CreateSeuratObject(Region_28_assay, assay = "ATAC", meta.data=md.Region_28)
Region_28

Region_29_assay <- CreateChromatinAssay(Region_29.counts, fragments = frags.Region_29)
Region_29 <- CreateSeuratObject(Region_29_assay, assay = "ATAC", meta.data=md.Region_29)
Region_29

#add information to identify dataset of origin

# Define base path
meta_base_path <- "/../test_output/"


#Region_28
Region_5$dataset <- 'Region_5'
Region_5 <- AddMetaData(object = Region_5, metadata = "Region", col.name = "Donor")
Region_5 <- AddMetaData(object = Region_5, metadata = "Region_5", col.name = "Region")

#meta.data from single processed data
df <- readRDS(file=file.path(meta_base_path,"Region_5/preprocessed.rds"))
meta.data <- df@meta.data
cell_ids <- rownames(meta.data)
Region_5 <- subset(Region_5, cells = cell_ids)
# Filter the data frame    
Region_5 <- AddMetaData(object = Region_5, metadata = meta.data$mtDNA_depth, col.name = "mtDNA_depth")
Region_5 <- AddMetaData(object = Region_5, metadata = meta.data$doublets, col.name = "doublets")

#Region_28
Region_10$dataset <- 'Region_10'
Region_10 <- AddMetaData(object = Region_10, metadata = "Region", col.name = "Donor")
Region_10 <- AddMetaData(object = Region_10, metadata = "Region_10", col.name = "Region")
#meta.data from single processed data
df <- readRDS(file=file.path(meta_base_path,"Region_10/preprocessed.rds"))
meta.data <- df@meta.data
cell_ids <- rownames(meta.data)
Region_10 <- subset(Region_10, cells = cell_ids)
# Filter the data frame    
Region_10 <- AddMetaData(object = Region_10, metadata = meta.data$mtDNA_depth, col.name = "mtDNA_depth")
Region_10 <- AddMetaData(object = Region_10, metadata = meta.data$doublets, col.name = "doublets")


#BB38_28
Region_28$dataset <- 'Region_28'
Region_28 <- AddMetaData(object = Region_28, metadata = "Region", col.name = "Donor")
Region_28 <- AddMetaData(object = Region_28, metadata = "Region_28", col.name = "Region")
#meta.data from single processed data
df <- readRDS(file=file.path(meta_base_path,"Region_28/preprocessed.rds"))
meta.data <- df@meta.data
cell_ids <- rownames(meta.data)
Region_28 <- subset(Region_28, cells = cell_ids)
# Filter the data frame    
Region_28 <- AddMetaData(object = Region_28, metadata = meta.data$mtDNA_depth, col.name = "mtDNA_depth")
Region_28 <- AddMetaData(object = Region_28, metadata = meta.data$doublets, col.name = "doublets")

#Region_29
Region_29$dataset <- 'Region_29'
Region_29 <- AddMetaData(object = Region_29, metadata = "BB29", col.name = "Donor")
Region_29 <- AddMetaData(object = Region_29, metadata = "Region_29", col.name = "Region")
#meta.data from single processed data
df <- readRDS(file=file.path(meta_base_path,"Region_29/preprocessed.rds"))
meta.data <- df@meta.data
cell_ids <- rownames(meta.data)
Region_29 <- subset(Region_29, cells = cell_ids)
# Filter the data frame    
Region_29 <- AddMetaData(object = Region_29, metadata = meta.data$mtDNA_depth, col.name = "mtDNA_depth")
Region_29 <- AddMetaData(object = Region_29, metadata = meta.data$doublets, col.name = "doublets")

Region_29

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = Region_5,
  y = list(Region_10, Region_28, Region_29),
  add.cell.ids = c("Region_5", "Region_10", "Region_28", "Region_29")
)
combined[["ATAC"]]

setwd("/../test_output")

granges(combined)
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

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

combined$nucleosome_group <- ifelse(combined$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragHisto <- FragmentHistogram(object = combined, group.by = 'nucleosome_group')

# Augment QC metrics that were computed by cellranger-atac
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
combined$pct_reads_in_DNase <- combined$DNase_sensitive_region_fragments / combined$passed_filters * 100
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments


#the filtering was already processed in the single experiment as I added the meta data - so no filtering here for the combined data set
#makes also sense as the data quality differs between samples 


combined

DefaultAssay(combined) <- "ATAC"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 10)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = "lsi", dims = 2:50)
combined <- FindNeighbors(combined, reduction = "lsi", dims = 2:50)
combined <- FindClusters(combined, resolution = 0.6, algorithm = 3)

saveRDS(combined, "merged.rds")


