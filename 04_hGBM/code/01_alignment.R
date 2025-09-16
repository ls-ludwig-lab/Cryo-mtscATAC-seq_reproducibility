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


# read in peak sets
patient1_primary <- read.table(
  file = "/../data/peaks.bed",
  col.names = c("chr", "start", "end")
)

patient2_recu <- read.table(
  file = "/../peaks.bed",
  col.names = c("chr", "start", "end")
)

patient1_recu <- read.table(
  file = "/../peaks.bed",
  col.names = c("chr", "start", "end")
)

patient2_primary <- read.table(
  file = "/../peaks.bed",
  col.names = c("chr", "start", "end")
)


# convert to genomic ranges
gr.patient1_primary <- makeGRangesFromDataFrame(patient1_primary)
gr.patient2_recu <- makeGRangesFromDataFrame(patient2_recu)
gr.patient1_recu <- makeGRangesFromDataFrame(patient1_recu)
gr.patient2_primary <- makeGRangesFromDataFrame(patient2_primary)


# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.patient1_primary, gr.patient2_recu, gr.patient1_recu, gr.patient2_primary))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# Create Fragment objects
#quantify our combined set of peaks we’ll need to create a Fragment object for each experiment. The Fragment class is a specialized class defined in Signac to hold all the information related to a single fragment file.

# load metadata
md.patient1_primary <- read.table(
  file = "/../singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row


md.patient2_recu <- read.table(
  file = "/../data/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.patient1_recu <- read.table(
  file = "/../singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.patient2_primary <- read.table(
  file = "/../singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row


# perform an initial filtering of low count cells
md.patient1_primary <- md.patient1_primary[md.patient1_primary$passed_filters > 500, ]
md.patient2_recu <- md.patient2_recu[md.patient2_recu$passed_filters > 500, ]
md.patient1_recu <- md.patient1_recu[md.patient1_recu$passed_filters > 500, ]
md.patient2_primary <- md.patient2_primary[md.patient2_primary$passed_filters > 500, ]

# create fragment objects
frags.patient1_primary <- CreateFragmentObject(
  path = "/../fragments.tsv.gz",
  cells = rownames(md.patient1_primary)
)


frags.patient2_recu <- CreateFragmentObject(
  path = "/../data/fragments.tsv.gz",
  cells = rownames(md.patient2_recu)
)

frags.patient1_recu<- CreateFragmentObject(
  path = "/../data/fragments.tsv.gz",
  cells = rownames(md.patient1_recu)
)

frags.patient2_primary <- CreateFragmentObject(
  path = "/../data/fragments.tsv.gz",
  cells = rownames(md.patient2_primary)
)


# Quantify peaks in each dataset
patient1_primary.counts <- FeatureMatrix(
  fragments = frags.patient1_primary,
  features = combined.peaks,
  cells = rownames(md.patient1_primary)
)


patient2_recu.counts <- FeatureMatrix(
  fragments = frags.patient2_recu,
  features = combined.peaks,
  cells = rownames(md.patient2_recu)
)


patient1_recu.counts <- FeatureMatrix(
  fragments = frags.patient1_recu,
  features = combined.peaks,
  cells = rownames(md.patient1_recu)
)


patient2_primary.counts <- FeatureMatrix(
  fragments = frags.patient2_primary,
  features = combined.peaks,
  cells = rownames(md.patient2_primary)
)

# Create the objects

patient1_primary_assay <- CreateChromatinAssay(patient1_primary.counts, fragments = frags.patient1_primary)
GBM_patient1_primary <- CreateSeuratObject(patient1_primary_assay, assay = "ATAC", meta.data=md.patient1_primary)


patient2_recu_assay <- CreateChromatinAssay(patient2_recu.counts, fragments = frags.patient2_recu)
GBM_patient2_recu <- CreateSeuratObject(patient2_recu_assay, assay = "ATAC", meta.data=md.patient2_recu)
GBM_patient2_recu


patient1_recu_assay <- CreateChromatinAssay(patient1_recu.counts, fragments = frags.patient1_recu)
GBM_patient1_recu <- CreateSeuratObject(patient1_recu_assay, assay = "ATAC", meta.data=patient1_recu)
GBM_patient1_recu


patient2_primary_assay <- CreateChromatinAssay(patient2_primary.counts, fragments = frags.patient2_primary)
GBM_patient2_primary <- CreateSeuratObject(patient2_primary_assay, assay = "ATAC", meta.data=md.patient2_primary)
GBM_patient2_primary


# add information to identify dataset of origin

#GBM patient1_primary
GBM_patient1_primary$dataset <- 'patient1_primary_single'
GBM_patient1_primary <- AddMetaData(object = GBM_patient1_primary, metadata = "patient1_primary", col.name = "Donor")
GBM_patient1_primary <- AddMetaData(object = GBM_patient1_primary, metadata = "primary", col.name = "Tumor")
GBM_patient1_primary <- AddMetaData(object = GBM_patient1_primary, metadata = "Donor1", col.name = "Patient")

#amulet
multi <- read.table("/../amulet/MultipletBarcodes_01.txt", sep = "\t", header = FALSE)

cell_ids <- as.data.frame(WhichCells(GBM_patient1_primary))

colnames(multi) <- c("cell_barcodes")
colnames(cell_ids) <- c("cell_barcodes")
#Identify common characters
common_chars <- intersect(multi$cell_barcodes, cell_ids$cell_barcodes)
cell_ids$doublet <- ifelse(cell_ids$cell_barcodes %in% common_chars, "multiplet", "singlet")
rownames(cell_ids) <- cell_ids$cell_barcodes
cell_ids <- cell_ids[, -1]
#adding it to the meta data
GBM_patient1_primary <- AddMetaData(GBM_patient1_primary, metadata = cell_ids, col.name = "doublets")
GBM_patient1_primary <- subset(GBM_patient1_primary, subset = doublets == "singlet")
GBM_patient1_primary
unique(GBM_patient1_primary$doublets)

#meta.data from single processed data
meta.data <- read.csv(file.path("/../output/meta.data.csv"),row.names=1)
cell_ids <- rownames(meta.data)
GBM_patient1_primary <- subset(GBM_patient1_primary, cells = cell_ids)
# Filter the data frame    
GBM_patient1_primary <- AddMetaData(object = GBM_patient1_primary, metadata = meta.data$mtDNA_depth, col.name = "mtDNA_depth")


#GBM patient2_recu
GBM_patient2_recu$dataset <- 'patient2_recu'
GBM_patient2_recu <- AddMetaData(object = GBM_patient2_recu, metadata = "patient2_recu", col.name = "Donor")
GBM_patient2_recu <- AddMetaData(object = GBM_patient2_recu, metadata = "relapse", col.name = "Tumor")
GBM_patient2_recu <- AddMetaData(object = GBM_patient2_recu, metadata = "Donor2", col.name = "Patient")

multi <- read.table("/../data/amulet/MultipletBarcodes_01.txt", sep = "\t", header = FALSE)

cell_ids <- as.data.frame(WhichCells(GBM_patient2_recu))

colnames(multi) <- c("cell_barcodes")
colnames(cell_ids) <- c("cell_barcodes")

# Identify common characters
common_chars <- intersect(multi$cell_barcodes, cell_ids$cell_barcodes)

cell_ids$doublet <- ifelse(cell_ids$cell_barcodes %in% common_chars, "multiplet", "singlet")
rownames(cell_ids) <- cell_ids$cell_barcodes
cell_ids <- cell_ids[, -1]

GBM_patient2_recu <- AddMetaData(GBM_patient2_recu, metadata = cell_ids, col.name = "doublets")
GBM_patient2_recu <- subset(GBM_patient2_recu, subset = doublets == "singlet")
GBM_patient2_recu
unique(GBM_patient2_recu$doublets)

#meta.data from single processed data
meta.data <- read.csv(file.path("/../output/meta.data.csv"),row.names=1)

cell_ids <- rownames(meta.data)
GBM_patient2_recu <- subset(GBM_patient2_recu, cells = cell_ids) 
GBM_patient2_recu <- AddMetaData(object = GBM_patient2_recu, metadata = meta.data$mtDNA_depth, col.name = "mtDNA_depth")



#GBM_patient1_recu
GBM_patient1_recu$dataset <- 'patient1_recu'
GBM_patient1_recu <- AddMetaData(object = GBM_patient1_recu, metadata = "patient1_recu", col.name = "Donor")
GBM_patient1_recu <- AddMetaData(object = GBM_patient1_recu, metadata = "relapse", col.name = "Tumor")
GBM_patient1_recu <- AddMetaData(object = GBM_patient1_recu, metadata = "Donor1", col.name = "Patient")

multi <- read.table("/../amulet/MultipletBarcodes_01.txt", sep = "\t", header = FALSE)

cell_ids <- as.data.frame(WhichCells(GBM_patient1_recu))

colnames(multi) <- c("cell_barcodes")
colnames(cell_ids) <- c("cell_barcodes")

# Identify common characters
common_chars <- intersect(multi$cell_barcodes, cell_ids$cell_barcodes)

cell_ids$doublet <- ifelse(cell_ids$cell_barcodes %in% common_chars, "multiplet", "singlet")
rownames(cell_ids) <- cell_ids$cell_barcodes
cell_ids <- cell_ids[, -1]

GBM_patient1_recu <- AddMetaData(GBM_patient1_recu, metadata = cell_ids, col.name = "doublets")
GBM_patient1_recu <- subset(GBM_patient1_recu, subset = doublets == "singlet")
GBM_patient1_recu
unique(GBM_patient1_recu$doublets)

#meta.data from single processed data
meta.data <- read.csv(file.path("/../output/meta.data.csv"),row.names=1)

cell_ids <- rownames(meta.data)
GBM_patient1_recu <- subset(GBM_patient1_recu, cells = cell_ids)   
GBM_patient1_recu <- AddMetaData(object = GBM_patient1_recu, metadata = meta.data$mtDNA_depth, col.name = "mtDNA_depth")


#GBM_patient2_primary
GBM_patient2_primary$dataset <- 'patient2_primary'
GBM_patient2_primary <- AddMetaData(object = GBM_patient2_primary, metadata = "patient2_primary", col.name = "Donor")
GBM_patient2_primary <- AddMetaData(object = GBM_patient2_primary, metadata = "primary", col.name = "Tumor")
GBM_patient2_primary <- AddMetaData(object = GBM_patient2_primary, metadata = "Donor2", col.name = "Patient")

multi <- read.table("/../data/amulet/MultipletBarcodes_01.txt", sep = "\t", header = FALSE)

cell_ids <- as.data.frame(WhichCells(GBM_patient2_primary))

colnames(multi) <- c("cell_barcodes")
colnames(cell_ids) <- c("cell_barcodes")

# Identify common characters
common_chars <- intersect(multi$cell_barcodes, cell_ids$cell_barcodes)

cell_ids$doublet <- ifelse(cell_ids$cell_barcodes %in% common_chars, "multiplet", "singlet")
rownames(cell_ids) <- cell_ids$cell_barcodes
cell_ids <- cell_ids[, -1]

GBM_patient2_primary <- AddMetaData(GBM_patient2_primary, metadata = cell_ids, col.name = "doublets")
GBM_patient2_primary <- subset(GBM_patient2_primary, subset = doublets == "singlet")
GBM_patient2_primary
unique(GBM_patient2_primary$doublets)

#meta.data from single processed data
meta.data <- read.csv(file.path("/../output/meta.data.csv"),row.names=1)

cell_ids <- rownames(meta.data)
GBM_patient2_primary <- subset(GBM_patient2_primary, cells = cell_ids)  
GBM_patient2_primary <- AddMetaData(object = GBM_patient2_primary, metadata = meta.data$mtDNA_depth, col.name = "mtDNA_depth")


# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = GBM_patient1_primary,
  y = list(GBM_patient2_recu, GBM_patient1_recu, GBM_patient2_primary),
  add.cell.ids = c("patient1_primary_single", "patient2_recu_single", "patient1_recu_single", "patient2_primary_single")
)
combined[["ATAC"]]

#ChromatinAssay data with 187802 features for 28612 cells
#Variable features: 0 
#Genome: 
#Annotation present: FALSE 
#Motifs present: FALSE 
#Fragment files: 6 

setwd("/../integration/output")

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
ggsave(plot = TSS, filename = "TSS.png", width = 7.75, height = 4.55, dpi=300)


combined$nucleosome_group <- ifelse(combined$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragHisto <- FragmentHistogram(object = combined, group.by = 'nucleosome_group')
ggsave(plot = FragHisto, filename = "FragHisto.png", width = 7.75, height = 4.55, dpi=300)

# Augment QC metrics that were computed by cellranger-atac
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
combined$pct_reads_in_DNase <- combined$DNase_sensitive_region_fragments / combined$passed_filters * 100
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments

Vln <- VlnPlot(combined, c("TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks", "peak_region_fragments"), pt.size = 0, ncol = 2, cols = "#72ACA7") 
ggsave(plot = Vln, filename = "Vln_beforeQC.png", width = 20, height = 20, dpi=300)

DensSca <- DensityScatter(combined, x = 'peak_region_fragments', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave(plot = DensSca, filename = "DensSca.png", width = 20, height = 20, dpi=300)

combined


#the filtering was already processed in the single experiment as I added the meta data - so no filtering here for the combined data set
#makes also sense as the data quality differs between samples 

Vln <- VlnPlot(combined, c("TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks", "peak_region_fragments", "mtDNA_depth"), pt.size = 0, ncol = 2, cols = "#72ACA7") 
ggsave(plot = Vln, filename = "Vln_afterQC.png", width = 20, height = 20, dpi=300)

combined

DefaultAssay(combined) <- "ATAC"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 10)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = "lsi", dims = 2:50)
combined <- FindNeighbors(combined, reduction = "lsi", dims = 2:50)
combined <- FindClusters(combined, resolution = 0.6, algorithm = 3)


Dim <- DimPlot(combined, group.by="dataset", raster = FALSE) 
ggsave(plot = Dim, filename = "Dim3_afterQC.png", width = 7.75, height = 4.55, dpi=300)


write.csv(meta.data, "meta.data_merged.csv")

saveRDS(combined, "merged.rds")


