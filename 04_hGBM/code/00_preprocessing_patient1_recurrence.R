# packages 
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(readxl)
library(RColorBrewer)
# Data loading
# Creating a Seurat Object

# load counts and metadata from cellranger-atac
counts <- Read10X_h5(filename = ".../data/filtered_peak_bc_matrix.h5")

metadata <- read.csv(
  file = ".../singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = ".../fragments.tsv.gz",
  min.cells = 5,
  min.features = 50
)


df <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
df

df[['peaks']]

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(df) <- annotations

# compute nucleosome signal score per cell
df <- NucleosomeSignal(object = df)

# compute TSS enrichment score per cell
df <- TSSEnrichment(object = df, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
df$pct_reads_in_peaks <- df$peak_region_fragments / df$passed_filters * 100
df$blacklist_ratio <- df$blacklist_region_fragments / df$peak_region_fragments
df$high.tss <- ifelse(df$TSS.enrichment > 2, 'High', 'Low')
TSS <- TSSPlot(df, group.by = 'high.tss') + NoLegend()
df$nucleosome_group <- ifelse(df$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragHisto <- FragmentHistogram(object = df, group.by = 'nucleosome_group')
FragHisto

df$pct_reads_in_peaks <- df$peak_region_fragments / df$passed_filters * 100
df$pct_reads_in_DNase <- df$DNase_sensitive_region_fragments / df$passed_filters * 100
df$blacklist_ratio <- df$blacklist_region_fragments / df$peak_region_fragments
Vln <- VlnPlot(df, c("TSS.enrichment", "nCount_peaks", "nucleosome_signal", "pct_reads_in_peaks"), pt.size = 0, ncol = 2, cols = "#72ACA7") 
Vln


DensSca <- DensityScatter(df, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)


# remove low-quality cells
df1 <- subset(df,
              subset = nCount_peaks > 1000 &
                TSS.enrichment > 1.0
)
df1

Vln <- VlnPlot(df1, c("TSS.enrichment", "nCount_peaks", "nucleosome_signal", "pct_reads_in_peaks"), pt.size = 0, ncol = 2, cols = "#72ACA7")

df1

# load mgatk output
mito.data <- ReadMGATK(dir = ".../final")

# create an assay
mito <- CreateAssayObject(counts = mito.data$counts)

# Subset to cell present in the scATAC-seq assat
mito <- subset(mito, cells = colnames(df1))

# add assay and metadata to the seurat object
df1[["mito"]] <- mito
df1 <- AddMetaData(df1, metadata = mito.data$depth, col.name = "mtDNA_depth")


features <- c("TSS.enrichment", "nCount_peaks", "nucleosome_signal", "pct_reads_in_peaks", "mtDNA_depth")
p <- VlnPlot(df1, features, pt.size = 0, ncol = 3, cols = "#72ACA7") 
p

# remove low-quality cells
# filter cells based on mitochondrial depth
df1 <- subset(df1, mtDNA_depth >= 5)
df1


p <- VlnPlot(df1, features, pt.size = 0, ncol = 3, cols = "#72ACA7") 
p

(df1, ".../preprocessed.rds")