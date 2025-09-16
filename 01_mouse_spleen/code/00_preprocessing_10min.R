# packages 
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Mmusculus.v79)
#BiocManager::install("EnsDb.Mmusculus.v79")


# Data loading
# Creating a Seurat Object
# load counts and metadata from cellranger-atac

setwd(".../outs")

counts <- Read10X_h5(filename = "filtered_peak_bc_matrix.h5")

metadata <- read.csv(
  file = "singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = "fragments.tsv.gz",
  min.cells = 5,
  min.features = 50
)


df <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
df

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

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

setwd("/output")

df$pct_reads_in_peaks <- df$peak_region_fragments / df$passed_filters * 100
df$pct_reads_in_DNase <- df$DNase_sensitive_region_fragments / df$passed_filters * 100
df$blacklist_ratio <- df$blacklist_region_fragments / df$peak_region_fragments
Vln <- VlnPlot(df, c("TSS.enrichment", "nCount_peaks", "nucleosome_signal", "pct_reads_in_peaks"), pt.size = 0, ncol = 2, cols = "#72ACA7") 
Vln


DensSca <- DensityScatter(df, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
DensSca

df
# remove low-quality cells
df1 <- subset(df,
              subset = nCount_peaks > 4000 &
                nCount_peaks < 25000 &
                TSS.enrichment > 2)
df1

Vln <- VlnPlot(df1, c("TSS.enrichment", "nCount_peaks", "nucleosome_signal", "pct_reads_in_peaks"), pt.size = 0, ncol = 2, cols = "#72ACA7")
Vln

#Adding Amulet /Doublet Information
multi <- read.table(".../MultipletBarcodes_01.txt", sep = "\t", header = FALSE)

cell_ids <- as.data.frame(WhichCells(df1))

colnames(multi) <- c("cell_barcodes")
colnames(cell_ids) <- c("cell_barcodes")

# Identify common characters
common_chars <- intersect(multi$cell_barcodes, cell_ids$cell_barcodes)

cell_ids$doublet <- ifelse(cell_ids$cell_barcodes %in% common_chars, "multiplet", "singlet")
rownames(cell_ids) <- cell_ids$cell_barcodes
cell_ids <- cell_ids[, -1]

df1 <- AddMetaData(df1, metadata = cell_ids, col.name = "doublets")
df1

df1 <- subset(df1, subset = doublets == "singlet")
df1
unique(df1$doublets)

# load mgatk output
mito.data <- ReadMGATK(dir = "../final")

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
df2 <- subset(df1, mtDNA_depth >= 5)
df2

saveRDS(df2, "preprocessed.rds")
