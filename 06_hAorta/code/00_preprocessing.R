# packages 
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)

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
  genome = 'hg38',
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

 
setwd(".../output")


ggsave(plot = FragHisto, filename = "FragHisto.png", width = 5, height = 5, dpi=300)

df$pct_reads_in_peaks <- df$peak_region_fragments / df$passed_filters * 100
df$pct_reads_in_DNase <- df$DNase_sensitive_region_fragments / df$passed_filters * 100
df$blacklist_ratio <- df$blacklist_region_fragments / df$peak_region_fragments
Vln <- VlnPlot(df, c("TSS.enrichment", "nCount_peaks", "nucleosome_signal", "pct_reads_in_peaks"), pt.size = 0, ncol = 2, cols = "#72ACA7") 
Vln
ggsave(plot = Vln, filename = "Vln.png", width = 8, height = 8, dpi=300)


DensSca <- DensityScatter(df, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
DensSca
ggsave(plot = DensSca, filename = "DensSca.png", width = 20, height = 20, dpi=300)


# remove low-quality cells
df1 <- subset(df,
              subset = nCount_peaks > 1000 &
                TSS.enrichment > 2 & 
                TSS.enrichment < 10 & 
                nucleosome_signal < 2
)
df1

Vln <- VlnPlot(df1, c("TSS.enrichment", "nCount_peaks", "nucleosome_signal", "pct_reads_in_peaks"), pt.size = 0, ncol = 2, cols = "#72ACA7")
Vln
ggsave(plot = Vln, filename = "Vln_afterQC.png", width = 20, height = 20, dpi=300)

#An object of class Seurat 
#133867 features across 9912 samples within 1 assay 
#Active assay: peaks (133867 features, 0 variable features)

#Adding Amulet /Doublet Information
multi <- read.table(".../data/huAorta_A4240806/MultipletBarcodes_01.txt", sep = "\t", header = FALSE)

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
mito.data <- ReadMGATK(dir = ".../data/final")

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
ggsave(plot = p, filename = "Vln_afterQC_mtDNA.png", width = 20, height = 20, dpi=300)

# remove low-quality cells
# filter cells based on mitochondrial depth
df2 <- subset(df1, mtDNA_depth >= 5)
df2
  
#An object of class Seurat 
#233088 features across 4548 samples within 2 assays 
#Active assay: peaks (100536 features, 0 variable features)
#1 other assay present: mito


p <- VlnPlot(df2, features, pt.size = 0, ncol = 3, cols = "#72ACA7") 
p
ggsave(plot = p, filename = "Vln_afterQC_mtDNA_QC.png", width = 20, height = 20, dpi=300)


DefaultAssay(df2) <- "peaks"
df2 <- RunTFIDF(df2)
df2 <- FindTopFeatures(df2, min.cutoff = 10)
df2 <- RunSVD(df2)
df2 <- RunUMAP(df2, reduction = "lsi", dims = 2:50)
df2 <- FindNeighbors(df2, reduction = "lsi", dims = 2:50)
df2 <- FindClusters(df2, resolution = 0.6, algorithm = 3)

p1 <- DimPlot(df2, label.size = 4, label=TRUE) +
  # -- Here you tell ggplot to eliminate all the axis info from your original plot -- #
  theme(aspect.ratio=1/1,
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),#legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  # -- Here you can add the title of the plot -- #
ggtitle("GBM N20_1201")
p1

ggsave(plot = p1, filename = "DimPlot0_6.png", width = 7.75, height = 4.55, dpi=300)


DefaultAssay(df2) <- "peaks"

# compute gene accessibility
gene.activities <- GeneActivity(df2)

# add to the Seurat object as a new assay
df2[['RNA']] <- CreateAssayObject(counts = gene.activities)

df2 <- NormalizeData(
  object = df2,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(df2$nCount_RNA)
)

DefaultAssay(df2) <- "RNA"

saveRDS(df2, "preprocessed.rds")
