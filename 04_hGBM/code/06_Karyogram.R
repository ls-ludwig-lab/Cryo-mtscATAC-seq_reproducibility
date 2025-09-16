library(Signac)
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(data.table)
library(GenomicRanges)
library(tidyverse)
library(ComplexHeatmap)
library(tidyheatmaps)
library(tibble)
library(circlize)

# The workflow shown here is representative for all Datasets


setwd(".../output")

df <- readRDS(file=".../integrated_annotated.rds")
df

df1 <- subset(df, subset = dataset == "patient1_primary")

head(Cells(df1))
# Modify the cell names to remove the "patient1_primary_single_" prefix
new_cell_names <- gsub(pattern = "^patient1_primary_single_", replacement = "", x = Cells(df1))
# Set the modified cell names back to the Seurat object
df1 <- RenameCells(object = df1, new.names = new_cell_names)
# Check the modified cell names
head(Cells(df1))
new_cell_names

meta.data <- as.data.frame(df1@meta.data)
head(meta.data)
dim(meta.data)
#[1] 6670   55
colnames(meta.data)


#load EpiAneufinder results
karyo <- as.matrix(read.table('/.../epiAneufinder_results/results_table.tsv'))

dim(karyo)
#look at the dimension of the epiAneufinder and change the colnames 
#also you will not need the first rows with meta data for the matrix --> delete
#[1] 25970  4090

colnames(karyo) <- gsub("cell.", "", colnames(karyo))
colnames(karyo) <- gsub(".1", "-1", colnames(karyo))
head(colnames(karyo))
dim1 <- dim(karyo)
dim1

class(karyo)
str(karyo)

#transforms the data frame
d1 <-as.data.frame(t(karyo)) 
# extracts the last row of d data frame (contains the bin row)
last_row <- d1[nrow(d1), ]
#binds both data frames
d1 <- rbind(last_row, d1)
#deletes last row, as the bin information is now on top, dont need it in the last row from d1 anymore
d1 <- d1[-nrow(d1), ]

#deletes the information regarding bin etc
d2 <- d1[-c(1,2,3,4), ]

#adds a column with the barcodes at the end of the table with
d2$cell_barcodes <- rownames(d2)

#extracts the cell barcoes
t3 <- d2$cell_barcodes

#changes the rownames to the cell barcodes 
rownames(d2) <- d2$cell_barcodes

meta.data$cell_ids <- rownames(meta.data)
#identifies common chracters of the cell barcodes in the karyogram data frame and from cell ids of the the seurat object
common_chars <- intersect(d2$cell_barcodes, meta.data$cell_ids)


# filters the row names in d2 based on the row names of cell_ids (seurat object)
d3 <- d2 %>%
  dplyr::filter(row.names(.) %in% row.names(meta.data))

# filters the row names in information (seurat object) based on the row names of the karyogram
meta.data1 <- meta.data %>%
  dplyr::filter(row.names(.) %in% row.names(d2))

dim(d3)


#depends on the reuslts output - need to change accordingly
#here check tht the cellbarcode and the rownames are the same (due to the data wrangling in the other steps)
head(d3[, 25975:25978])
#                   25975 25976 25977      cell_barcodes
#ACAGCGCCAGAAAGAG-1     1     1     1 ACAGCGCCAGAAAGAG-1
#TGCTTCGAGGATGTCG-1     1     1     1 TGCTTCGAGGATGTCG-1
#CTTGTCGGTAACCCAT-1     1     1     1 CTTGTCGGTAACCCAT-1
#CCTTGGTTCTTTATCG-1     1     1     1 CCTTGGTTCTTTATCG-1

d3 <- d3[,-c(25978)]

#matching the rows based on the rownames 
d3 <- d3[match(row.names(meta.data1), row.names(d3)), ]

t4 <- cbind(meta.data1, d3)

first_column <- t4 %>% dplyr::select(10)

d3 <- as.matrix(d3)

d3 <- apply(d3, 2, as.numeric)
dim(d3)
#[1]  3226 25977
d4 <- as.data.frame(t(d1[-c(5:3254), ]))

dim(d4)


desired_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                   "chr20", "chr21", "chr22")

d4$seq <- factor(d4$seq, levels = desired_order)
column_order <- which(levels(d4$seq) %in% d4$seq)

column_ha = HeatmapAnnotation(Chromosome = d4$seq,
                              show_legend = c(FALSE),
                              col = list(Chromosome = c("chr1" = "grey90",
        "chr2" = "grey90",
        "chr3" = "grey90",
        "chr4" = "grey90",
        "chr5" = "grey90",
        "chr6" = "grey90",
        "chr7" = "grey90",
        "chr8" = "grey90",
        "chr9" = "grey90",
        "chr10" = "grey90",
        "chr11" = "grey90",
        "chr12" = "grey90",
        "chr13" = "grey90",
        "chr14" = "grey90",
        "chr15" = "grey90",
        "chr16" = "grey90",
        "chr17" = "grey90",
        "chr18" = "grey90",
        "chr19" = "grey90",
        "chr20" = "grey90",
        "chr21" = "grey90",
        "chr22" = "grey90")))

dim(d3)
dim(t4)

abundance <- table(t4$possible_Celltype_broader)
t4 <- t4 %>%
  mutate(mtDNA_log_transformed = log(mtDNA_depth), base =10)

col_fun_mtDNA = colorRamp2(c(5,100), hcl_palette = "Viridis")
col_fun = colorRamp2(c(0,1, 2), hcl_palette = "RdBu")
col_fun = colorRamp2(c(0,1,2), c("red", "white", "blue"))

row_ha = rowAnnotation("cluster" = t4$possible_Celltype_broader,
                      "mtDNA_depth log10" = t4$mtDNA_depth,
                      col = list(
cluster = c("malignant_other" = "#829BD4",
            "Microglia" = "#D9488B",
            "NPC1-OPC-like" = "lightslateblue",
            "OPC-like" = "dodgerblue", 
            "Oligodendrocytes" =  "dodgerblue4",
            "Neuronal Cells" = "orange",
            "AC-MEs-like" = "#A7BAF2",
            "DLX-GAD-high Neuronal Cells" = "sienna2",
            "OPC-like" = "dodgerblue3",
            "AC-like" = "paleturquoise3",
            "Endothelial Cells" = "#400135",
            "Microglia I" = "#F294AD",
            "Pericytes" = "#73026B",
            "T Cells" ="#8C233F"),
"mtDNA_depth log10" = col_fun_mtDNA))


pdf(file="ComplexHeatKaryo_Cluster.pdf", width = 25, height = 20)
Heatmap(d3, 
        show_column_names = FALSE,
        show_row_names = FALSE,
        show_row_dend = FALSE,
        top_annotation = column_ha,
        right_annotation = row_ha,
        row_split=t4$possible_Celltype,
        column_split=d4$seq, 
        cluster_columns = FALSE,
        col = col_fun,
        use_raster = TRUE,
        column_title = "CNVs patient1_primary",
    width = unit(15, "cm"),
    height = unit(16, "cm"),
   raster_quality = 20)
dev.off()
