library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(BuenColors)


# The workflow shown here is representative for all Datasets


setwd("...")
output_dir <- c("...")

df <- readRDS(file=".../integrated_annotated.rds")

cnv = read.table(".../epiAneufinder_results/results_table.tsv")

cnv = cnv[cnv$seq %in% paste0("chr", as.character(1:22)),]

df_1 <- subset(df, subset = dataset == "patient1_primary")

# Check the current cell names
head(Cells(df_1))
new_cell_names <- gsub(pattern = "^patient1_primary_single_", replacement = "", x = Cells(df_1))
# Set the modified cell names back to the Seurat object
df_1 <- RenameCells(object = df_1, new.names = new_cell_names)
# Check the modified cell names
head(Cells(df_1))

seq = cnv[, 1:3]
cnv = cnv[, -(1:3)]

colnames(cnv) <- gsub("cell.", "", colnames(cnv))
colnames(cnv) <- gsub(".1", "-1", colnames(cnv))
dim(cnv)

head(rownames(cnv))
cells.sample = Cells(df_1)

good.cells = intersect(colnames(cnv), cells.sample)
cnv = cnv[, good.cells]
dim(cnv)
rownames(cnv) = paste0(seq$seq, "_", seq$start, "-", seq$end)

# Each column represents a CELL 
# The values represent the fraction of CNV values in that column that are either 0 or 2.
score_all = data.frame(Genome_CNVscore = apply(cnv, 2, function(x) sum(x %in% c(0,2))))/nrow(cnv)

#                    Genome_CNVscore
#ACAGCGCCAGAAAGAG-1     0.069869500
#TGCTTCGAGGATGTCG-1     0.045039843
#this means that cell 1 has a proportion of 0.069 of bins containing 0 or 2

scores_chr = sapply(1:22, function(y) {
  apply(cnv[grep(paste0("chr", as.character(y), "_"), rownames(cnv)),], 2, function(x) sum(x %in% c(0,2))/(sum(grepl(paste0("chr", as.character(y), "_"), rownames(cnv)))))
})



colnames(scores_chr) = paste0("Chr", as.character(1:22), "_CNVscore")

#                   Chr1_CNVscore Chr2_CNVscore Chr3_CNVscore Chr4_CNVscore Chr5_CNVscore Chr6_CNVscore Chr7_CNVscore
#ACAGCGCCAGAAAGAG-1   0.000000000    0.03795812    0.00000000   0.020075963    0.00000000     0.4152542    0.00000000
#TGCTTCGAGGATGTCG-1   0.030166592    0.00000000    0.01565030   0.004883342    0.25366999     0.0000000    0.08609272


scores_df = cbind(score_all, scores_chr)
head(scores_df)

meta.data <- df_1@meta.data

head(meta.data)
head(scores_df)

# Reorder df2 to match the rownames of df1
scores_df <- scores_df[match(rownames(meta.data), rownames(scores_df)), ]

# Now you can safely use cbind()
df_combined <- cbind(meta.data, scores_df)

# Print result
head(df_combined)

df_1 = AddMetaData(df_1, scores_df)


# --- 1) Helpers ---------------------------------------------------------------

# Return row indices for a chromosome
chr_rows <- function(chr, cnv_matrix) {
  grep(paste0("^chr", chr, "_"), rownames(cnv_matrix))
}

# Compute per-cell summary stats for one chromosome
summarize_chr <- function(chr, cnv_matrix) {
  r <- chr_rows(chr, cnv_matrix)
  X <- cnv_matrix[r, , drop = FALSE]
  n_bins <- nrow(X)
  
  # counts
  n_loss <- colSums(X == 0, na.rm = TRUE)
  n_neut <- colSums(X == 1, na.rm = TRUE)
  n_gain <- colSums(X == 2, na.rm = TRUE)
  
  # fractions
  f_loss <- n_loss / n_bins
  f_gain <- n_gain / n_bins
  
  # continuous gradients (choose what you prefer):
  # (a) mean-centered value in [-1, +1]
  mean_centered <- colMeans(X - 1, na.rm = TRUE)   # -1=all loss, 0=neutral, +1=all gain
  # (b) net fraction difference in [-1, +1]
  net_fraction  <- f_gain - f_loss                 # -1=all loss, +1=all gain
  
  # discrete call with a sensible default threshold on |net_fraction|
  # tweak 'thr' as you like (e.g., 0.10 or 0.15)
  thr <- 0.10
  call <- ifelse(net_fraction >  thr, "gain",
                 ifelse(net_fraction < -thr, "loss", "neutral"))
  
  data.frame(
    cell = colnames(X),
    chr = paste0("chr", chr),
    n_bins = n_bins,
    n_loss = n_loss, n_neut = n_neut, n_gain = n_gain,
    f_loss = f_loss, f_gain = f_gain,
    mean_centered = mean_centered,
    net_fraction = net_fraction,
    call = call,
    row.names = NULL,
    check.names = FALSE
  )
}

# --- 2) Summaries for chr7 & chr10 -------------------------------------------

sum_chr7  <- summarize_chr(7,  cnv)
sum_chr10 <- summarize_chr(10, cnv)

# Merge side-by-side for quick inspection
cnv_chr7_chr10 <- merge(
  sum_chr7[,  c("cell","mean_centered","net_fraction","call")],
  sum_chr10[, c("cell","mean_centered","net_fraction","call")],
  by = "cell", suffixes = c("_chr7", "_chr10")
)

head(cnv_chr7_chr10)
rownames(cnv_chr7_chr10) <- cnv_chr7_chr10$cell

df_1 = AddMetaData(df_1, cnv_chr7_chr10)
df_1

meta.data <- df_1@meta.data

# order so high values are plotted last (on top)
meta.data <- meta.data[order(meta.data$mean_centered_chr7), ]

# Split data into NA and non-NA
df_na    <- meta.data[is.na(meta.data$mean_centered_chr7), ]
df_nonna <- meta.data[!is.na(meta.data$mean_centered_chr7), ]



p <- ggplot() +
  # background layer: NA cells in grey
  geom_point(data = df_na, aes(x = UMAP_1, y = UMAP_2),
             color = "grey95", size = 1, alpha = 0.9) +
  # foreground: colored cells
  geom_point(data = df_nonna, aes(x = UMAP_1, y = UMAP_2, color = mean_centered_chr7),
             size = 0.5) +
  scale_color_gradientn(colors = rev(jdb_palette("solar_extra"))) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    axis.line = element_blank()
  )

p
ggsave(p, file = paste0(output_dir, "CNV_7.pdf"),
       width = 5, height = 5, limitsize = FALSE)



# order so LOW values are plotted last (on top)
meta.data <- meta.data[order(meta.data$mean_centered_chr10, decreasing = TRUE), ]

# Split data into NA and non-NA
df_na    <- meta.data[is.na(meta.data$mean_centered_chr10), ]
df_nonna <- meta.data[!is.na(meta.data$mean_centered_chr10), ]


p10 <- ggplot() +
  # background: NA cells
  geom_point(data = df_na, aes(x = UMAP_1, y = UMAP_2),
             color = "grey95", size = 1, alpha = 0.9) +
  # foreground: ordered so lows are drawn on top
  geom_point(data = df_nonna, aes(x = UMAP_1, y = UMAP_2, color = mean_centered_chr10),
             size = 0.5) +
  scale_color_gradientn(colors =  rev(jdb_palette("solar_extra"))) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    axis.line = element_blank()
  )

p10
ggsave(p10, file = paste0(output_dir, "CNV_10.pdf"),
       width = 5, height = 5, limitsize = FALSE)

