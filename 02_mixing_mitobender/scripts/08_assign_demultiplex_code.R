library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
"%ni%" <- Negate("%in%")

setwd()
source(".../08_demultiplexing_code.R")

output_dir <- ".../rds_files"
df1 <- assign_pull(readRDS(".../mgatk_mixing.rds"))

write.table(df1, file = file.path(output_dir, "f0p5pct_data_meta.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

df1 <- readRDS(readRDS(".../mgatk_mixing.rds"))
df3 <- readRDS(".../mitoBender/mitoBender_mgatk_mixing.rds")


rowRanges(df3) <- rowRanges(df1)
df3 <- assign_pull(df3)
head(df3)
write.table(df3, file = file.path(output_dir, "f0p5pct_mitoBender_data_meta.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
