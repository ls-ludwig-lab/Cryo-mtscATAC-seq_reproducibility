library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

# From manual curation of previous heatmap
# homoplasmic variants
#Donor specific variants identified in Heatmap - are not displayed here due to human genotyping
#Donor specific variants are homoplasmic, haplotype specific variants, furthermore homoplasmic variants could also be germline variants which are family or person specific and therefore due to data protection of patients not displayed


donor0_vars <- c("variant_1", "variant_2", "variant_3","...")
donor1_vars <- c("variant_4", "variant_5", "variant_6","...")

assign_pull <- function(SE){
  cov <- assays(SE)[["coverage"]]+ 0.001
  ref_allele <- as.character(rowRanges(SE)$refAllele)
  
  getMutMatrix <- function(letter){
    mat <- (assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov
    rownames(mat) <- paste0(as.character(1:dim(mat)[1]), toupper(ref_allele), ">", letter)
    return(mat[toupper(ref_allele) != letter,])
  }
  
  mmat <- rbind(getMutMatrix("A"), getMutMatrix("C"), getMutMatrix("G"), getMutMatrix("T"))
  df <- data.frame(
    cell_id = colnames(mmat),
    d0_mean = round(colMeans(mmat[donor0_vars,]), 3),
    d1_mean = round(colMeans(mmat[donor1_vars,]), 3),
    mean_cov = round(colMeans(cov), 1)
  )
  df$assign <- ifelse(
    df$d0_mean > 0.95, "donor0", 
    ifelse(df$d1_mean > 0.95, "donor1", 
           ifelse(df$mean_cov < 5, "Low_coverage", "Collision")))
  df
  
}

