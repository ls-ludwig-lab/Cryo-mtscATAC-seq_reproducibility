# Load required libraries
library(Signac)              
library(Seurat)              
library(ggplot2)             
library(SummarizedExperiment) 
library(Matrix)              
library(data.table)          


# Setup
set.seed(257)                 # Reproducibility
"%ni%" <- Negate("%in%")      # Define custom 'not in' operator

# ============================================
# Define mitochondrial DNA variant calling function
# ============================================
call_mutations_mgatk <- function(SE, stabilize_variance = TRUE, low_coverage_threshold = 10) {
  
  # --- Extract key coverage and allele info ---
  cov <- assays(SE)[["coverage"]]                     # Coverage matrix
  ref_allele <- toupper(as.character(rowRanges(SE)$refAllele)) # Reference alleles
  
  # --- Process mutation for one alternate base (A, C, G, T) ---
  process_letter <- function(letter) {
    print(letter)
    
    # Filter positions where reference allele is not the alternate letter
    boo <- ref_allele != letter & ref_allele != "N"
    pos <- start(rowRanges(SE))
    variant_name <- paste0(as.character(pos), ref_allele, ">", letter)[boo]
    nucleotide   <- paste0(ref_allele, ">", letter)[boo]
    position_filt <- pos[boo]
    
    # --- Single-cell functions ---
    getMutMatrix <- function(letter) {
      mat <- ((assays(SE)[[paste0(letter, "_counts_fw")]] + 
                 assays(SE)[[paste0(letter, "_counts_rev")]]) / cov)[boo, ]
      rownames(mat) <- variant_name
      as(mat, "dgCMatrix")
    }
    
    # Bulk allele frequency across all cells
    getBulk <- function(letter) {
      vec <- (Matrix::rowSums(assays(SE)[[paste0(letter, "_counts_fw")]] +
                                assays(SE)[[paste0(letter, "_counts_rev")]]) /
                Matrix::rowSums(cov))[boo]
      vec
    }
    
    # Row variance helper
    rowVars <- function(x, ...) {
      Matrix::rowSums((x - Matrix::rowMeans(x, ...))^2, ...) / (dim(x)[2] - 1)
    }
    
    # Replace NA/NaN with zeros
    update_missing_w_zero <- function(vec) {
      ifelse(is.na(vec) | is.nan(vec), 0, vec)
    }
    
    # --- Strand correlation per mutation ---
    dt <- merge(
      data.table(Matrix::summary(assays(SE)[[paste0(letter, "_counts_fw")]][boo, ])),
      data.table(Matrix::summary(assays(SE)[[paste0(letter, "_counts_rev")]][boo, ])),
      by.x = c("i", "j"), by.y = c("i", "j"),
      all = TRUE
    )[x.x > 0 | x.y > 0]
    
    dt$x.x <- update_missing_w_zero(dt$x.x)
    dt$x.y <- update_missing_w_zero(dt$x.y)
    
    dt2 <- data.table(
      variant  = variant_name[dt[[1]]],
      cell_idx = dt[[2]],
      forward  = dt[[3]],
      reverse  = dt[[4]]
    )
    
    cor_dt <- dt2[, .(cor = cor(forward, reverse, method = "pearson", use = "pairwise.complete")),
                  by = variant]
    
    cor_vec_val <- cor_dt$cor
    names(cor_vec_val) <- as.character(cor_dt$variant)
    
    # --- Single-cell allele frequency matrix ---
    mat <- getMutMatrix(letter)
    mmat <- sparseMatrix(
      i = c(summary(mat)$i, dim(mat)[1]),
      j = c(summary(mat)$j, dim(mat)[2]),
      x = c(update_missing_w_zero(summary(mat)$x), 0)
    )
    
    # --- Bulk mean frequency ---
    mean <- update_missing_w_zero(getBulk(letter))
    
    # --- Variance stabilization (replace low-coverage cells with variant mean) ---
    if (stabilize_variance) {
      idx_mat <- which(data.matrix(cov[boo, ] < low_coverage_threshold), arr.ind = TRUE)
      idx_mat_mean <- mean[idx_mat[, 1]]
      
      ones <- 1 - sparseMatrix(
        i = c(idx_mat[, 1], dim(mmat)[1]),
        j = c(idx_mat[, 2], dim(mmat)[2]),
        x = 1
      )
      
      means_mat <- sparseMatrix(
        i = c(idx_mat[, 1], dim(mmat)[1]),
        j = c(idx_mat[, 2], dim(mmat)[2]),
        x = c(idx_mat_mean, 0)
      )
      
      mmat2 <- mmat * ones + means_mat
      variance <- rowVars(mmat2)
      rm(mmat2, ones, means_mat, idx_mat, idx_mat_mean)
    } else {
      variance <- rowVars(mmat)
    }
    
    # --- Detection statistics ---
    detected <- (assays(SE)[[paste0(letter, "_counts_fw")]][boo, ] >= 2) +
      (assays(SE)[[paste0(letter, "_counts_rev")]][boo, ] >= 2)
    
    # --- Summarize variant-level statistics ---
    var_summary_df <- data.frame(
      position            = position_filt,
      nucleotide          = nucleotide,
      variant             = variant_name,
      vmr                 = variance / (mean + 1e-11),
      mean                = round(mean, 7),
      variance            = round(variance, 7),
      n_cells_conf_detected = Matrix::rowSums(detected == 2),
      n_cells_over_5      = Matrix::rowSums(mmat >= 0.05),
      n_cells_over_10     = Matrix::rowSums(mmat >= 0.10),
      n_cells_over_20     = Matrix::rowSums(mmat >= 0.20),
      strand_correlation  = cor_vec_val[variant_name],
      mean_coverage       = Matrix::rowMeans(cov)[boo],
      stringsAsFactors    = FALSE,
      row.names           = variant_name
    )
    
    # Return as new SummarizedExperiment
    SummarizedExperiment(
      rowData = var_summary_df,
      colData = colData(SE),
      assays  = list(allele_frequency = mmat, coverage = cov[boo, ])
    )
  }
  
  # Run for all possible alternate bases
  SummarizedExperiment::rbind(
    process_letter("A"),
    process_letter("C"),
    process_letter("G"),
    process_letter("T")
  )
}

# ============================================
# Import mgatk output
# ============================================
NB01_mgatk <- readRDS("cryo_mtscATAC_NB01.rds")
NB02_mgatk <- readRDS("cryo_mtscATAC_NB02.rds")
GN01_mgatk <- readRDS("cryo_mtscATAC_GN01.rds")

# Run variant calling
NB01_mgatk <- call_mutations_mgatk(NB01_mgatk)
NB02_mgatk <- call_mutations_mgatk(NB02_mgatk)
GN01_mgatk <- call_mutations_mgatk(GN01_mgatk)

# ============================================
# Load Seurat objects
# ============================================
NB <- readRDS("path_to/Harmony_Integrated_Object.rds")
GN <- readRDS("path_to/GN01.rds")

# ============================================
# Match variant calling output to Seurat objects
# ============================================

# --- NB01 ---
df <- subset(NB, Sample == "NB01")
new_cell_names <- gsub("^NB01_", "", Cells(df))
df <- RenameCells(df, new.names = new_cell_names)
NB01_mgatk <- NB01_mgatk[, colnames(df)]
saveRDS(NB01_mgatk, "Variant_Calling_NB01.rds")

# --- NB02 ---
df <- subset(NB, Sample == "NB02")
new_cell_names <- gsub("^NB02_", "", Cells(df))
df <- RenameCells(df, new.names = new_cell_names)
NB02_mgatk <- NB02_mgatk[, colnames(df)]
saveRDS(NB02_mgatk, "Variant_Calling_NB02.rds")

# --- GN01 ---
overlap_cells <- intersect(colnames(GN01_mgatk), colnames(GN))
GN01_mgatk <- GN01_mgatk[, overlap_cells]
saveRDS(GN01_mgatk, "Variant_Calling_GN01.rds")










