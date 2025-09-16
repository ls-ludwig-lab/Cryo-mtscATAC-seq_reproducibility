library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

# Inputs
datasets <- list(
  list(file = "filtered_0p5pct.rds", label = "0p5pct")
)

setwd("...")


y_vars <- c("nCount_peaks", "nFeature_peaks", "mtDNA_depth", "TSS.enrichment")

# Helper: read + collect meta for one dataset
read_meta <- function(d) {
  obj <- readRDS(paste0(
    ".../",
    d$file
  ))
  md <- obj@meta.data
  md$dataset <- d$label
  md
}

# Preload all metadata once
meta_all <- bind_rows(lapply(datasets, read_meta))

# Loop over QC metrics and plot violin + inner boxplot, faceted by dataset
for (y_var in y_vars) {
  # keep only rows with finite values for the metric
  dat <- meta_all %>%
    dplyr::filter(is.finite(.data[[y_var]]))
  head(dat)
  p <- ggplot(dat, aes(x = dataset, y = .data[[y_var]])) +
    geom_violin(trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.75) +  # boxplot inside violin
    scale_color_manual(values = "grey60") +
    scale_fill_manual(values = "white") +
    scale_y_log10() +
    labs(x = NULL, y = y_var) +
    theme_classic() +
    theme(
      axis.title.y = element_text(size = 16),
      axis.text.y  = element_text(size = 12),
      axis.text.x  = element_text(size = 11, angle = 45, hjust = 1),
      strip.text   = element_text(size = 13),
      legend.position = "none"
    )
  
  print(p)
  ggsave(paste0("violin_box_", y_var, "_by_dataset.pdf"), plot = p, width = 3, height = 6)
}
