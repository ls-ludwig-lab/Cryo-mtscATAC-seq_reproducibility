
# Load required libraries
library(Signac)          
library(Seurat)          
library(ggplot2)         
library(SummarizedExperiment) 
library(dplyr)           


# Setup
set.seed(257)                 # Ensure reproducibility
"%ni%" <- Negate("%in%")      # Define custom 'not in' operator

# Load Variant Calling Objects
NB01_Variants <- readRDS("Variant_Calling_NB01.rds")
NB02_Variants <- readRDS("Variant_Calling_NB02.rds")

# ============================================
# Filter for informative variants
# ============================================
Primary <- NB01_Variants
misc_df_primary <- data.frame(rowData(Primary))
vars_primary <- misc_df_primary %>%
  filter(
    n_cells_conf_detected >= 5,
    strand_correlation > 0.65,
    log10(vmr) > -2,
    mean_coverage >= 5
  ) %>%
  pull(variant)

Second_look <- NB02_Variants
misc_df_second_look <- data.frame(rowData(Second_look))
vars_second_look <- misc_df_second_look %>%
  filter(
    n_cells_conf_detected >= 5,
    strand_correlation > 0.65,
    log10(vmr) > -2,
    mean_coverage >= 5
  ) %>%
  pull(variant)

# ============================================
# Heteroplasmy correlation analysis
# ============================================
# Variants present in either dataset
vars_both <- unique(c(vars_primary, vars_second_look))
prim <- Primary[vars_both, ]
sec  <- Second_look[vars_both, ]

# Annotate whether each variant was called in primary or second dataset
rowData(prim)$called_in_primary  <- rownames(prim) %in% vars_primary
rowData(sec)$called_in_relapse   <- rownames(sec) %in% vars_second_look

# Build dataframe of mean allele frequencies across both datasets
mean_df <- data.frame(
  Primary     = rowData(prim)$mean,
  Second_look = rowData(sec)$mean,
  variant     = rownames(prim)
) %>%
  mutate(log10_FC = log10(Second_look / Primary)) %>%
  arrange(desc(log10_FC))

# Correlations between allele frequencies
cor(log10(mean_df$Second_look), log10(mean_df$Primary))
cor(mean_df$Second_look, mean_df$Primary)

# Prepare for plotting
mean_df$Primary_pct     <- mean_df$Primary * 100
mean_df$Second_look_pct <- mean_df$Second_look * 100
label_data <- subset(mean_df, variant %in% c("953T>C", "16063T>C", "11825G>A"))

# Scatterplot of allele frequencies
ggplot() +
  geom_point(data = mean_df, aes(x = Primary_pct, y = Second_look_pct, color = "grey"), size = 3) +
  geom_point(data = subset(mean_df, variant == "953T>C"),
             aes(x = Primary_pct, y = Second_look_pct, color = "953T>C"), size = 4) +
  geom_point(data = subset(mean_df, variant == "16063T>C"),
             aes(x = Primary_pct, y = Second_look_pct, color = "16063T>C"), size = 4) +
  geom_point(data = subset(mean_df, variant == "11825G>A"),
             aes(x = Primary_pct, y = Second_look_pct, color = "11825G>A"), size = 4) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_text(
    data = label_data,
    aes(x = Primary_pct, y = Second_look_pct, label = variant, color = variant),
    hjust = -0.1, vjust = -0.5, size = 5
  ) +
  scale_color_manual(values = c(
    "953T>C"    = "#76C3A9",
    "16063T>C"  = "#76C3A9",
    "11825G>A"  = "#355C4C",
    "grey"      = "grey"
  )) +
  labs(x = "Primary Heteroplasmy %", y = "Second_look Heteroplasmy %") +
  xlim(0, 2.5) + ylim(0, 2.5) +
  theme(
    panel.background = element_blank(),
    panel.border     = element_blank(),
    axis.text.x      = element_text(size = 16),
    axis.text.y      = element_text(size = 16),
    axis.title.x     = element_text(size = 18),
    axis.title.y     = element_text(size = 18),
    legend.position  = "none",
    axis.line        = element_line(colour = "black", linewidth = 0.9)
  )


# ============================================
# Plot informative variants on UMAP embedding
# ============================================
NB <- readRDS("path_to/Harmony_Integrated_Object.rds")
DefaultAssay(NB) <- "alleles"

# --- Variant 953T>C ---
allele_values <- GetAssayData(NB, assay = "alleles", slot = "data")["953T>C", ]
NB$binary_953T_C <- ifelse(allele_values > 0.2, 1, 0)

FeaturePlot(
  NB,
  features  = "binary_953T_C",
  order     = TRUE,
  cols      = c("grey85", "darkred"),
  ncol      = 3,
  pt.size   = 0.8,
  reduction = "umap",
  split.by  = "Dataset"
) & NoAxes()

# --- Variant 11825G>A ---
allele_values <- GetAssayData(NB, assay = "alleles", slot = "data")["11825G>A", ]
NB$binary_11825G_A <- ifelse(allele_values > 0.2, 1, 0)

FeaturePlot(
  NB,
  features  = "binary_11825G_A",
  order     = TRUE,
  cols      = c("grey85", "darkred"),
  ncol      = 3,
  pt.size   = 0.8,
  reduction = "umap",
  split.by  = "Dataset"
) & NoAxes()









