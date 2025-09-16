library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(readxl)
library(RColorBrewer)
library(SummarizedExperiment)

# Define color schemes for plotting cell types and dataset regions
output_dir <- ".../Donor1_0p5_single/"

df <- readRDS(file.path(output_dir,file="mitoFilter_df_annotated.rds"))

variants <- readRDS(".../variants_donor1_mitoBender.rds")

misc_df <- data.frame(rowData(variants))

vars <- misc_df %>%
  dplyr::filter(
    n_cells_conf_detected >= 5,
    strand_correlation > 0.65,
    log10(vmr) > -2,
    mean_coverage >= 5
  ) %>%
  pull(variant)
vars
variants <- variants[vars,]
variants

df
# Subset to HQ cells that exist so far
df <- df[,colnames(variants)]
df

assay_data <- assay(variants, "allele_frequency") 
assay_data

# Subset the assay_data to keep only the rows of interest
df[["alleles"]] <- CreateAssayObject(counts = assay_data)


alleles_counts <- df@assays[["alleles"]]@counts
head(alleles_counts)

# Extract metadata (cell types)
metadata_df <- df@meta.data[, "possible_Celltype", drop = FALSE]
head(metadata_df)

DimPlot(df, label = TRUE)

head(metadata_df)
unique(metadata_df$possible_Celltype)

alleles_counts_df <- as.data.frame(t(alleles_counts[, rownames(metadata_df)]))
colnames(alleles_counts_df) <- ifelse(grepl("^[0-9]", colnames(alleles_counts_df)),
                                        paste0("X", colnames(alleles_counts_df)),
                                        colnames(alleles_counts_df))
  
# Merge metadata and alleles into one dataframe
mtDNA_vars_celltype <- cbind(metadata_df, alleles_counts_df)
mtDNA_vars_celltype
  
# Compute Kruskal-Wallis p-values for each variant vs. cell type presence
type <- "possible_Celltype"
bias <- data.frame(possible_Celltype = NA, mut = NA, bias = NA)
bias 

for(possible_Celltype in unique(mtDNA_vars_celltype$possible_Celltype)) {
    y <- as.numeric(mtDNA_vars_celltype$possible_Celltype == possible_Celltype)
    for(mut in grep("^X", colnames(mtDNA_vars_celltype), value = TRUE)) {
      tmp <- data.frame(possible_Celltype, mut,
                        bias = kruskal.test(mtDNA_vars_celltype[[mut]] ~ y)$p.value)
      bias <- rbind(tmp, bias)
    }
  }
  
# Adjust p-values using Benjamini-Hochberg correction
bias$kruskal_pvalue_adj <- p.adjust(bias$bias, method = "BH")
bias <- na.omit(bias)
  
# Calculate mean heteroplasmy per mutation per cell type
celltypes <- unique(mtDNA_vars_celltype$possible_Celltype)
mutation_cols <- grep("^X", colnames(mtDNA_vars_celltype), value = TRUE)
  
heteroplasmy_matrix <- do.call(cbind, lapply(celltypes, function(ct) {
    df_ct <- mtDNA_vars_celltype[mtDNA_vars_celltype$possible_Celltype == ct, mutation_cols, drop = FALSE]
    colMeans(df_ct)
  }))
colnames(heteroplasmy_matrix) <- celltypes


heteroplasmy_celltype <- reshape2::melt(heteroplasmy_matrix, value.name = "Heteroplasmy_celltype")
heteroplasmy_celltype$Celltype <- heteroplasmy_celltype$Var2
heteroplasmy_celltype$mut <- heteroplasmy_celltype$Var1
heteroplasmy_celltype$Var1 <- NULL
heteroplasmy_celltype$Var2 <- NULL
  
# Count number of cells with heteroplasmy > 0.2 per mutation and cell type
n_cells_matrix <- do.call(cbind, lapply(celltypes, function(ct) {
    df_ct <- mtDNA_vars_celltype[mtDNA_vars_celltype$possible_Celltype == ct, mutation_cols, drop = FALSE]
    colSums(df_ct > 0.2)
  }))
colnames(n_cells_matrix) <- celltypes
  
n_cells_per_celltype <- reshape2::melt(n_cells_matrix, value.name = "n_cells_celltype")
n_cells_per_celltype$Celltype <- n_cells_per_celltype$Var2
n_cells_per_celltype$mut <- n_cells_per_celltype$Var1
n_cells_per_celltype$Var1 <- NULL
n_cells_per_celltype$Var2 <- NULL
head(n_cells_per_celltype)
unique(n_cells_per_celltype$Celltype)

colnames(bias)[colnames(bias) == "possible_Celltype"] <- "Celltype"
colnames(heteroplasmy_celltype)[colnames(heteroplasmy_celltype) == "possible_Celltype"] <- "Celltype"
colnames(n_cells_per_celltype)[colnames(n_cells_per_celltype) == "possible_Celltype"] <- "Celltype"

# Merge p-values, heteroplasmy and cell count data into one summary dataframe
complete_bias_df <- merge(bias, merge(heteroplasmy_celltype, n_cells_per_celltype, by=c("Celltype","mut")), by=c("Celltype","mut"))
complete_bias_df

  # Export the summary statistics as CSV for this dataset
write.csv(complete_bias_df, paste0(output_dir, "_Celltype_bias_mitoFilter.csv"))
 
celltype_colors <- c(
  "Oligodendrocytes"    = "deepskyblue3",
  "Oligodendrocytes I" = "deepskyblue2",
  "OPCs"                = "mediumpurple1",
  "Microglia"           = "#F294AD",
  "Lymphocytes" = "plum1",
  "Astrocytes"          = "darkseagreen2",
  "Astrocytes I"          = "darkseagreen3",
  "Astrocyte-like neural stem cells"          = "darkseagreen",
  "Glia Cells"         = "skyblue",
  "Vascular Cells"      = "#400135",  
  "Ependymal Cells"     = "#6B2747",  
  "non neuronal mesoderm" = "#6B2759",
  "Projection Nerons (cortical glutameric lineage)"             = "orange1",
  "Excitatory Neurons (pyramidal, glutamatergic)"             = "orange2"
)
   

unique(complete_bias_df$Celltype)
  
# Create lineage bias scatter plot
p <- ggplot(complete_bias_df, aes(
  x = -log10(kruskal_pvalue_adj),
  y = Heteroplasmy_celltype * 100,
  col = Celltype, label = mut
)) +
  geom_point(aes(size = n_cells_celltype, alpha = 0.5)) +
  ylim(0,80) +
  scale_color_manual(values = celltype_colors) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    #legend.position = "none",
    strip.background = element_rect(fill = "white"),
    axis.line = element_line(color = "black", size = 0.8),  # adds axis lines
    axis.title.x = element_text(color = "black", size = 20, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, face = "bold"),
    axis.text = element_text(size = 15)
  ) +
  xlim(0, 10) +
  labs(x = "Cluster Bias", y = "Bulk Heteroplasmy per Celltype")

p


ggsave(p, file = paste0(output_dir, "Celltype_bias_donor1_mitoFilter_allCells.pdf"),
       width = 20, height = 20, limitsize = FALSE)

