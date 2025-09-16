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
library(svglite)
library(MetBrewer)

#make directories accordingly
#functions
'%ni%' <- Negate('%in%')

setwd(".../output/")

df <- readRDS(file="preprocessed.rds")
variants <- readRDS("variants_ncof5.rds")

# Subset to HQ cells that exist so far
variants <- variants[,colnames(df)]

# top50 variants by mean heteroplasmy defined in ComplexHeatmap
mutation_columns <- c( "16390G>A", "2221C>T" , "16391G>A", "5540G>A",  "16291C>T", "13676A>G", "7362G>A",  "11711G>A", "5703G>A" , "15755T>C", "5979G>A", 
                       "16292C>T", "189A>G" ,  "709G>A" ,  "2992G>A"  ,"5112G>A"  ,"2702G>A",  "5814T>C" , "5329T>C" , "8141G>A"  ,"564G>A"   ,"10310G>A",
                       "66G>A"  ,  "7950T>C" , "2571G>A",  "5553T>C"  ,"3591G>A",  "16366C>T" ,"16261C>T" ,"879T>C" ,  "385A>G" ,  "16293A>G" ,"195T>C"  ,
                       "10373G>A" ,"16327C>T" ,"12439T>C", "14384G>A" ,"368A>G" ,  "199T>C"  , "10290G>A", "3380G>A" , "1485G>A" , "249A>G"  , "15729T>C",
                       "12093T>C", "16223C>T", "16166A>G", "294T>C" ,  "2275T>C" , "6943T>C")


assay_data <- assay(variants, "allele_frequency") 
assay_data

# Subset the assay_data to keep only the rows of interest
assay_data <- assay_data[rownames(assay_data) %in% mutation_columns, ]

# Add the new assay to the Seurat object
df[["alleles"]] <- CreateAssayObject(counts = assay_data)

alleles_counts <- df@assays[["alleles"]]@counts
head(alleles_counts)

# Extract metadata (cell types)
metadata_df <- df@meta.data[, "peaks_snn_res.0.5", drop = FALSE]
head(metadata_df)

metadata_df <- metadata_df %>%
  mutate(Celltype = case_when(
    peaks_snn_res.0.5 == "0" ~ "SMC",
    peaks_snn_res.0.5 == "1" ~ "SMC I",
    peaks_snn_res.0.5 == "2" ~ "Fibroblasts",
    peaks_snn_res.0.5 == "3" ~ "Macrophages",
    peaks_snn_res.0.5 == "4" ~ "NK/T Cells",
    peaks_snn_res.0.5 == "5" ~ "Endothelial Cells",
    TRUE ~ "Unknown"
  ))
head(metadata_df)


clustering.color <- c("SMC" = "deepskyblue3",
                      "SMC I" = "skyblue2",
                      "Fibroblasts" = "darkseagreen3",
                      "Macrophages" = "plum2",
                      "NK/T Cells" = "lightsalmon",
                      "Endothelial Cells" = "brown2")

# 1. Ensure column names of mat match rownames of df
stopifnot(all(colnames(alleles_counts) %in% rownames(metadata_df)))

# 2. Reorder matrix columns to match data frame rownames
mat_ordered <- alleles_counts[, rownames(metadata_df)]

# (Optional) 3. Combine — e.g., if you want to turn the matrix into a data frame and cbind it
# (assuming mat rows match something meaningful to df like features or variants)
alleles_counts_df <- as.data.frame(t(mat_ordered))  # transpose to make samples as rows
head(alleles_counts_df)

new_column_names <- colnames(alleles_counts_df)
new_column_names <- ifelse(grepl("^[0-9]", new_column_names), paste0("X", new_column_names), new_column_names)
colnames(alleles_counts_df) <- new_column_names
head(alleles_counts_df)

# Identify columns with at least one value > 0.2
rows_to_keep <- rowSums(alleles_counts_df > 0.2) > 0
rows_to_keep
# Subset the matrix to keep only relevant columns
alleles_counts_df <- alleles_counts_df[rows_to_keep, ]

head(metadata_df)
metadata_df <- metadata_df[rows_to_keep, ]


# Combine with df (side-by-side)
mtDNA_vars_celltype <- cbind(metadata_df, alleles_counts_df)
mtDNA_vars_celltype$patient <- "hAorta"


patient <- "hAorta"

dir_csv <- (".../output/")
dir_plots <- (".../output/")

head(mtDNA_vars_celltype)

head(mtDNA_vars_celltype)
###############
## lineage ##
###############

type <- "Celltype"


if(type %in% c("Celltype")){
  #calculate lineage bias using the kruskal test
  bias <- data.frame(Celltype = NA, mut = NA, bias = NA)
  for(Celltype in unique(mtDNA_vars_celltype$Celltype)){
    y <- as.numeric(mtDNA_vars_celltype$Celltype == Celltype)
    for(mut in colnames(mtDNA_vars_celltype[grepl("X",colnames(mtDNA_vars_celltype))])){
      tmp <- data.frame(Celltype,mut,bias = kruskal.test(mtDNA_vars_celltype[[mut]] ~ y)[[3]])
      bias <- rbind(tmp,bias)
    }
  }
  print("Lineage bias calculated")
} 
  

#calculate the adjusted p value
bias$kruskal_pvalue_adj <- p.adjust(bias$bias, method = "BH")
  

#evaluate the heteroplasmy of each variant per cell type (not per dataset as this might be misleading)
  heteroplasmy_celltype <- sapply(unique(mtDNA_vars_celltype$Celltype), function(x) {colMeans(mtDNA_vars_celltype[grepl("X",colnames(mtDNA_vars_celltype))][which(mtDNA_vars_celltype$Celltype == x),])})
  heteroplasmy_celltype <- reshape2::melt(heteroplasmy_celltype, value.name = "Heteroplasmy_celltype")
  heteroplasmy_celltype$Celltype <- heteroplasmy_celltype$Var2
  heteroplasmy_celltype$mut <- heteroplasmy_celltype$Var1
  heteroplasmy_celltype$Var1 <- NULL
  heteroplasmy_celltype$Var2 <- NULL
  
  #check number of cells with a heteroplasmy over 10% per celltype 
  n_cells_per_celltype <- sapply(unique(mtDNA_vars_celltype$Celltype), function(x) {colSums(mtDNA_vars_celltype[grepl("X",colnames(mtDNA_vars_celltype))][which(mtDNA_vars_celltype$Celltype == x),]>0.2)})
  n_cells_per_celltype <- reshape2::melt(n_cells_per_celltype, value.name = "n_cells_celltype")
  n_cells_per_celltype$Celltype <- n_cells_per_celltype$Var2
  n_cells_per_celltype$mut <- n_cells_per_celltype$Var1
  n_cells_per_celltype$Var1 <- NULL
  n_cells_per_celltype$Var2 <- NULL
  
  #merge datasets
  heteroplamy_number <- merge(heteroplasmy_celltype, n_cells_per_celltype, by= c("Celltype","mut"))
  
  complete_bias_df <- merge(bias, heteroplamy_number, by= c("Celltype","mut"))
  
  #write csv
  write.csv(complete_bias_df, paste0(dir_csv, "Celltype_bias.csv"))
  print("Saved file for Celltypes")
  
  #plot the data
  pMut_mean <- ggplot(complete_bias_df, aes(x = -log10(kruskal_pvalue_adj), y = Heteroplasmy_celltype*100, col = Celltype, label = mut)) +
    geom_point(aes(size = n_cells_celltype), alpha = 0.5) + 
    #scale_size_continuous(range = c(1,50))+
    scale_y_log10(breaks = c(0.01,0.01, 0.05, 0.1, 0.2,0.5)*100) +
    #geom_text_repel(aes(label=ifelse(-log10(kruskal_pvalue_adj)>0.9*max(-log10(kruskal_pvalue_adj)),as.character(mut),'')),hjust=0,vjust=0, size = 35/.pt) +
    #labs(x = "-log10 p cluster bias", y = "Pseudobulk AF% (log10 scale)") +
    #geom_text()+
    scale_color_manual(values = c(clustering.color)) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.2)
    ) +
    #theme(strip.background =element_rect(fill="white"))+
    #xlim(0,max(-log10(complete_bias_df$kruskal_pvalue_adj)))+
    labs(x ="Cluster Bias", y = "Bulk Heteroplasmy per Celltype")+ 
    theme(axis.title.x = element_text(color="black", size=5, face="bold"),
          axis.title.y = element_text(color="black", size=5, face="bold"),
          axis.text = element_text(size = 5))
  
  cowplot::ggsave2(pMut_mean, file = paste0(dir_plots,"_lineage_bias.pdf"), width = 5, height = 3, limitsize = FALSE)
  pMut_mean
  print(paste0(" Saved plot for lineages"))

 