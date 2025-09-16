library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(readxl)
library(RColorBrewer)
library(dplyr)
library(tidyr)

setwd(".../integration/output")

d <- readRDS(file=".../integration/output/integrated_annotated.rds")

unique(d$dataset)

# Define the regions you want to subset
datasets <- c("FF_10min", "FF30_min", "Cells")  
unique_idents <- c("FF_10min_", "FF_30min_", "Cells_")  

#Looping with seq_along(datasets)
#Ensures we correctly map each dataset to its corresponding unique_idents value.

#Proper dataset subsetting (dataset == region)
#Filters correctly for each dataset instead of the incorrect dataset == dataset.

#Using the correct prefix in gsub()
#Now removes idents, etc.

# Loop over datasets and their corresponding unique identifiers
for (i in seq_along(datasets)) {
  region <- datasets[i]  # Get current dataset name
  prefix <- unique_idents[i]  # Get corresponding prefix
  
  # Subset the dataset
  df <- subset(d, subset = dataset == region)
  
  # Modify the cell names by removing the specific prefix
  new_cell_names <- gsub(pattern = paste0("^", prefix), replacement = "", x = Cells(df))
  
  # Rename the cells
  df <- RenameCells(object = df, new.names = new_cell_names)
  
  # Assign the modified Seurat object to a variable with the region name
  assign(paste0("df_", gsub("[^a-zA-Z0-9]", "_", region)), df)
  
  # Print first few cell names to check
  print(head(Cells(df)))
}



# Access modified Seurat objects using subset_list[["BB28_5"]], etc.

#coverage mtDNA chromsome

# -- Functions -- #
pull_coverage <- function(SE, cells, resolution = 5){
  zoo::rollmean(rowMeans(assays(SE)[["coverage"]]), resolution)
}

mito.SE <- readRDS(".../fresh_isolated/data/final/CR-ATAC_A4275_mSpleen_freshIsolated_mtMask.rds")

mito.SE <- mito.SE[,colnames(df_Cells)]
# --- 1. Visualize mitochondrial coverage --- #
cov_df <- data.frame(
  pos = zoo::rollmean(1:16299, 5),
  mito = pull_coverage(mito.SE))

cov_df$Sample <- c("Cells")


mito.SE1 <- readRDS(".../FF_10min/data/mgatk_CR-ATAC_mSpleen_FF_10min_mtMask/final/CR-ATAC_mSpleen_FF_10min_mtMask.rds")
mito.SE1 <- mito.SE1[,colnames(df_FF_10min)]
cov_df1 <- data.frame(
  pos = zoo::rollmean(1:16299, 5),
  mito = pull_coverage(mito.SE1))

cov_df1$Sample <- c("FF_10min")


mito.SE2 <- readRDS(".../FF_30min/data/mgatk_CR-ATAC_mSpleen_FF_30min_mtMask/final/CR-ATAC_mSpleen_FF_30min_mtMask.rds")
mito.SE2 <- mito.SE2[,colnames(df_FF30_min)]

cov_df2 <- data.frame(
  pos = zoo::rollmean(1:16299,5),
  mito = pull_coverage(mito.SE2))

cov_df2$Sample <- c("FF_30min")




c3 <- rbind(cov_df,cov_df1, cov_df2)


dataset_color <- c("FF_10min" = "palevioletred2",
                   "FF_30min" = "darkorange",
                   "Cells" = "darkslategray3" )

p4 <- ggplot(c3, aes(x =pos, y = mito, group =Sample)) +
  geom_line(aes(color = Sample)) +
  scale_y_continuous(breaks = c(5,10,15,20,25)) +
  scale_color_manual(values=dataset_color) +
  theme_classic() +
  theme(legend.position="bottom",
        axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=20),
        axis.line.x = element_line(),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=20))

p4

ggsave(plot = p4, filename = "mean_cov_mtDNA_subset_high_conf.pdf", width = 10, height = 4, dpi=300)



