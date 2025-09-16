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



setwd(".../output/")

d <- readRDS(file="integrated_annotated.rds")

unique(d$dataset)

# Define the regions you want to subset
datasets <- c("patient1_primary", "patient2_recu", "patient1_recu", "patient2_primary")  
unique_idents <- c("patient1_primary_single_", "patient2_recu_single_", "patient1_recu_single_", "patient2_primary_single_")  

#Looping with seq_along(datasets)
#Ensures we correctly map each dataset to its corresponding unique_idents value.

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

mito.SE <- readRDS("mgatk_output.rds")

mito.SE <- mito.SE[,colnames(df_patient1_primary)]
# --- 1. Visualize mitochondrial coverage --- #
cov_df <- data.frame(
  pos = zoo::rollmean(1:16569, 5),
  mito = pull_coverage(mito.SE))

cov_df$Sample <- c("patient1_primary")


mito.SE1 <- readRDS("/mgatk_output.rds")
mito.SE1 <- mito.SE1[,colnames(df_patient1_recu)]
cov_df1 <- data.frame(
  pos = zoo::rollmean(1:16569, 5),
  mito = pull_coverage(mito.SE1))

cov_df1$Sample <- c("patient1_recu")


mito.SE2 <- readRDS("mgatk_output.rds")
mito.SE2 <- mito.SE2[,colnames(df_patient2_primary)]

cov_df2 <- data.frame(
  pos = zoo::rollmean(1:16569,5),
  mito = pull_coverage(mito.SE2))

cov_df2$Sample <- c("patient2_primary")


mito.SE3 <- readRDS("/mgatk_output.rds")

mito.SE3 <- mito.SE3[,colnames(df_patient2_recu)]

cov_df3 <- data.frame(
  pos = zoo::rollmean(1:16569,5),
  mito = pull_coverage(mito.SE3))

cov_df3$Sample <- c("patient2_recu")

c3 <- rbind(cov_df,cov_df1, cov_df2, cov_df3)

c3$Sample <- factor(c3$Sample, levels = c("patient1_primary", "patient1_recu", "patient2_primary", "patient2_recu"))

color.palette.Donor <-c("steelblue1", "hotpink2",
                        "steelblue3","hotpink4")

p4 <- ggplot(c3, aes(x =pos, y = mito, group =Sample)) +
  geom_line(aes(color = Sample)) +
  scale_y_continuous(breaks = c(5,10,20,30,40,50,60)) +
  scale_color_manual(values=dataset.color)+ 
  theme_classic() +
  theme(legend.position="bottom",
        axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=20),
        axis.line.x = element_line(),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=20))

p4

ggsave(plot = p4, filename = "mean_cov_mtDNA_subset_high_conf.pdf", width = 10, height = 4, dpi=300)


mean(cov_df$mito)
#[1] 23.04134
 mean(cov_df1$mito)
#[1] 14.04495
mean(cov_df2$mito)
#[1] 22.40815
mean(cov_df3$mito)
#[1] 19.77207
