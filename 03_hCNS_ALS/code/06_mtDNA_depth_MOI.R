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



setwd(".../Region/output/")
d <- readRDS(file="integrated_annotated.rds")

# Define the regions you want to subset
regions <- c("Region_5", "Region_10", "Region_28", "Region_29")  # Add more as needed

#Now you have df_Region_5, df_Region_10, etc. as separate variables
for (region in regions) {
  # Subset the dataset
  df <- subset(d, subset = dataset == region)
  
  # Modify the cell names
  new_cell_names <- gsub(pattern = paste0("^", region, "_"), replacement = "", x = Cells(df))
  
  # Rename the cells
  df <- RenameCells(object = df, new.names = new_cell_names)
  
  # Assign the modified Seurat object to a variable with the region name
  assign(paste0("df_", gsub("[^a-zA-Z0-9]", "_", region)), df)
  
  # Print first few cell names to check
  print(head(Cells(df)))
}



# Access modified Seurat objects using subset_list[["Region_5"]], etc.


color.palette.Region <-c("firebrick", "orange2","dodgerblue3","plum")
df$Region <- factor(df$Region, levels = c("Region_5", "Region_10", "Region_28", "Region_29"))



#coverage mtDNA chromsome

#load for this the mgatk output rds file

# -- Functions -- #
pull_coverage <- function(SE, cells, resolution = 5){
  zoo::rollmean(rowMeans(assays(SE)[["coverage"]]), resolution)
}

mito.SE <- readRDS(".../Region_5_mtMask.rds")

mito.SE <- mito.SE[,colnames(df_Region_5)]
# --- 1. Visualize mitochondrial coverage --- #
cov_df <- data.frame(
  pos = zoo::rollmean(1:16569, 5),
  mito = pull_coverage(mito.SE))

cov_df$Region <- c("Region_5")


mito.SE1 <- readRDS(".../Region_10_mtMask.rds")

mito.SE1 <- mito.SE1[,colnames(df_Region_10)]
cov_df1 <- data.frame(
  pos = zoo::rollmean(1:16569, 5),
  mito = pull_coverage(mito.SE1))

cov_df1$Region <- c("Region_10")


mito.SE2 <- readRDS(".../Region_28_mtMask.rds")

mito.SE2 <- mito.SE2[,colnames(df_Region_28)]
cov_df2 <- data.frame(
  pos = zoo::rollmean(1:16569,5),
  mito = pull_coverage(mito.SE2))

cov_df2$Region <- c("Region_28")


mito.SE3 <- readRDS(".../Region_29_mtMask.rds")
mito.SE3 <- mito.SE3[,colnames(df_Region_29)]
cov_df3 <- data.frame(
  pos = zoo::rollmean(1:16569,5),
  mito = pull_coverage(mito.SE3))

cov_df3$Region <- c("Region_29")

c3 <- rbind(cov_df,cov_df1, cov_df2, cov_df3)

c3$Region <- factor(c3$Region, levels = c("Region_5", "Region_10", "Region_28", "Region_29"))


p4 <- ggplot(c3, aes(x =pos, y = mito, group =Region)) +
  geom_line(aes(color = Region)) +
  scale_y_continuous(breaks = c(5,10,20,30,40,50,60)) +
  scale_color_manual(values=color.palette.Region)+ 
  theme_classic() +
  theme(legend.position="bottom",
        axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=20),
        axis.line.x = element_line(),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=20))

p4

ggsave(plot = p4, filename = "mean_cov_mtDNA_subset_high_conf.pdf", width = 10, height = 4, dpi=300)


