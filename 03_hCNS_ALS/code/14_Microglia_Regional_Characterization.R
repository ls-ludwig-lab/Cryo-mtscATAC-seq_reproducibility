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
library(ggrepel)


#### Microglia subset analysis ####
#############################

setwd(".../test_output")

df <- readRDS(file="Microglia_subset.rds")

DefaultAssay(df) <- 'ATAC'
gene.activities <- GeneActivity(df)

# add to the Seurat object as a new assay
df[['RNA']] <- CreateAssayObject(counts = gene.activities)

df <- NormalizeData(
  object = df,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(df$nCount_RNA)
)

DefaultAssay(df) <- "ATAC"

comparisons <- list(
  "Region_10_vs_Region_5" = c("Region_5", "Region_10"),
  "Region_28_vs_Region_5" = c("Region_5", "Region_28"),
  "Region_29_vs_Region_5" = c("Region_5", "Region_29")
)


Idents(df) <- "Region"

results <- lapply(comparisons, function(regions) {
  FindMarkers(
    object = df,
    ident.1 = regions[1],
    ident.2 = regions[2],
    test.use = "wilcox",
    min.pct = 0.1
  )
})

results

da_peaks <-results[["Region_29_vs_Region_5"]]

ggplot(data=da_peaks, aes(x=avg_log2FC, y=-log10(p_val))) + geom_point() 

# 1. Add new column to data frame
da_peaks$diffexpressed <- "Not significant"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as up-regulated
da_peaks$diffexpressed[da_peaks$avg_log2FC > 0.6 & da_peaks$p_val_adj < 0.05] <- "Up-regulated in Region 5"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as down-regulated
da_peaks$diffexpressed[da_peaks$avg_log2FC < -0.6 & da_peaks$p_val_adj < 0.05] <- "Down-regulated in Region 5"


# 3. Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
da_peaks$delabel <- NA
da_peaks$delabel[da_peaks$diffexpressed != "Not significant"] <- rownames(da_peaks)[da_peaks$diffexpressed != "Not significant"]

# 4. Add to plot
ggplot(data=da_peaks, aes(x=avg_log2FC, y=-log10(p_val_adj), label = delabel)) +
  geom_point(aes(col=diffexpressed)) +
  geom_text() +
  xlab(expression("average log"[2]*"FC")) +
  ylab(expression("-log"[10]*"FDR")) +
  labs(color = "Differentially expressed genes") +
  scale_color_manual(values = c("dodgerblue3","gray50", "firebrick3"))

# 5. Use repel
p1 <- ggplot(data=da_peaks, aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(aes(col=diffexpressed)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
        panel.background = element_rect(fill ="white",
                                        color = "white"),
        legend.position = "none") +
  xlab(expression("average log"[2]*"FC")) +
  ylab(expression("-log"[10]*"FDR")) +
  labs(color = "Differentially expressed genes") +
  scale_color_manual(values = c("dodgerblue3","gray50","firebrick3"))+
  geom_vline(xintercept=c(-0.6, 0.6), linetype='dotted', col="black") +
  geom_hline(yintercept=-log10(0.05), linetype='dotted', col="black")

p1 

ggsave(plot = p1, filename = "DA_peaks_Clone_Microglia_Region5vs29.pdf", width = 5, height = 5, dpi = 300)


open_Region5 <- rownames(da_peaks[da_peaks$avg_log2FC > 1.5, ])
open_Region28<- rownames(da_peaks[da_peaks$avg_log2FC < -1.5, ])

closest_Region5 <- ClosestFeature(df, regions = open_Region5)
closest_Region5
write.csv(closest_Region5, "closest_Region5.csv")
open_Region28 <- ClosestFeature(df, regions = open_Region28)
open_Region28
write.csv(open_Region28, "open_Region29.csv")
