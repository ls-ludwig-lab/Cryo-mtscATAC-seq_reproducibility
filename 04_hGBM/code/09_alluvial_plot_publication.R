## ---- Setup ----
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggalluvial)  # make sure installed

setwd(".../output")

## ---- Load ----
df <- readRDS(file="integrated_annotated.rds")
meta.data <- df@meta.data
meta.data
colnames(df@meta.data)

write.csv(meta.data, "meta.data.csv")

df <- read.csv("meta.data.csv", check.names = FALSE)

## Optional sanity checks
# unique(df$Celltype_simplified)
# unique(df$Patient)
# unique(df$Donor)

## ---- Subset to one patient ----
df1 <- subset(df, Patient == "Donor1")

## ---- Proportions per Donor ----
# counts cell types within each Donor
ab <- table(df1$Celltype_simplified, df1$Donor)
ab
unique(df$Celltype_simplified)

# convert to proportions (column-wise)
prop <- sweep(x = ab, MARGIN = 2, STATS = colSums(ab), FUN = "/")

## Build long df from your 'prop' table
df.plot <- as.data.frame(prop) |>
  dplyr::rename(CellType = Var1, Sample = Var2, CellProp = Freq) |>
  # collapse any accidental duplicates
  dplyr::group_by(Sample, CellType) |>
  dplyr::summarise(CellProp = sum(CellProp), .groups = "drop")
  # drop zero rows (not needed by ggalluvial and avoids weirdness)

## Order + colors (keep your Microglia 1 spelling)
celltype_order <- c(
  "malignant_other","NPC1-OPC-like","AC-MEs-like","OPC-like","OPC-like 2",
  "AC-like","Oligodendrocytes","Neuronal Cells","DLX-GAD-high Neuronal Cells",
  "Endothelial Cells","Pericytes", "Microglia","Microglia 1","T Cells"
)

unique(df1$Celltype_simplified)

Celltype.color <- c(
  "malignant_other" = "#829BD4",
  "NPC1-OPC-like" = "lightslateblue",
  "AC-MEs-like" = "#A7BAF2",
  "OPC-like" = "dodgerblue",
  "OPC-like 2" = "dodgerblue3",
  "AC-like" = "paleturquoise3",
  "Oligodendrocytes" = "dodgerblue4",
  "Neuronal Cells" = "orange",
  "DLX-GAD-high Neuronal Cells" = "sienna2",
  "Endothelial Cells" = "#400135",
  "Pericytes" = "#73026B",
  "Microglia" = "#D9488B",
  "Microglia 1" = "#F294AD",
  "T Cells" = "#8C233F"
)

df.plot <- df.plot |>
  mutate(
    CellType = factor(CellType, levels = celltype_order),
    Sample = factor(Sample)
  )

## (Optional) sanity check for duplicates
# df.plot %>% count(Sample, CellType) %>% filter(n > 1)

## Plot (switch to geom_alluvium)
p1 <- ggplot(df.plot,
             aes(x = Sample, stratum = CellType, alluvium = CellType,
                 y = CellProp, fill = CellType)) +
  geom_alluvium(width = 0.6) +
  geom_stratum(width = 0.6, color = "grey40") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = Celltype.color, drop = FALSE) +
  labs(y = "Cell Proportion (%)", x = NULL, fill = "Cell Type") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

p1
ggsave(plot = p1, filename = "PropCells_Allu_Donor1.pdf",
       width = 3, height = 6, dpi = 300)



## ---- Subset to one patient ----
df1 <- subset(df, Patient == "Donor2")

## ---- Proportions per Donor ----
# counts cell types within each Donor
ab <- table(df1$Celltype_simplified, df1$Donor)
ab
unique(df$Celltype_simplified)

# convert to proportions (column-wise)
prop <- sweep(x = ab, MARGIN = 2, STATS = colSums(ab), FUN = "/")

## Build long df from your 'prop' table
df.plot <- as.data.frame(prop) |>
  dplyr::rename(CellType = Var1, Sample = Var2, CellProp = Freq) |>
  # collapse any accidental duplicates
  dplyr::group_by(Sample, CellType) |>
  dplyr::summarise(CellProp = sum(CellProp), .groups = "drop")
# drop zero rows (not needed by ggalluvial and avoids weirdness)

## Order + colors (keep your Microglia 1 spelling)
celltype_order <- c(
  "malignant_other","NPC1-OPC-like","AC-MEs-like","OPC-like","OPC-like 2",
  "AC-like","Oligodendrocytes","Neuronal Cells","DLX-GAD-high Neuronal Cells",
  "Endothelial Cells","Pericytes", "Microglia","Microglia 1","T Cells"
)

unique(df1$Celltype_simplified)

Celltype.color <- c(
  "malignant_other" = "#829BD4",
  "NPC1-OPC-like" = "lightslateblue",
  "AC-MEs-like" = "#A7BAF2",
  "OPC-like" = "dodgerblue",
  "OPC-like 2" = "dodgerblue3",
  "AC-like" = "paleturquoise3",
  "Oligodendrocytes" = "dodgerblue4",
  "Neuronal Cells" = "orange",
  "DLX-GAD-high Neuronal Cells" = "sienna2",
  "Endothelial Cells" = "#400135",
  "Pericytes" = "#73026B",
  "Microglia" = "#D9488B",
  "Microglia 1" = "#F294AD",
  "T Cells" = "#8C233F"
)

df.plot <- df.plot |>
  mutate(
    CellType = factor(CellType, levels = celltype_order),
    Sample = factor(Sample)
  )

## (Optional) sanity check for duplicates
# df.plot %>% count(Sample, CellType) %>% filter(n > 1)

## Plot (switch to geom_alluvium)
p1 <- ggplot(df.plot,
             aes(x = Sample, stratum = CellType, alluvium = CellType,
                 y = CellProp, fill = CellType)) +
  geom_alluvium(width = 0.6) +
  geom_stratum(width = 0.6, color = "grey40") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = Celltype.color, drop = FALSE) +
  labs(y = "Cell Proportion (%)", x = NULL, fill = "Cell Type") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

p1
ggsave(plot = p1, filename = "PropCells_Allu_Donor2.pdf",
       width = 3, height = 6, dpi = 300)
