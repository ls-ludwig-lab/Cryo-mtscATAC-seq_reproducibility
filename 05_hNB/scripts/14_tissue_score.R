# load packages
library(ggplot2)
library(tidyverse)
library(dplyr)

set.seed(257)

# Load Tissue Score Data 
Tissue_Score <- read.table("Tissue_Score.tsv")


# ============================================
# Calculate Tissue Score (TS)
# ============================================
# Formula: TS = average of 
# - Intact nuclei (I)
# - Adjusted Necrosis: 1 - Necrosis_N (higher necrosis reduces score)
# - Adjusted Freezing Artefacts: 1 - Freezing_Artefacts_FA (higher artefacts reduce score)

Tissue_Score$Tissue_Score_TS <- round(
  (Tissue_Score$Intact_Nuclei_I + (1 - Tissue_Score$Necrosis_N) + (1 - Tissue_Score$Freezing_Artefacts_FA)) / 3,
  7
)

# ============================================
# Define Color Palette
# ============================================
color_palette <- c(
  "NB01" = "#95BBA1", "NB02" = "#658381", "GN01" = "#b37bb3",
  "NB03" = "#F1656E", "NB04" = "#934258", "NB05" = "#A8D7E3",
  "NB06" = "#616887", "NB07" = "#F9C071", "NB08" = "#E6A476",
  "NB09" = "#F8C0BE", "NB10" = "#FF6699", "NB11" =  "#CCCCFF", "NB12" = "#333399"
)

# ============================================
# Correlation Analysis: Tissue Score vs Fragments_Reads_Peaks
# ============================================

# Pearson correlation
cor_test_fragments <- cor.test(Tissue_Score$Tissue_Score_TS, Tissue_Score$Fragments_Reads_Peaks, method = "pearson")
print(cor_test_fragments)

r_value <- round(cor_test_fragments$estimate, 2)
p_value <- signif(cor_test_fragments$p.value, 3)

# Plot with regression line and correlation annotation
ggplot(Tissue_Score, aes(x = Tissue_Score_TS, y = Fragments_Reads_Peaks, color = Identifier)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.5, fill = "grey85") +
  geom_point(size = 5) +
  scale_color_manual(values = color_palette) +
  # geom_text(aes(label = Pseudonym), size = 4, vjust = -1, show.legend = FALSE) +  # Optional: label points
  labs(x = "Tissue Score (TS)", y = "Fragments Reads Peaks") +
  annotate("text", x = 0.1, y = max(Tissue_Score$Fragments_Reads_Peaks, na.rm = TRUE), 
           label = paste0("r = ", r_value, ", p = ", p_value), color = "red", size = 6, hjust = 0) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )

# ============================================
# Correlation Analysis: Tissue Score vs TSS Enrichment
# ============================================

# Pearson correlation
cor_test_tss <- cor.test(Tissue_Score$Tissue_Score_TS, Tissue_Score$TSS, method = "pearson")
print(cor_test_tss)

r_value <- round(cor_test_tss$estimate, 2)
p_value <- signif(cor_test_tss$p.value, 3)

# Plot with regression line and correlation annotation
ggplot(Tissue_Score, aes(x = Tissue_Score_TS, y = TSS, color = Identifier)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.5, fill = "grey85") +
  geom_point(size = 5) +
  scale_color_manual(values = color_palette) +
  labs(x = "Tissue Score (TS)", y = "TSS") +
  ylim(0, 5) +
  annotate("text", x = 0.1, y = 5, label = paste0("r = ", r_value, ", p = ", p_value),
           color = "red", size = 6, hjust = 0) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )