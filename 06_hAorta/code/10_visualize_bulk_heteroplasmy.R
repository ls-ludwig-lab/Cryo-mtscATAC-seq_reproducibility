library(data.table)
library(dplyr)
library(stringr)
library(Signac)
library(SummarizedExperiment)
library(ggplot2)
library(ComplexHeatmap)
library(tibble)
library(tidyr)
library(viridis)

output_dir <- ".../output/"


#variants called in this working environment

primary <- readRDS(file=".../output/variants.rds")
primary

df <- readRDS(file="preprocessed.rds")
primary <- primary[,colnames(df)]

misc_df_primary <- data.frame(rowData(primary))
vars_primary <- misc_df_primary %>%  dplyr::filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 5) %>% pull(variant)


# Export the variants called in either culture
length(vars_primary)

prim <- primary[vars_primary,]

# Annotate where the variants came from 
rowData(prim)$called_in_primary <- rownames(prim) %in% vars_primary

# Make plot of allele frequencies

mean_df <- data.frame(
  primary = rowData(prim)$mean,
  v2 = rownames(prim))

mean_df$pct <- mean_df$primary*100

mean_df$v2

mean_df
setwd(".../output/")

# Assuming `mean_df` is already sorted and mutated as you did
mean_df <- mean_df %>%
  arrange(desc(pct)) %>%
  mutate(v2 = factor(v2, levels = v2))

# Get the x-axis position of the 51st variant
vline_pos <- 51

# Plot
p1 <- ggplot(mean_df, aes(x = v2, y = pct)) +
  geom_point(aes(color = pct), size = 2) +
  geom_vline(xintercept = vline_pos, linetype = "dashed", color = "black") +
  scale_color_viridis(name = "Heteroplasmy", option = "D", direction = -1) +
  labs(x = "mtDNA Variant", y = "Mean Heteroplasmy", title = "Pseudobulk Heteroplasmy Events") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

p1

# Save the plot
output_file <- file.path(output_dir, "pseudobulk_heteroplasmy_events.pdf")
ggsave(filename = output_file, plot = p1, width = 5, height = 5)


