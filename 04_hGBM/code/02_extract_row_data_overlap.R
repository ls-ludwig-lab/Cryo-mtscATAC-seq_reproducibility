library(data.table)
library(dplyr)
library(stringr)
library(Signac)
library(SummarizedExperiment)
library(ggplot2)
library(ComplexHeatmap)
library(tibble)
library(tidyr)
library(circlize)


# The workflow shown here is representative for all Datasets


#variants called in this working environment

primary <- readRDS(file=".../output/primary.rds")
primary

relapse <- readRDS(file=".../output/primary.rds")
relapse

misc_df_primary <- data.frame(rowData(primary))
vars_primary <- misc_df_primary %>%  filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 5) %>% pull(variant)


misc_df_relapse <- data.frame(rowData(relapse))
vars_relapse <- misc_df_relapse %>%  filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 & mean_coverage >= 5) %>% pull(variant)

# Export the variants called in either culture
length(vars_primary)
length(vars_relapse)
vars_both <- unique(c(vars_primary, vars_relapse))
length(vars_both)

prim <- primary[vars_both,]
rel <- relapse[vars_both,]

# Annotate where the variants came from 
rowData(prim)$called_in_primary <- rownames(prim) %in% vars_primary
rowData(rel)$called_in_relapse <- rownames(rel) %in% vars_relapse

# Make plot of allele frequencies

mean_df <- data.frame(
  primary = rowData(prim)$mean,
  relapse = rowData(rel)$mean, 
  variant = rownames(prim),
  v2 = rownames(rel)
) %>% mutate(log10_FC = log10(relapse/primary))  %>% arrange(desc(log10_FC))

cor(log10(mean_df$relapse), log10(mean_df$primary))
cor((mean_df$relapse), (mean_df$primary))

order <- mean_df %>% 
  arrange(desc(log10_FC))
setwd("/.../called_variants")


label_data <- subset(mean_df, variant %in% c("12865A>G","12008G>A"))
p1 <- ggplot() +
  geom_point(data = mean_df, aes(x = primary, y = relapse, color = "grey"), size = 2) +
  geom_point(data = subset(mean_df, variant == "12865A>G"), aes(x = primary, y = relapse, color = "12865A>G"), size = 2) +
  geom_point(data = subset(mean_df, variant == "12008G>A"), aes(x = primary, y = relapse, color = "12008G>A"), size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  # Add labels with matching colors
  geom_text(data = label_data, aes(x = primary, y = relapse, label = variant, color = variant), 
            hjust = -0.1, vjust = -0.5, size = 3) +  # Adjust `hjust` and `vjust` for positioning
  scale_color_manual(values = c("12865A>G" = "#8C233F", 
                                 "12008G>A" = "#5B79BE", 
                                "grey" = "grey")) +
  labs(x = "Primary Heteroplasmy %", y = "Relapse Heteroplasmy %") + 
  xlim(0,0.01) +
  ylim(0,0.01) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size = 14),  # Adjust the size of x-axis text
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),  # Adjust the size of x-axis label
    axis.title.y = element_text(size = 16),
    legend.position = "none",
    #plot.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  )
p1

cowplot::ggsave2(p1, file = "called_variants_ncof5.pdf", width = 5, height = 5)




