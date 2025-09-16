library(dplyr)
library(SummarizedExperiment)
library(ggplot2)
library(tidyr)
# Set input and output directories
output_dir <- "/.../output_plots"
# 1. Set input directory first
input_dir <- "/.../output_rds/"

source(".../01_read_in_variants.R")

# Example dataframe: replace with your actual df
df <- summary_df  # assuming this is your full table with all columns

# If needed: convert Age to numeric
df$Age_numeric <- as.numeric(gsub(",", ".", gsub(" years", "", df$Age)))

df

# Optional: check data types
str(df)
cor_test <- cor.test(df$Age_numeric, df$n_unique_variants_union, method = "pearson")
print(cor_test)


ggplot(df, aes(x = Age_numeric, y = n_unique_variants_union)) +
  geom_point(aes(color = project, shape = Patient), size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Mean mtDNA depth per cell", y = "Number of detected variants") +
  theme_minimal()


plot_df <- df %>%
  dplyr::select(Age_numeric, Patient, n_rows, n_unique_variants_union) %>%
  pivot_longer(cols = c(n_rows, n_unique_variants_union), 
               names_to = "metric", values_to = "value")

ggplot(plot_df, aes(x =Age_numeric, y = value, color = metric, shape = metric)) +
  geom_point(size = 4) +
  labs(x = "Sample", y = "Number of Variants") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Make sure Age is numeric
df$Age_numeric <- as.numeric(gsub(",", ".", gsub(" years", "", df$Age)))

unique(df$file)

sample.color <- c("ALS" = "#8C918D", "GBM" = "#424543", "pediatric neoplasm" = "#C5CDC6")

color.file <- c("GN01" = "#B37BB3",
                "NB01" = "#76C3A9",
                "NB02" = "#355C4C",      
                "ALS_1MNC" = "orange2",    
                "ALS_MO" = "dodgerblue3",     
                "ALS_CSC" = "plum",     
                "ALS_PFC" = "firebrick",      
                "Patient1_primary" ="#5B79BE", 
                "Patient1_recurrence" = "#8C233F",
                "Patient2_primary" ="#A7BAF2", 
                "Patient2_recurrence" ="#c1121f")


p <- ggplot(df, aes(x = Age_numeric, y = n_rows)) +
  geom_point(aes(fill = sample, shape = Patient, colour = project), 
             size = 5, 
             stroke = 0.4) +     # optional: line thickness
  scale_fill_manual(values = color.file) +
  scale_color_manual(values = sample.color) +
  scale_shape_manual(values = 21:25) +  # optional: define shape range for up to 5 Patients
  labs(x = "Age", y = "# mtDNA variants") +
  theme_minimal(base_size = 14)

p
ggsave(file.path(output_dir, "variants_per_age.pdf"), plot = p, width = 10, height = 8, dpi = 300)

