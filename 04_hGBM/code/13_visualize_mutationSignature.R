library(data.table)
library(dplyr)
library(BuenColors)
library(SummarizedExperiment)
library(ggplot2)
# Simple reverse complement function
reverse_complement <- function(s){
  chartr("ATGC","TACG",s)
}

# Process 3 digit signature based on letters
ref_all <- fread(".../chrM_refAllele.txt")
colnames(ref_all) <- c("pos", "ref")
ref_all$ref <- toupper(ref_all$ref)
l <- as.character(ref_all$ref)

# Gs happen to be at the first and last position
ref_all$three <- paste0(c("G", l[-length(l)]), l, c(l[-1], "G"))

# Remove Ns
ref_all <- ref_all[!grepl("N", ref_all$three),]

# Make every possible mutation
ref_all_long <- rbind(ref_all,ref_all, ref_all,ref_all)
ref_all_long$alt <- rep(c("A", "C", "G", "T"), each = dim(ref_all)[1])
ref_all_long <- ref_all_long[ref_all_long$ref != ref_all_long$alt,]

# add some meta data
ref_all_long$variant <- paste0(as.character(ref_all_long$pos), ref_all_long$ref, ">", ref_all_long$alt)
ref_all_long$change <- paste0(ref_all_long$ref, ref_all_long$alt)
ref_all_long$change_rc <- reverse_complement(paste0(ref_all_long$ref, ref_all_long$alt))

# A/G rich strand is "heavy" -- https://en.wikipedia.org/wiki/Heavy_strand
table(ref_all$ref) # so the reference strand is light (more C/T)
ref_all_long$strand <- ifelse(ref_all_long$ref %in% c("C","T"), "L", "H")

# Change to C/T as ref allele
ref_all_long$rc3 <- reverse_complement(ref_all_long$three)
ref_all_long$three_plot <- ifelse(ref_all_long$strand == "L", ref_all_long$three, ref_all_long$rc3)
ref_all_long$group_change <- ifelse(ref_all_long$strand == "L", ref_all_long$change, ref_all_long$change_rc)

# Annotate with called variants
called_variants1 <- rowData(readRDS("ncof1_variants.rds"))$variant
called_variants2 <- rowData(readRDS("ncof1_variants.rds"))$variant
called_variants3 <- rowData(readRDS("ncof1_variants.rds"))$variant
called_variants4 <-rowData(readRDS("ncof1_variants.rds"))$variant


#Looping/ Lapply through all different Regions -> are there differences of the profile in each Region?
#Explanation of Lapply
  
#Input Preparation: A list (called_variants_list) is created to store all called_variants objects.
#lapply Loop: Each iteration processes one called_variants object:
#  Updates ref_all_long$called based on the current called_variants.
#Computes the proportions and fold change (prop_df).
#Generates and saves a plot for the current iteration.
#Saving Outputs: Plots are saved with unique file names based on their iteration index.
#Output Handling: Results are stored in a list containing both plots and dataframes for further use.

# List of called_variants objects
called_variants_list <- list(called_variants1, called_variants2, called_variants3, called_variants4)


# Perform the process for each called_variants and save the plots
results <- lapply(seq_along(called_variants_list), function(i) {
  ref_all_long$called <- ref_all_long$variant %in% called_variants_list[[i]]
  
  # Compute changes in expected/observed
  total <- dim(ref_all_long)[1]
  total_called <- sum(ref_all_long$called)
  prop_df <- ref_all_long %>% 
    dplyr::group_by(three_plot, group_change, strand) %>%
    dplyr::summarise(
      n = dplyr::n(),
      observed_prop_called = sum(called) / total_called, 
      expected_prop = n / total,
      .groups = "drop"
    ) %>%
    dplyr::mutate(fc_called = observed_prop_called / expected_prop)
  
  
  prop_df$change_plot <- paste0(prop_df$group_change, "_", prop_df$three_plot)
  
  # Combine all variants across all datasets first to find the global max
  combined_prop_df <- do.call(rbind, lapply(seq_along(called_variants_list), function(i) {
    ref_all_long$called <- ref_all_long$variant %in% called_variants_list[[i]]
    total <- nrow(ref_all_long)
    total_called <- sum(ref_all_long$called)
    ref_all_long %>% 
      group_by(three_plot, group_change, strand) %>%
      summarise(
        observed_prop_called = sum(called) / total_called, 
        expected_prop = n() / total,
        .groups = "drop"
      ) %>%
      mutate(
        fc_called = observed_prop_called / expected_prop,
        change_plot = paste0(group_change, "_", three_plot)
      )
  }))
  
  YMAX <- ceiling(max(combined_prop_df$fc_called, na.rm = TRUE))
  
  # Visualize
  p <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_called)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "grey", size = 1)) +
    scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 1, linetype = 2, color = "black") +
    coord_cartesian(ylim = c(0, YMAX)) +
    labs(x = "Change in nucleotide", y = "Substitution Rate (Expected / Observed)")
  
  
  # Display the plot in RStudio
  print(p)
  
  # Save the plot
  output_file <- paste0(".../mito_signature_", i, ".pdf")
  cowplot::ggsave2(p, file = output_file, width = 6, height = 5)
  
  return(list(plot = p, prop_df = prop_df))
})




