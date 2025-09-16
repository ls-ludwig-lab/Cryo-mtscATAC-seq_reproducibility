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
ref_all <- fread(".../data/final/chrM_refAllele.txt")
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
setwd(".../output/")

called_variants <- rowData(readRDS("variants_ncof1.rds"))$variant


ref_all_long$called <- ref_all_long$variant %in% called_variants

# Compute changes in expected/observed
total <- dim(ref_all_long)[1]
total_called <- sum(ref_all_long$called)
prop_df <- ref_all_long %>% group_by(three_plot, group_change, strand) %>%
  summarize(observed_prop_called = sum(called)/total_called, expected_prop = n()/total, n = n()) %>%
  mutate(fc_called = observed_prop_called/expected_prop)

prop_df$change_plot <- paste0(prop_df$group_change, "_", prop_df$three_plot)

# Visualize
p1 <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_called)) +
  geom_bar(stat = "identity", position = "dodge")  +
  theme(axis.title.x=element_blank(),
        axis.text.x =element_text(angle=45, h =1),
        panel.background = element_blank(),
        axis.line = element_line(color = "grey", size = 1))+
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide", y = "Substitution Rate (Expected / Observed)")
p1
cowplot::ggsave2(p1, file = "all_mito_signature.pdf", width = 10, height = 5)


