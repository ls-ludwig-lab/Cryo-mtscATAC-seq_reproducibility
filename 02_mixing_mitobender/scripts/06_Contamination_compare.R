library(dplyr)
library(data.table)
library(stringr)
library(BuenColors)
library(ggrastr)

#Donor specific variants identified in Heatmap - are not displayed here due to human genotyping
#Donor specific variants are homoplasmic, haplotype specific variants, furthermore homoplasmic variants could also be germline variants which are family or person specific and therefore due to data protection of patients not displayed


donor0_vars <- c("variant_1", "variant_2", "variant_3","...")
donor1_vars <- c("variant_4", "variant_5", "variant_6","...")


# Function to parse out the parts of the string as necessary
substrRight <- function(x, n =1){
  substr(x, nchar(x)-n+1, nchar(x))
}

substrLeft <- function(x, n = 1){
  substr(x, 1, nchar(x)-n)
}

# Pull the letter for each cell line
# For donor0 variants, its pos_donor1>donor0
donor0temp <- str_split_fixed(donor0_vars, ">", 2)
donor0var_df <- data.frame(
  pos = as.numeric(substrLeft(donor0temp[,1])),
  donor1 = substrRight(donor0temp[,1]),
  donor0 = donor0temp[,2]
)

# For donor1 variants, its pos_donor0>donor1
donor1temp <- str_split_fixed(donor1_vars, ">", 2)

donor1var_df <- data.frame(
  pos = as.numeric(substrLeft(donor1temp[,1])),
  donor0 = substrRight(donor1temp[,1]),
  donor1 = donor1temp[,2]
)

# Define a global variable of all of the variants 
df2 <- rbind(donor1var_df, donor0var_df)

# Function to get the per barcode counts per cell type (donor1/donor0) by alt letter
# Counts base-specific reads for donor0/donor1 variants using extractme()
extractme <- function(cell, letter, SE){
  idxx <- df2[which(df2[,cell] == letter), "pos"]
  Matrix::colSums(assays(SE)[[paste0(letter, "_counts_fw")]][idxx, ] + assays(SE)[[paste0(letter, "_counts_rev")]][idxx, ])
}

# Given a summarized experiment from mgatk, compute the essentials for ultimately determining contamination
process_SE_contamination <- function(SE, library){
  df3 <- data.frame(
    #Barcode-wise summary
    barcode = colnames(SE),
    #Total base counts per barcode, Sum of extractme("donor0", base, SE) for A, C, G, T bases.
    donor0 = extractme("donor0", "A", SE) +  extractme("donor0", "C", SE) +  extractme("donor0", "G", SE) +  extractme("donor0", "T", SE),
    donor1 = extractme("donor1", "A", SE) +  extractme("donor1", "C", SE) +  extractme("donor1", "G", SE) +  extractme("donor1", "T", SE), 
    #depth column is taken from colData(SE)$depth.
    depth = colData(SE)$depth
  )
  #Calculates the fraction of the less abundant population for each cell
  # minor fraction (as %) of the total signal per cell
  #If a cell’s minor fraction is very low → likely pure assignment to one donor.
  #If a cell’s minor fraction is high → suggests possible contamination, doublets, or mixed populations.
  
  #or in other word:The percentage of reads coming from the less abundant mitochondrial donor in this cell (either donor0 or donor1).

  df3 <-  df3 %>% mutate(minor_population = pmin(donor0/(donor1 + donor0 + 0.001)*100 ,donor1/(donor1 + donor0 + 0.001)*100), library = library)
  df3                        
}

# Function to estimate contamination based on reads
estimate_contamination <- function(df3, rm_doublets = FALSE){
  if(rm_doublets){
    #If rm_doublets=TRUE, it removes any barcodes where the minor population fraction is ≥5% (a typical cutoff for "likely doublets").
    df3 <- df3[df3$minor_population < 5,]
  }
  #Calculates the fraction of total reads that are “contaminant” across all barcodes
  sum(pmin(df3$donor0, df3$donor1)) / sum(df3$donor0 + df3$donor1)*100
}



#Plotting:
#X-axis: number of mtDNA reads (A/C/G/T combined) that match donor0-specific alleles in this cell (each dot is a cell)
#Y-axis: the same but for donor1 specifi
#Color:(minor population): % of reads that cell that came from the less abundant donor

#Bottom-right	Likely donor0-only cells (high donor0 signal, low donor1)
#Top-left	Likely donor1-only cells (high donor1 signal, low donor0)
#Diagonal line	Similar counts for donor0 and donor1 → likely doublets or ambient mtDNA
#Color intensity	Warm colors (red/orange) = higher contamination (minor_population closer to 50%)

#estimate_contamination(df3, rm_doublets=FALSE): contamination in raw data.
#estimate_contamination(df3, rm_doublets=TRUE): contamination after removing barcodes likely to be doublets.

#output mgatk
fix0p5pct1 <- process_SE_contamination(readRDS(".../mgatk.rds") ,"unfiltered")
fix0p5pct1

head(fix0p5pct1)
barcodes <- fix0p5pct1$barcode

#output mitobender
fix0p5pct <- readRDS(".../mitoBender.rds")
fix0p5pct <- fix0p5pct[, colnames(fix0p5pct) %in% barcodes]
fix0p5pct <- process_SE_contamination(fix0p5pct, "filtered")
head(fix0p5pct)

get_density <- function(x, y, n = 100) {
  #Uses kernel density estimation to get a smoothed 2D density grid.
  dens <- MASS::kde2d(x, y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  density <- dens$z[cbind(ix, iy)]
  return(density)
}

fix0p5pct$density <- get_density(log10(fix0p5pct$donor0 + 1), log10(fix0p5pct$donor1 + 1))
head(fix0p5pct)
fix0p5pct1$density <- get_density(log10(fix0p5pct1$donor0 + 1), log10(fix0p5pct1$donor1 + 1))
head(fix0p5pct1)

estimate_contamination(fix0p5pct, FALSE)
estimate_contamination(fix0p5pct1, FALSE)

estimate_contamination(fix0p5pct, TRUE)
estimate_contamination(fix0p5pct1, TRUE)

# Create a named vector with the contamination rates
contam_rates <- c(
  fix0p5pct1_noDoublets = estimate_contamination(fix0p5pct1, FALSE),
  fix0p5pct_noDoublets = estimate_contamination(fix0p5pct, FALSE),
  fix0p5pct1_withDoublets = estimate_contamination(fix0p5pct1, TRUE),
  fix0p5pct_withDoublets = estimate_contamination(fix0p5pct, TRUE)
)

# Print the named vector
print(contam_rates)
write.csv(contam_rates, "contam_rates_SVZ_unfiltered.csv")

# Create the data frame
df <- data.frame(
  Condition = rep(c("0,5%", "0,5% + mitoBender"), 2),
  Contamination = contam_rates,
  Doublets = rep(c("all cells", "doublets removed"), each = 2)
)
df
color.palette_fixation <- c("0,5%"= 'grey90', 
                            "0,5% + mitoBender"= 'lightblue')

# Plot the bar diagram
p1 <- ggplot(df, aes(x = Doublets, y = Contamination, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge2(padding = 0.2), width = 0.9,color = "grey50") +
  scale_fill_manual(values = color.palette_fixation) +
  labs(
    y = "Contamination (%)",
    x = "Condition"
  ) +
  theme_minimal() +
  ylim(0,12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  )
p1
cowplot::ggsave2(p1, width = 2, height = 5, dpi=300,  file = "barplot_mixing_mixing.pdf")
