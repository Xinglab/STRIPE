#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.03
# Supplementary Figure 6

# Minor haplotype expression ratios for genes with heterozygous pathogenic variants in previously diagnosed patients. 
# The minor haplotype expression ratio is calculated as the haplotype with fewer reads divided by the total number of 
# reads across both haplotypes. Minor haplotypes are organized based on the class of variants that they carry. 

# =====================================================================================================================
#                                                      LIBRARIES 
# =====================================================================================================================

# Load required libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir <- "/mnt/isilon/lin_lab_share/STRIPE"
outfile <- file.path(workdir, "manuscript/Supplementary_Figures/Figure_S6/Figure_S6.pdf")

# Read in Table_S6.txt and compute haplotype expression ratios for the minor haplotype
inDF <- read.table(file.path(workdir, "manuscript/Supplementary_Tables/Table_S6/Table_S6.txt"), sep = "\t", 
    header = TRUE) %>% mutate(minor_class = ifelse(hap1_read_count > hap2_read_count, hap2_variant_class,
    hap1_variant_class), minor_ratio = ifelse(hap1_read_count > hap2_read_count, hap2_read_count/(hap1_read_count +
    hap2_read_count), hap1_read_count/(hap1_read_count + hap2_read_count)))

# Plot minor haplotype expression ratios for different classes of pathogenic variants
palette <- setNames(c("#F7A63D", "#63BA96", "#4C78B9", "#88CCEE", "#CE90BE", "#954492", "#E83578", "#223671"), 
    c("Missense", "Synonymous", "Splice donor", "Splice acceptor", "Inframe deletion", "Frameshift duplication", 
    "Nonsense", "Structural deletion"))

p <- ggplot(inDF %>% arrange(minor_ratio) %>% mutate(sample_rank = as.integer(rownames(.)))) + geom_segment(aes(x = sample_rank, xend = sample_rank, y = 0, 
    yend = minor_ratio * 100), color = "black", linewidth = 0.25) + geom_point(aes(x = sample_rank, y = minor_ratio * 100, color = minor_class), stroke = NA, size = 2) +
    theme_classic() + xlab("Known disease-causing alleles\non minor haplotypes") + ylab("Minor haplotype\nexpression ratio (%)") + theme(panel.background = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(color = "black", size = 6), axis.ticks.y = element_line(color = "black", 
    linewidth = 0.25), axis.title = element_text(color = "black", size = 7), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.text = element_text(color = "black",
    size = 6), legend.title = element_blank(), legend.key.size = unit(0.3, "cm")) + coord_cartesian(ylim = c(0, 50)) + scale_color_manual(values = palette)

ggsave(outfile, plot = p, width = 3.75, height = 1.5)
