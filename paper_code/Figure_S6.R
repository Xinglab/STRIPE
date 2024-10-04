#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.03
# Supplementary Figure 6

# Minor haplotype expression ratios for genes with heterozygous pathogenic variants in previously diagnosed patients.

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
medianDF <- aggregate(inDF$minor_ratio, list(inDF$minor_class), FUN = median) %>% setNames(c("minor_class", "median_ratio"))
classLevels <- arrange(medianDF, -median_ratio) %>% pull(minor_class)
inDF$minor_class <- factor(as.character(inDF$minor_class), levels = classLevels)

# Plot minor haplotype expression ratios for different classes of pathogenic variants
palette <- setNames(c("#F7A63D", "#63BA96", "#4C78B9", "#88CCEE", "#CE90BE", "#954492", "#E83578", "#223671"), 
    c("Missense", "Synonymous", "Splice donor", "Splice acceptor", "Inframe deletion", "Frameshift duplication", 
    "Nonsense", "Structural deletion"))

p <- ggplot() + geom_boxplot(data = inDF, aes(x = minor_ratio * 100, y = minor_class), linewidth = 0.5, width = 0.5, alpha = 0.5, color = "black") + 
    geom_jitter(data = inDF, aes(x = minor_ratio * 100, y = minor_class, color = minor_class), stroke = NA, alpha = 0.8, width = 0, height = 0.1, size = 2)  + 
    theme_bw() + xlab("Minor haplotype expression ratio (%)") + theme(panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title.x = element_text(color = "black", size = 7), axis.title.y = element_blank(), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    legend.position = "none") + coord_cartesian(xlim = c(0, 50)) + scale_color_manual(values = palette)

ggsave(outfile, plot = p, width = 3.5, height = 2)
