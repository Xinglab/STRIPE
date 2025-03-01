#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.01
# Supplementary Figure 4a

# Number of known pathogenic variants within exonic regions that were detected and phased by STRIPE using information
# from patient TEQUILA-seq reads alone.

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
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Supplementary_Figures/Figure_S4/Figure_S4a.pdf")

# Construct dataframe summarizing number of known pathogenic variants within exonic regions that were
# detected and phased by STRIPE
df <- tibble(
    Class = rep(c("Synonymous", "Nonsense", "Inframe deletion", "Missense"), each = 2),
    Category = rep(c("Detected, not phased", "Detected and phased"), 4),
    Count = c(0, 1, 0, 1, 2, 1, 3, 15)
)
df$Category <- factor(df$Category, levels = c("Detected, not phased", "Detected and phased"))
df$Class <- factor(df$Class, levels = c("Synonymous", "Nonsense", "Inframe deletion", "Missense"))

# Plot number of known pathogenic variants within exonic regions that were detected and phased by STRIPE
p <- ggplot(df, aes(x = Count, y = Class, fill = Category)) + geom_bar(stat = "identity") + theme_classic() + xlim(0, 20) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.text = element_text(color = "black", size = 6), axis.title.x = element_text(color = "black", size = 7), 
    axis.title.y = element_blank(), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    axis.line = element_line(color = "black", linewidth = 0.25), legend.title = element_blank(), 
    legend.text = element_text(color = "black", size = 6), legend.key.size = unit(0.4, 'cm')) + xlab("Number of variants") + 
    scale_fill_manual(values = c("#FFE699", "#00BA38"))

# Save p to outfile
ggsave(outfile, plot = p, width = 4, height = 1.5)
