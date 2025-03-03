#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.01
# Supplementary Figure 5a

# Expression ratios for haplotypes carrying pathogenic variants based on TEQUILA-seq data on previously diagnosed patients. 
# The expression ratio for a haplotype is calculated as the number of reads assigned to that haplotype divided by the total 
# number of reads assigned across both haplotypes. Variants located in the -3 to +6 region of the donor splice site or the 
# -20 to +3 region of the acceptor splice site were classified as "splice site region" variants. For plotting purposes, we
# represented the complex NUBPL haplotype (NM_025152.3:c.166G>A, p.Gly56Arg + c.815-27T>C) in patient PMD-P07 by the missense
# variant alone.

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
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Supplementary_Figures/Figure_S5/Figure_S5.pdf")

# Construct a dataframe with read counts for haplotypes carrying pathogenic variants in previously diagnosed patients.
df <- tibble(
    Sample_ID = c("CDG-P02", "CDG-P03", "CDG-P04", "CDG-P05", "CDG-P06", "CDG-P07", "CDG-P08", "CDG-P09", "CDG-P10",
        "CDG-P11", "CDG-P12", "CDG-P13", "PMD-P03", "PMD-P04", "PMD-P05", "PMD-P07", "PMD-P08", "PMD-P09"),
    Count1 = c(534, 407, 1788, 1012, 1489, 279, 401, 782, 3137, 3587, 8380, 2575, 622, 2228, 3329, 1017, 830, 1217),
    Count2 = c(928, 534, 2507, 303, 1519, 64, 483, 1592, 7431, 4936, 8309, 2862, 462, 1038, 4231, 1933, 859, 0),
    Class1 = c("Splice site region", "Missense", "Splice site region", "Missense", "Splice site region", "Missense",
        "Missense", "Splice site region", "Missense", "Missense", "Structural deletion", "Missense", 
        "Structural deletion", "Inframe deletion", "Missense", "Missense", "Splice site region", "Missense"),
    Class2 = c("Splice site region", "Missense", "Missense", "Splice site region", "Inframe deletion", "Nonsense",
        "Missense", "Inframe deletion", "Missense", "Missense", "Missense", "No pathogenic variant", "Missense",
        "Splice site region", "Splice site region", "Splice site region", "No pathogenic variant", "Structural deletion")
) %>% mutate(Ratio1 = Count1/(Count1+Count2), Ratio2 = Count2/(Count1+Count2), Index = pmin(Count1, Count2)/(Count1+Count2))
df$Sample_ID <- factor(df$Sample_ID, levels = df %>% arrange(Index) %>% pull(Sample_ID))
df <- bind_rows(df %>% select(Sample_ID, Class1, Ratio1) %>% setNames(c("Sample_ID", "Class", "Ratio")), 
    df %>% select(Sample_ID, Class2, Ratio2) %>% setNames(c("Sample_ID", "Class", "Ratio"))) %>% 
    mutate(Class = factor(Class, levels = c("Structural deletion", "Nonsense", "Splice site region",
    "Inframe deletion", "Missense", "No pathogenic variant")))

# Plot expression ratios for haplotypes carrying pathogenic variants based on TEQUILA-seq data on previously diagnosed patients
p <- ggplot(df, aes(x = Sample_ID, y = Ratio * 100, fill = Class)) + geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.25) +
    theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.text.x = element_text(color = "black", size = 6, angle = -30, hjust = 0), axis.title.y = element_text(color = "black", size = 7), 
    axis.text.y = element_text(color = "black", size = 6), axis.ticks = element_line(color = "black", linewidth = 0.25), axis.title.x = element_blank(), 
    axis.line = element_line(color = "black", linewidth = 0.25), legend.title = element_blank(), legend.text = element_text(color = "black", 
    size = 6), legend.key.size = unit(0.4, 'cm'), legend.position = "bottom") + ylab("Haplotype expression ratio (%)") + 
    scale_fill_manual(values = c("#446C68", "#DC626E", "#388BBF", "#FCAD65", "#66C2A5", "#FFFFFF"))

# Save p to outfile
ggsave(outfile, plot = p, width = 3.75, height = 2.5)