#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.02.28
# Figure 1c

# Overview of individuals that underwent TEQUILA-seq. We sequenced skin fibroblast RNAs from 88 individuals, including 20 unaffected 
# healthy controls and 68 rare disease patients who had been clinically referred for evaluation of CDG (congenital disorders of 
# glycosylation) or PMD (primary mitochondrial diseases). Of the 68 rare disease patients, 22 had known genetic causes from prior 
# testing while the remaining 46 were genetically undiagnosed. 

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
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_1/Figure_1c.pdf")

# Construct dataframe summarizing number of individuals in cohort by disease group and genetic diagnosis category
df <- tibble(
    Category = rep(c("Control", "Known genetic cause", "Variant(s) of uncertain\nsignificance",
        "Single variant,\nrecessive condition", "No candidate genes"), each = 3),
    Group = rep(c("Control", "CDG", "PMD"), 5),
    Count = c(20, 0, 0, 0, 13, 9, 0, 1, 3, 0, 4, 0, 0, 11, 27)
)
df$Category <- factor(df$Category, levels = rev(c("Control", "Known genetic cause", "Variant(s) of uncertain\nsignificance",
    "Single variant,\nrecessive condition", "No candidate genes")))
df$Group <- factor(df$Group, levels = c("Control", "CDG", "PMD"))

# Plot number of individuals in cohort by disease group and genetic diagnosis category
p <- ggplot(df, aes(x = Count, y = Category, fill = Group)) + geom_bar(stat = "identity") + theme_classic() + xlim(0, 40) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.text = element_text(color = "black", size = 6), axis.title.x = element_text(color = "black", size = 7), 
    axis.title.y = element_blank(), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    axis.line = element_line(color = "black", linewidth = 0.25), legend.title = element_blank(), 
    legend.text = element_text(color = "black", size = 6), legend.position = c(0.8, 0.75), legend.key.size = unit(0.4, 'cm')) + 
    xlab("Number of individuals") + scale_fill_manual(values = c("#BEBEBE", "#3188BD", "#AADDA3"))

# Save p to outfile
ggsave(outfile, plot = p, width = 3, height = 1.75)
