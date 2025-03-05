#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.05
# Figure 7e

# (e) Relative abundances of N-tetrasaccharide, a biomarker for ALG1-CDG, in serum samples from patient CDG-P22, 3 
# confirmed cases of ALG1-CDG, and 98 age-matched controls, as measured by MALDI-TOF mass spectrometry (at 1124 m/z).

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
infile <- file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_7/maldi-tof_ms.txt")
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_7/Figure_7e.pdf")

# Read in infile as a dataframe
inDF <- read.table(infile, sep = "\t", header = FALSE)
colnames(inDF) <- c("Value", "Group")
inDF$Group <- factor(inDF$Group, levels = c("ref", "ALG1", "case"))
levels(inDF$Group) <- c("Control", "ALG1-CDG", "CDG-P22")

# Plot relative abundances (1124 m/z) of N-tetrasaccharide in patient CDG-P22, confirmed cases of ALG1-CDG, and age-matched controls
p <- ggplot() + geom_jitter(data = inDF, aes(x = factor(1), y = Value, fill = Group, color = Group), shape = 21, width = 0.2, height = 0, size = 1) + theme_classic() +
    ylab("Relative abundance\nN-tetrasaccharide") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6), axis.title.y = element_text(color = "black", size = 7), axis.ticks.y = element_line(color = "black", linewidth = 0.25), 
    axis.line = element_line(color = "black", linewidth = 0.25), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
    legend.title = element_blank(), legend.text = element_text(color = "black", size = 6)) + scale_fill_manual(values = c("#BEBEBE", "#FDD976", "#F967A2")) +
    scale_color_manual(values = c(NA, NA, "#000000")) + ylim(0, 8)

ggsave(outfile, plot = p, width = 2, height = 1.5)
