#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.05
# Supplementary Figure 14

# Relative abundances of Neu5Ac1Hex1GlcNAc1-Asn, a biomarker for NGLY1 deficiency, in urine samples from patient CDG-P16, 
# 4 confirmed cases of NGLY1 deficiency, and 145 age-matched controls, as measured by MALDI-TOF mass spectrometry (at 990 m/z).

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
infile <- file.path(workdir, "manuscript/Revisions/20250228/Supplementary_Figures/Figure_S14/maldi-tof_ms.txt")
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Supplementary_Figures/Figure_S14/Figure_S14.pdf")

# Read in infile as a dataframe
inDF <- read.table(infile, sep = "\t", header = FALSE)
colnames(inDF) <- c("Value", "Group")
inDF$Group <- factor(inDF$Group, levels = c("ref", "NGLY1", "case"))
levels(inDF$Group) <- c("Control", "NGLY1 deficiency", "CDG-P16")

# Plot relative abundances (990 m/z) of Neu5Ac1Hex1GlcNAc1-Asn in patient CDG-P16, confirmed cases of NGLY1 deficiency, and age-matched controls
p <- ggplot() + geom_jitter(data = inDF, aes(x = factor(1), y = Value, fill = Group, color = Group), shape = 21, width = 0.2, height = 0, size = 1) + theme_classic() +
    ylab("Relative abundance\nNeu5Ac1Hex1GlcNAc1-Asn") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6), axis.title.y = element_text(color = "black", size = 7), axis.ticks.y = element_line(color = "black", linewidth = 0.25), 
    axis.line = element_line(color = "black", linewidth = 0.25), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
    legend.title = element_blank(), legend.text = element_text(color = "black", size = 6)) + scale_fill_manual(values = c("#BEBEBE", "#FDD976", "#F967A2")) +
    scale_color_manual(values = c(NA, NA, "#000000")) + ylim(0, 15)

ggsave(outfile, plot = p, width = 2.25, height = 1.5)
