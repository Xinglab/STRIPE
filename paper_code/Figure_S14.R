#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.01.31
# Supplementary Figure 14

# Relative abundances of Neu5Ac1Hex1GlcNAc1-Asn, a biomarker for NGLY1 deficiency, in urine samples from patient FCDGC-02003, 
# confirmed cases of NGLY1 deficiency, and age-matched controls, as measured by MALDI-TOF mass spectrometry (at 990 m/z).

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
infile <- file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S14/maldi-tof_ms.txt")
outfile <- file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S14/Figure_S14.pdf")

# Read in infile as a dataframe
inDF <- read.table(infile, sep = "\t", header = FALSE)
colnames(inDF) <- c("Value", "Group")
inDF$Group <- factor(inDF$Group, levels = c("ref", "NGLY1", "case"))
levels(inDF$Group) <- c("Control", "NGLY1\ndeficiency", "Patient\nFCDGC-02003")
summaryDF <- inDF %>% group_by(Group) %>% summarise(Value = mean(Value, na.rm = T))

# Plot relative abundances (990 m/z) of Neu5Ac1Hex1GlcNAc1-Asn in patient FCDGC-02003, confirmed cases of NGLY1 deficiency, and age-matched controls
p <- ggplot() + geom_bar(data = summaryDF, aes(x = Group, y = Value, fill = Group), linewidth = 0.5, stat = "identity", width = 0.75) + 
    geom_point(data = inDF, aes(x = Group, y = Value, color = Group), stroke = NA, size = 0.75, alpha = 0.5, 
    position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0)) + theme_classic() + 
    ylab("Relative abundance\nNeu5Ac1Hex1GlcNAc1-Asn") + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title.y = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    axis.line = element_line(color = "black", linewidth = 0.25), axis.title.x = element_blank(), legend.position = "none") + ylim(0, 15) +
    scale_fill_manual(values = c("#B0D405", "#97DDE9", "#F8A500")) + scale_color_manual(values = c("black", "black", "black"))

ggsave(outfile, plot = p, width = 2.5, height = 2)
