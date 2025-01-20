#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.01.03
# Supplementary Figure 7

# (a) Distribution of usage frequency shifts (relative to tissue-matched GTEx controls) for splice junctions linked to haplotypes
# carrying pathogenic variants based on patient TEQUILA-seq reads. A splice junction is considered "proximal" if it is within 
# 200 bp of the pathogenic variant residing on the same haplotype. Variants located in the -3 to +6 region of the donor splice site 
# or the -20 to +3 region of the acceptor splice site were classified as "splice site region" variants. (b) Same as in (a) but organized 
# based on whether a variant is predicted to impact splicing (i.e. SpliceAI > 0.1) rather than variant class. P values were calculated 
# using a Wilcoxon rank-sum test.

# =====================================================================================================================
#                                                      LIBRARIES 
# =====================================================================================================================

# Load required libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(cowplot))

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir <- "/mnt/isilon/lin_lab_share/STRIPE"
outfile <- file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S7/Figure_S7.pdf")

# Read in Supplementary Table 7 as a dataframe
inDF <- read.table(file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Tables/Table_S7/Table_S7.txt"), sep = "\t", header = TRUE)
inDF <- inDF %>% mutate(SpliceAI = ifelse(Individual == "CDG-132-1" & Variant_Class == "Splice donor", 0.99, ifelse(Individual == "CDG-132-1" & 
    Variant_Class == "Inframe deletion", 0, SpliceAI))) %>% mutate(Variant_Class = ifelse(grepl("Splice", Variant_Class), "Splice site region",
    ifelse(grepl("Frameshift", Variant_Class), "Frameshift variant", Variant_Class)))
inDF <- inDF %>% mutate(Junction_Category = ifelse(Junction_Distance <= 200, "Proximal", "Distal"))

classLevels <- inDF %>% filter(Junction_Category == "Proximal") %>% mutate(Shift = abs(Usage - GTEx_Usage)) %>% group_by(Variant_Class) %>% 
    summarise(Mean_Shift = mean(Shift)) %>% arrange(Mean_Shift) %>% pull(Variant_Class) %>% as.character
inDF$Variant_Class <- factor(as.character(inDF$Variant_Class), levels = classLevels)
inDF$Junction_Category <- factor(as.character(inDF$Junction_Category), levels = c("Proximal", "Distal"))

p1 <- ggplot() + geom_boxplot(data = inDF, aes(x = abs(Usage - GTEx_Usage) * 100, y = Variant_Class, fill = Variant_Class), linewidth = 0.5, color = "black", width = 0.5, 
    outlier.alpha = 0.8, outlier.stroke = NA) + theme_classic() + xlab("Absolute shift in splice junction usage frequency (%)") + 
    theme(panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), axis.title.x = element_text(color = "black", size = 7), 
    axis.title.y = element_blank(), axis.ticks = element_line(color = "black", linewidth = 0.25), axis.line = element_line(color = "black", linewidth = 0.25),
    strip.text = element_text(color = "black", size = 6), legend.position = "none") + facet_wrap(. ~ Junction_Category, nrow = 1) +
    scale_fill_manual(values = c("#DC626E", "#FCAD65", "#66C2A5", "#CBB3D4", "#388BBF", "#446C68"))

inDF$Category <- factor(ifelse(inDF$SpliceAI > 0.1, "SpliceAI > 0.1", "SpliceAI <= 0.1"), levels = c("SpliceAI <= 0.1", "SpliceAI > 0.1"))

# Wilcoxon rank-sum test to assess whether usage frequency shifts for junctions from haplotypes with nearby splice-disrupting variants
# are similar with those from haplotypes with nearby non-splice-disrupting variants
#
#   wilcox.test(filter(inDF, Category == "SpliceAI > 0.1" & Junction_Category == "Proximal") %>% mutate(Shift = abs(Usage - GTEx_Usage)) %>% pull(Shift),
#       filter(inDF, Category == "SpliceAI <= 0.1" & Junction_Category == "Proximal") %>% mutate(Shift = abs(Usage - GTEx_Usage)) %>% pull(Shift))
#       
#
#   p-value = 2.4e-04

# Wilcoxon rank-sum test to assess whether usage frequency shifts for junctions from haplotypes with distal splice-disrupting variants
# are similar with those from haplotypes with distal non-splice-disrupting variants
#
#   wilcox.test(filter(inDF, Category == "SpliceAI > 0.1" & Junction_Category == "Distal") %>% mutate(Shift = abs(Usage - GTEx_Usage)) %>% pull(Shift),
#       filter(inDF, Category == "SpliceAI <= 0.1" & Junction_Category == "Distal") %>% mutate(Shift = abs(Usage - GTEx_Usage)) %>% pull(Shift))
#       
#
#   p-value = 0.14

p2 <- ggplot() + geom_boxplot(data = inDF %>% drop_na, aes(x = abs(Usage - GTEx_Usage) * 100, y = Category, fill = Category), linewidth = 0.5, color = "black", width = 0.5, 
    outlier.alpha = 0.8, outlier.stroke = NA) + theme_classic() + xlab("Absolute shift in splice junction usage frequency (%)") + 
    theme(panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), axis.title.x = element_text(color = "black", size = 7), 
    axis.title.y = element_blank(), axis.ticks = element_line(color = "black", linewidth = 0.25), axis.line = element_line(color = "black", linewidth = 0.25),
    strip.text = element_text(color = "black", size = 6), legend.position = "none") + facet_wrap(. ~ Junction_Category, nrow = 1) +
    scale_fill_manual(values = c("#D8DBE3", "#CF5444"))

# Assemble p1 and p2 onto the same plotting grid
p <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1.75, 1), labels = c("a", "b"), label_size = 8, align = "v")
ggsave(outfile, plot = p, width = 3.5, height = 2.75)