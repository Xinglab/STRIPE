#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.02
# Supplementary Figure 6

# (a) Distances of splice junctions to known pathogenic variants in cis and their haplotype-resolved usage frequency shifts
# relative to tissue-matched GTEx controls based on TEQUILA-seq data on previously diagnosed individuals. Variants located in 
# the -3 to +6 region of the donor splice site or the -20 to +3 region of the acceptor splice site were classified as "splice 
# site region" variants. For plotting purposes, we represented the complex NUBPL haplotype (NM_025152.3:c.166G>A, p.Gly56Arg + 
# c.815-27T>C) in patient PMD-P07 by the missense variant alone. (b) Same as in (a) but organized based on whether a variant 
# is predicted to impact splicing (i.e. SpliceAI > 0.1) rather than variant class.

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
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Supplementary_Figures/Figure_S6/Figure_S6.pdf")

# Read in Supplementary Table 7 as a dataframe
inDF <- read.table(file.path(workdir, "manuscript/Revisions/20250228/Supplementary_Tables/Table_S7/Table_S7.txt"), sep = "\t", header = TRUE) %>%
    mutate(Shift = abs(Usage - GTEx_Usage))
labelDF <- table(inDF$Variant_Class) %>% data.frame %>% mutate(Label = paste(Var1, " (n = ", as.character(Freq), ")", sep = ""),
    Var1 = factor(Var1, levels = c("Structural deletion", "Nonsense", "Splice site region", "Inframe deletion", "Missense"))) %>%
    arrange(Var1)
inDF <- inDF %>% mutate(Variant_Class = factor(Variant_Class, levels = levels(labelDF$Var1)), Junction_Distance = pmax(1, Junction_Distance))
levels(inDF$Variant_Class) <- labelDF$Label

p1 <- ggplot(inDF, aes(x = Junction_Distance, y = Shift * 100, color = Variant_Class)) + geom_point(stroke = NA, alpha = 0.75, size = 1) + 
    facet_wrap(. ~ Variant_Class, nrow = 2) + scale_x_continuous(trans = "log10", breaks = c(1, 10, 100, 1000, 10000, 1e5, 1e6)) + theme_classic() +
    xlab("Distance between splice junction and known pathogenic variant in cis (bp)") + ylab("Absolute shift in splice junction usage frequency (%)") + 
    theme(panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), axis.title = element_text(color = "black", size = 7), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), axis.line = element_line(color = "black", linewidth = 0.25), 
    strip.text = element_text(color = "black", size = 6), legend.position = "none") + scale_color_manual(values = c("#446C68", "#DC626E", "#388BBF", 
    "#FCAD65", "#66C2A5"))

# Wilcoxon rank-sum test to assess whether usage frequency shifts for junctions from haplotypes with splice-disrupting variants within 100 bp away
# are similar with those from haplotypes with nearby non-splice-disrupting variants
#
#   wilcox.test(inDF %>% drop_na %>% filter(Junction_Distance <= 100 & SpliceAI > 0.1) %>% pull(Shift),
#       inDF %>% drop_na %>% filter(Junction_Distance <= 100 & SpliceAI <= 0.1) %>% pull(Shift))
#       
#   p-value = 6.4e-05

inDF2 <- inDF %>% drop_na %>% mutate(Group = ifelse(SpliceAI > 0.1, "SpliceAI > 0.1", "SpliceAI <= 0.1"))
labelDF2 <- table(inDF2$Group) %>% data.frame %>% mutate(Label = paste(Var1, " (n = ", as.character(Freq), ")", sep = ""))
p2 <- ggplot(, aes(x = Junction_Distance, y = Shift * 100, 
    color = Group)) + geom_point(stroke = NA, alpha = 0.75, size = 1) + facet_wrap(. ~ Group, nrow = 2) + scale_x_continuous(trans = "log10", breaks = c(1, 10, 
    100, 1000, 10000, 1e5, 1e6)) + theme_classic() + geom_vline(xintercept = 100, linetype = "dashed", color = "black", linewidth = 0.25) + 
    xlab("Variant-to-junction distance (bp)") + ylab("Absolute shift in splice junction usage frequency (%)") + 
    theme(panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), axis.title = element_text(color = "black", size = 7), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), axis.line = element_line(color = "black", linewidth = 0.25), strip.text = element_text(color = "black", size = 6), 
    legend.position = "none") + scale_color_manual(values = c("#225FA8", "#F868A2"))

# Assemble p1 and p2 onto the same plotting grid
p <- plot_grid(p1, p2, nrow = 1, rel_widths = c(2.5, 1), labels = c("a", "b"), label_size = 8, align = "h")
ggsave(outfile, plot = p, width = 5.5, height = 2.75)