#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.04
# Supplementary Figure 8

# (a) Outlier shift values for haplotype-resolved splice junctions in genes with heterozygous or hemizygous pathogenic variants 
# in previously diagnosed patients. For a given splice junction, an outlier shift value is computed as the absolute difference 
# between its haplotype-resolved usage frequency in a sample and its average usage frequency across tissue-matched GTEx controls. 
# Splice junctions supported by at least 20 haplotype-resolved reads in a given sample are displayed, and splice junctions with
# outlier shift values less than 5% are greyed out. (b) Same as in (a) but organized based on whether a variant is predicted to 
# impact splicing (SpliceAI score > 0.1) rather than variant class. 

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
outfile <- file.path(workdir, "manuscript/Supplementary_Figures/Figure_S8/Figure_S8.pdf")

# Read in Supplementary Table 7 as a dataframe
inDF <- read.table(file.path(workdir, "manuscript/Supplementary_Tables/Table_S7/Table_S7.txt"), sep = "\t", header = TRUE)

# Flag junctions in inDF with absolute shift values < 5%
inDF$Minor <- abs(inDF$hap_junction_usage_shift) < 0.05

# =====================================================================================================================
#                                                      PANEL A
# =====================================================================================================================

medianDF <- aggregate(abs(inDF$hap_junction_usage_shift), list(inDF$variant_class), FUN = median) %>% setNames(c("variant_class", "median_shift"))
classLevels <- arrange(medianDF, median_shift) %>% pull(variant_class)
inDF$variant_class <- factor(as.character(inDF$variant_class), levels = classLevels)

palette <- setNames(c("#F7A63D", "#63BA96", "#4C78B9", "#88CCEE", "#CE90BE", "#954492", "#E83578", "#223671"), 
    c("Missense", "Synonymous", "Splice donor", "Splice acceptor", "Inframe deletion", "Frameshift duplication", 
    "Nonsense", "Structural deletion"))
p1 <- ggplot() + geom_boxplot(data = inDF, aes(x = abs(hap_junction_usage_shift) * 100, y = variant_class), linewidth = 0.5, width = 0.5, alpha = 0.5,
    color = "black", outlier.shape = NA) + geom_jitter(data = inDF %>% filter(Minor), aes(x = abs(hap_junction_usage_shift) * 100, y = variant_class), 
    color = "#BFBFBF", stroke = NA, alpha = 0.5, width = 0, height = 0.1) + geom_jitter(data = inDF %>% filter(!Minor), 
    aes(x = abs(hap_junction_usage_shift) * 100, y = variant_class, color = variant_class), stroke = NA, alpha = 0.8, width = 0, height = 0.1) +
    theme_bw() + xlab("Outlier shift in splice junction usage frequency (%)") + theme(panel.background = element_blank(), 
    axis.text = element_text(color = "black", size = 6), axis.title.x = element_text(color = "black", size = 7), axis.title.y = element_blank(), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), legend.position = "none") + scale_color_manual(values = palette) + 
    geom_vline(xintercept = 5, linewidth = 0.25, linetype = "dashed", color = "black") + xlim(0, 100)

# =====================================================================================================================
#                                                      PANEL B
# =====================================================================================================================

inDF$Category <- factor(ifelse(inDF$spliceai > 0.1, "Splice-disrupting", "Not splice-disrupting"),
    levels = c("Not splice-disrupting", "Splice-disrupting"))

# Wilcoxon rank-sum test to assess whether outlier shift values for junctions from haplotypes with splice-disrupting variants
# are similar with those from haplotypes with non-splice-disrupting variants
#
#   wilcox.test(filter(inDF, Category == "Splice-disrupting") %>% pull(hap_junction_usage_shift) %>% abs,
#       filter(inDF, Category == "Not splice-disrupting") %>% pull(hap_junction_usage_shift) %>% abs)
#
#   p-value = 1.190e-06

p2 <- ggplot() + geom_boxplot(data = inDF %>% drop_na, aes(x = abs(hap_junction_usage_shift) * 100, y = Category), linewidth = 0.5, width = 0.5, alpha = 0.5,
    color = "black", outlier.shape = NA) + geom_jitter(data = inDF %>% drop_na %>% filter(Minor), aes(x = abs(hap_junction_usage_shift) * 100, y = Category), 
    color = "#BFBFBF", stroke = NA, alpha = 0.5, width = 0, height = 0.1) + geom_jitter(data = inDF %>% drop_na %>% filter(!Minor), 
    aes(x = abs(hap_junction_usage_shift) * 100, y = Category), color = "#DF5C79", stroke = NA, alpha = 0.8, width = 0, height = 0.1) +
    theme_bw() + xlab("Outlier shift in splice junction usage frequency (%)") + theme(panel.background = element_blank(), 
    axis.text = element_text(color = "black", size = 6), axis.title.x = element_text(color = "black", size = 7), axis.title.y = element_blank(), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), legend.position = "none") + geom_vline(xintercept = 5, linewidth = 0.25, 
    linetype = "dashed", color = "black") + xlim(0, 100)

# Assemble p1 and p2 onto the same plotting grid
p <- plot_grid(p1, p2, ncol = 1, rel_heights = c(2.5, 1), labels = c("a", "b"), label_size = 8, align = "v")
ggsave(outfile, plot = p, width = 3.5, height = 3.5)
