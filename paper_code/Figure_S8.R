#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.04
# Supplementary Figure 8

# (a) Outlier shift values for haplotype-resolved splice junctions in genes with heterozygous or hemizygous pathogenic variants 
# in previously diagnosed patients. For a given splice junction, an outlier shift value is computed as the absolute difference 
# between its haplotype-resolved usage frequency in a sample and its average usage frequency across tissue-matched GTEx controls. 
# Splice junctions supported by at least 20 haplotype-resolved reads in a given sample are displayed and organized based on whether 
# they are proximal (i.e. within 100 bp) or distal to known disease-causing variants in cis. (b) Same as in (a) but 
# organized based on whether a variant is predicted to impact splicing (i.e. SpliceAI score > 0.1) rather than variant class. 

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
inDF <- bind_rows(
    inDF %>% filter(variant_class == "Structural deletion") %>% separate(variant_grch38, c(NA, "Variant_Start", "Variant_End")) %>%
        separate(hap_junction_coord, c(NA, "Junction_Start", "Junction_End")) %>% select(hap_junction_usage_shift, variant_class,
        spliceai, Variant_Start, Variant_End, Junction_Start, Junction_End),
    inDF %>% filter(variant_class != "Structural deletion") %>% separate(variant_grch38, c(NA, "Variant_Start", NA, NA)) %>% 
        mutate(Variant_End = Variant_Start) %>% separate(hap_junction_coord, c(NA, "Junction_Start", "Junction_End")) %>% 
        select(hap_junction_usage_shift, variant_class, spliceai, Variant_Start, Variant_End, Junction_Start, Junction_End)
) %>% mutate(Distance_Group = ifelse(pmax(as.integer(Variant_Start), as.integer(Junction_Start)) - pmin(as.integer(Variant_End),
    as.integer(Junction_End)) <= 100, "Proximal (i.e. within 100 bp)", "Distal"))

# =====================================================================================================================
#                                                      PANEL A
# =====================================================================================================================

classLevels <- inDF %>% filter(Distance_Group == "Proximal (i.e. within 100 bp)") %>% mutate(hap_junction_usage_shift = abs(hap_junction_usage_shift)) %>%
    group_by(variant_class) %>% summarise(Median_Shift = mean(hap_junction_usage_shift)) %>% arrange(Median_Shift) %>% pull(variant_class) %>% as.character
inDF$variant_class <- factor(as.character(inDF$variant_class), levels = classLevels)
inDF$Distance_Group <- factor(as.character(inDF$Distance_Group), levels = c("Distal", "Proximal (i.e. within 100 bp)"))

p1 <- ggplot() + geom_boxplot(data = inDF, aes(x = abs(hap_junction_usage_shift) * 100, y = variant_class, fill = Distance_Group), 
    linewidth = 0.5, color = "black", width = 0.5, outlier.alpha = 0.8, outlier.stroke = NA) + theme_bw() +
    xlab("Outlier shift in splice junction usage frequency (%)") + theme(panel.background = element_blank(), 
    axis.text = element_text(color = "black", size = 6), axis.title.x = element_text(color = "black", size = 7), axis.title.y = element_blank(), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), legend.position = "none") + scale_fill_manual(values = c("#7E5D81", "#BDCEDB")) +
    coord_cartesian(xlim = c(0, 100))

# =====================================================================================================================
#                                                      PANEL B
# =====================================================================================================================

inDF$Category <- factor(ifelse(inDF$spliceai > 0.1, "SpliceAI > 0.1", "SpliceAI <= 0.1"), levels = c("SpliceAI <= 0.1", "SpliceAI > 0.1"))

# Wilcoxon rank-sum test to assess whether outlier shift values for junctions from haplotypes with nearby splice-disrupting variants
# are similar with those from haplotypes with nearby non-splice-disrupting variants
#
#   wilcox.test(filter(inDF, Category == "SpliceAI > 0.1" & Distance_Group == "Proximal (i.e. within 100 bp)") %>% pull(hap_junction_usage_shift) %>% 
#       abs, filter(inDF, Category == "SpliceAI <= 0.1" & Distance_Group == "Proximal (i.e. within 100 bp)") %>% pull(hap_junction_usage_shift) %>% abs)
#       
#
#   p-value = 4.093e-09

p2 <- ggplot() + geom_boxplot(data = inDF %>% drop_na, aes(x = abs(hap_junction_usage_shift) * 100, y = Category, fill = Distance_Group), 
    linewidth = 0.5, color = "black", width = 0.5, outlier.alpha = 0.8, outlier.stroke = NA) + theme_bw() +
    xlab("Outlier shift in splice junction usage frequency (%)") + theme(panel.background = element_blank(), 
    axis.text = element_text(color = "black", size = 6), axis.title.x = element_text(color = "black", size = 7), axis.title.y = element_blank(), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), legend.position = "bottom", legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 6)) + scale_fill_manual(values = c("#7E5D81", "#BDCEDB")) + coord_cartesian(xlim = c(0, 100))

# Assemble p1 and p2 onto the same plotting grid
p <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1.5, 1), labels = c("a", "b"), label_size = 8, align = "v")
ggsave(outfile, plot = p, width = 3.5, height = 4.5)
