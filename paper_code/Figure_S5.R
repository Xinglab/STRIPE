#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.10.08
# Supplementary Figure 5

# Number of pathogenic variants that reside on gene haplotypes with or without unusually lower expression relative to
# their opposite haplotypes, organized by variant class (Methods). Variants located in the -3 to +6 region of the donor 
# splice site or the -20 to +3 region of the acceptor splice site were classified as "splice site region" variants. For 
# plotting purposes, we represented the complex NUBPL haplotype (NM_025152.3:c.166G>A, p.Gly56Arg + c.815-27T>C) in 
# patient PMD-P07 by the missense variant alone.

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

# Construct a dataframe summarizing the number of variants that reside on gene haplotypes with or without unusually lower expression
# relative to their opposite haplotypes (Methods, based on Table S6):

df <- tibble(
    Class = c("Structural\ndeletion", "Nonsense", "Splice site\nregion", "Inframe\ndeletion", "Missense"),
    Count1 = c(1, 1, 3, 0, 1),
    Count2 = c(1, 0, 7, 3, 16)
) %>% mutate(Class = factor(Class, levels = Class)) %>% gather("Category", "Count", -Class) %>%
mutate(Category = ifelse(Category == "Count1", "Yes", "No")) %>% mutate(Category = factor(Category, levels = c("Yes", "No"))) %>%
mutate(PTV = Class %in% c("Structural\ndeletion", "Nonsense", "Splice site\nregion"))

p <- ggplot(df, aes(x = Class, y = Count, fill = Category)) + geom_bar(stat = "identity") + theme_classic() + ylab("Number of variants") +
    theme(panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), axis.title.y = element_text(color = "black", size = 7), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), axis.line = element_line(color = "black", linewidth = 0.25), 
    legend.position = "bottom", legend.text = element_text(color = "black", size = 6), legend.title = element_text(color = "black", size = 7),
    legend.key.size = unit(0.5, "cm"), axis.title.x = element_blank()) + guides(fill = guide_legend(title = "On gene haplotype with unusually lower\nexpression relative to opposite haplotype")) +
    ylim(0, 20) + scale_fill_manual(values = c("#cb2031", "#cbcbcb"))

ggsave(outfile, plot = p, width = 3.75, height = 2.5)

# Fisher exact test to assess whether PTVs are more like to be found on haplotypes showing significantly lower expression compared to the opposite haplotype
# compared to non-PTVs. Here, we define PTVs as variants that are nonsense, frameshift, splice site region, or structural deletion

fisher.test(df %>% group_by(PTV, Category) %>% summarize(Count = sum(Count)) %>% pull(Count) %>% matrix(., nrow = 2, ncol = 2, byrow = TRUE))$p.value # 0.025
fisher.test(df %>% group_by(PTV, Category) %>% summarize(Count = sum(Count)) %>% pull(Count) %>% matrix(., nrow = 2, ncol = 2, byrow = TRUE))$estimate # 10.9
