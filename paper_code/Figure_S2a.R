#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.02.28
# Supplementary Figure 2a

# Genes represented in the CDG-466 (left) and PMD-359 (right) gene panels. Both panels include known disease genes as well 
# as genes that may be implicated in disease-associated pathways. 

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
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Supplementary_Figures/Figure_S2/Figure_S2a.pdf")

# Construct dataframe summarizing number of genes in the CDG-466 and PMD-359 gene panels
df <- tibble(
    Panel = rep(c("CDG-466", "PMD-359"), each = 6),
    Category = rep(c("Known CDG gene", "Glycosyltransferase/sulfotransferase", "Glycoside hydrolase",
        "Other glycosylation-related gene", "Known PMD gene", "Other nuclear-encoded mitochondrial gene"), 2),
    Count = c(183, 149, 80, 54, 0, 0, 0, 0, 0, 0, 292, 67)
)
df$Panel <- factor(df$Panel, levels = c("CDG-466", "PMD-359"))
levels(df$Panel) <- c("CDG-466 gene panel", "PMD-359 gene panel")
df$Category <- factor(df$Category, levels = c("Known CDG gene", "Glycosyltransferase/sulfotransferase", "Glycoside hydrolase",
        "Other glycosylation-related gene", "Known PMD gene", "Other nuclear-encoded mitochondrial gene"))
df <- df %>% group_by(Panel) %>% mutate(Fraction = Count/sum(Count), Max = cumsum(Fraction)) %>%
    mutate(Min = c(0, head(Max, -1)))

# Plot number of genes in the CDG-466 and PMD-359 gene panels
p <- ggplot(df, aes(ymax = Max, ymin = Min, xmax = 4, xmin = 3, fill = Category)) + geom_rect() + coord_polar(theta = "y") +
    xlim(c(2, 4)) + facet_wrap(. ~ Panel, nrow = 1) + theme_classic() + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_blank(),
    axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "none",
    strip.text = element_text(color = "black", size = 6)) + scale_fill_manual(values = c("#005579", "#8CB4CD", "#3E84AA",
    "#5FC7CF", "#8AC28E", "#CFE289"))

# Save p to outfile
ggsave(outfile, plot = p, width = 3, height = 1.75)