#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.01
# Figure 1d

# Comparison of target gene abundances (CDG-466 gene panel) between untargeted long-read RNA-seq and TEQUILA-seq data for fibroblast cell 
# line CDG-152-1. Only target genes with TPM ≥ 1 in untargeted data are displayed.

# =====================================================================================================================
#                                                      LIBRARIES 
# =====================================================================================================================

# Load required libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))

# =====================================================================================================================
#                                                   HELPER FUNCTIONS
# =====================================================================================================================

PullFeature <- function(infoString, featureName) {
    # Function designed to pull out the value for featureName in infoString
    present <- unlist(lapply(strsplit(infoString, "; "), function(x) lapply(strsplit(x, " "), "[[", 1) == featureName))
    return(ifelse(sum(present) == 1, unlist(lapply(strsplit(unlist(strsplit(infoString, "; "))[present], " "), "[[", 2)), NA))
}

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir <- "/mnt/isilon/lin_lab_share/STRIPE"
outfile <- file.path(workdir, "manuscript/Main_Figures/Figure_1/Figure_1d.pdf")

# Read in target gene abundances for untargeted long-read RNA-seq data on fibroblast cell line CDG-152-1
cdg.genes <- read.table(file.path(workdir, "CDG/references/target_genes.txt"))$V1
untargeted.map <- read.table(file.path(workdir, "CDG/CDG-152-1/RNA/stripe_qc/quality_control/stringtie/CDG-152-1_ONT.gtf"),
    sep = "\t", header = FALSE) %>% filter(V3 == "transcript") %>% mutate(Gene_ID = unlist(lapply(V9, function(x) 
    PullFeature(x, "gene_id"))), Reference_Gene_ID = unlist(lapply(V9, function(x) PullFeature(x, "ref_gene_id")))) %>%
    select(Gene_ID, Reference_Gene_ID) %>% drop_na()
untargeted.data <- read.table(file.path(workdir, "CDG/CDG-152-1/RNA/stripe_qc/quality_control/stringtie/CDG-152-1_ONT.gene_abundance.tsv"),
    sep = "\t", header = TRUE) %>% filter(!(Gene.ID %in% untargeted.map$Gene_ID)) %>% separate(Gene.ID, c("Gene_ID", NA), sep = "\\.") %>% 
    filter(Gene_ID %in% cdg.genes) %>% select(Gene_ID, TPM) %>% mutate(Library = "Untargeted")

# Read in gene abundances for TEQUILA-seq data on fibroblast cell line CDG-152-1
targeted.map <- read.table(file.path(workdir, "CDG/CDG-152-1/RNA/stripe/quality_control/stringtie/CDG-152-1_TEQUILA.gtf"),
    sep = "\t", header = FALSE) %>% filter(V3 == "transcript") %>% mutate(Gene_ID = unlist(lapply(V9, function(x) 
    PullFeature(x, "gene_id"))), Reference_Gene_ID = unlist(lapply(V9, function(x) PullFeature(x, "ref_gene_id")))) %>%
    select(Gene_ID, Reference_Gene_ID) %>% drop_na()
targeted.data <- read.table(file.path(workdir, "CDG/CDG-152-1/RNA/stripe/quality_control/stringtie/CDG-152-1_TEQUILA.gene_abundance.tsv"),
    sep = "\t", header = TRUE) %>% filter(!(Gene.ID %in% untargeted.map$Gene_ID)) %>% separate(Gene.ID, c("Gene_ID", NA), sep = "\\.") %>% 
    filter(Gene_ID %in% cdg.genes) %>% select(Gene_ID, TPM) %>% mutate(Library = "TEQUILA")

# Filter for genes with TPM >= 1 in untargeted long-read RNA-seq
outDF <- bind_rows(untargeted.data, targeted.data) %>% mutate(Library = factor(Library, levels = c("Untargeted", "TEQUILA"))) %>%
    spread(Library, TPM) %>% filter(Untargeted >= 1)
    
# Compute Spearman's correlation in target gene abundances between untargeted long-read RNA-seq and TEQUILA-seq data
corrDF <- outDF %>% summarise(Corr = cor(Untargeted, TEQUILA, method = "spearman")) %>% mutate(Label = paste("Spearman's r = ", 
    format(round(Corr, digits = 2), nsmall = 2), sep = ""))

# Plot target gene abundances between untargeted long-read RNA-seq and TEQUILA-seq data
p <- ggplot(outDF, aes(x = log2(Untargeted + 1), y = log2(TEQUILA + 1))) + theme_classic() + geom_segment(aes(x = 0, xend = 15, y = 0, yend = 15), 
    linetype = "dashed", linewidth = 0.25) + geom_point(stroke = NA, size = 1, alpha = 0.5, color = "#0571B0") + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    plot.title = element_text(color = "black", size = 7, hjust = 0.5)) + xlab("Untargeted, log2(TPM+1)") + ylab("TEQUILA, log2(TPM+1)") + 
    ggtitle("Gene abundance (CDG-466)") + geom_text(data = corrDF, mapping = aes(x = 5.25, y = 3.5, label = Label), size = 5/14*6, hjust = 0)

# Save p to outfile
ggsave(outfile, plot = p, width = 2, height = 2)
