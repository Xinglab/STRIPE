#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.01
# Supplementary Figure 3

# (a, b) Comparison of (a) target gene abundances and (b) usage frequencies for splice junctions in target genes between 
# untargeted long-read RNA-seq and TEQUILA-seq (CDG-466 gene panel) data for fibroblast cell line CDG-152-1. Genes with 
# TPM ≥ 1 in untargeted data are displayed in (a), and splice junctions that have a usage frequency ≥ 5% and are supported 
# by ≥ 5 reads in untargeted data are displayed in (b). (c, d) Same as in (a, b) but between untargeted long-read RNA-seq 
# and TEQUILA-seq (PMD-359 gene panel) data for fibroblast cell line E1877. 

# =====================================================================================================================
#                                                      LIBRARIES 
# =====================================================================================================================

# Load required libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(cowplot))

# =====================================================================================================================
#                                                   HELPER FUNCTIONS
# =====================================================================================================================

PullFeature <- function(infoString, featureName) {
    # Function designed to pull out the value for featureName in infoString
    present <- unlist(lapply(strsplit(infoString, "; "), function(x) lapply(strsplit(x, " "), "[[", 1) == featureName))
    return(ifelse(sum(present) == 1, unlist(lapply(strsplit(unlist(strsplit(infoString, "; "))[present], " "), "[[", 2)), NA))
}

ComputePSI <- function(countDF) {
    # Function designed to compute usage frequencies for splice junctions in countDF
    startCoverage <- select(countDF, !(End)) %>% group_by(Chrom, Start) %>% mutate_if(is.numeric, sum) %>% ungroup %>% select(-c(Chrom, Start))
    endCoverage <- select(countDF, !(Start)) %>% group_by(Chrom, End) %>% mutate_if(is.numeric, sum) %>% ungroup %>% select(-c(Chrom, End))
    totalCoverage <- startCoverage + endCoverage - select(countDF, -c(Chrom, Start, End))
    return(bind_cols(countDF %>% select(c(Chrom, Start, End)), select(countDF, -c(Chrom, Start, End)) / mutate(totalCoverage, across(everything(), 
        ~ ifelse(. < 20, NA, .)))))
}

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir <- "/mnt/isilon/lin_lab_share/STRIPE"
outfile <- file.path(workdir, "manuscript/Supplementary_Figures/Figure_S3/Figure_S3.pdf")

# =====================================================================================================================
#                                                       PANEL A
# =====================================================================================================================

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
p1 <- ggplot(outDF, aes(x = log2(Untargeted + 1), y = log2(TEQUILA + 1))) + theme_classic() + geom_segment(aes(x = 0, xend = 15, y = 0, yend = 15), 
    linetype = "dashed", linewidth = 0.25) + geom_point(stroke = NA, size = 1, alpha = 0.5, color = "#0571B0") + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    plot.title = element_text(color = "black", size = 7, hjust = 0.5)) + xlab("Untargeted, log2(TPM+1)") + ylab("TEQUILA-seq, log2(TPM+1)") + 
    ggtitle("Gene abundance (CDG-466)") + geom_text(data = corrDF, mapping = aes(x = 5.25, y = 3.5, label = Label), size = 5/14*6, hjust = 0)

# =====================================================================================================================
#                                                       PANEL B
# =====================================================================================================================

# Calculate usage frequencies for splice junctions mapping to target genes in untargeted long-read RNA-seq and TEQUILA-seq data on fibroblast cell line CDG-152-1
untargeted.junctions <- read.table(file.path(workdir, "CDG/CDG-152-1/RNA/stripe_qc/quality_control/CDG-152-1_ONT.junctions.tsv"),
    sep = "\t", header = TRUE) %>% rename(Untargeted = Read_Counts)
targeted.junctions <- read.table(file.path(workdir, "CDG/CDG-152-1/RNA/stripe/quality_control/CDG-152-1_TEQUILA.junctions.tsv"),
    sep = "\t", header = TRUE) %>% rename(TEQUILA = Read_Counts)
countDF <- full_join(untargeted.junctions, targeted.junctions, by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)
psiDF <- ComputePSI(countDF)

# Only keep junctions supported by at least 5 reads and have a usage frequency >= 5% in untargeted long-read RNA-seq
outDF <- psiDF[countDF$Untargeted >= 5, ] %>% filter(Untargeted >= 0.05) %>% drop_na %>% tibble

# Compute Spearman's correlation in junction usage frequencies between untargeted long-read RNA-seq and TEQUILA-seq
corrDF <- outDF %>% summarise(Corr = cor(Untargeted, TEQUILA, method = "spearman")) %>%
    mutate(Label = paste("Spearman's r = ", format(round(Corr, digits = 2), nsmall = 2), sep = ""))

# Plot usage frequencies for splice junctions mapping to target genes between untargeted long-read RNA-seq and TEQUILA-seq
p2 <- ggplot(outDF, aes(x = Untargeted * 100, y = TEQUILA * 100)) + theme_classic() + geom_point(stroke = NA, size = 1, alpha = 0.5, color = "#0571B0") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.text = element_text(color = "black", size = 6), axis.title = element_text(color = "black", size = 7), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), plot.title = element_text(color = "black", size = 7, hjust = 0.5)) + 
    xlab("Untargeted (%)") + ylab("TEQUILA-seq (%)") + ggtitle("Splice junction usage (CDG-466)") + coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    geom_text(data = corrDF, mapping = aes(x = 35, y = 10, label = Label), size = 5/14*6, hjust = 0)

# =====================================================================================================================
#                                                       PANEL C
# =====================================================================================================================

# Read in target gene abundances for untargeted long-read RNA-seq data on fibroblast cell line E1877
pmd.genes <- read.table(file.path(workdir, "PMD/references/target_genes.txt"))$V1
untargeted.map <- read.table(file.path(workdir, "PMD/E1877/RNA/stripe_qc/quality_control/stringtie/E1877_ONT.gtf"),
    sep = "\t", header = FALSE) %>% filter(V3 == "transcript") %>% mutate(Gene_ID = unlist(lapply(V9, function(x) 
    PullFeature(x, "gene_id"))), Reference_Gene_ID = unlist(lapply(V9, function(x) PullFeature(x, "ref_gene_id")))) %>%
    select(Gene_ID, Reference_Gene_ID) %>% drop_na()
untargeted.data <- read.table(file.path(workdir, "PMD/E1877/RNA/stripe_qc/quality_control/stringtie/E1877_ONT.gene_abundance.tsv"),
    sep = "\t", header = TRUE) %>% filter(!(Gene.ID %in% untargeted.map$Gene_ID)) %>% separate(Gene.ID, c("Gene_ID", NA), sep = "\\.") %>% 
    filter(Gene_ID %in% pmd.genes) %>% select(Gene_ID, TPM) %>% mutate(Library = "Untargeted")

# Read in gene abundances for TEQUILA-seq data on fibroblast cell line E1877
targeted.map <- read.table(file.path(workdir, "PMD/E1877/RNA/stripe/quality_control/stringtie/E1877_TEQUILA.gtf"),
    sep = "\t", header = FALSE) %>% filter(V3 == "transcript") %>% mutate(Gene_ID = unlist(lapply(V9, function(x) 
    PullFeature(x, "gene_id"))), Reference_Gene_ID = unlist(lapply(V9, function(x) PullFeature(x, "ref_gene_id")))) %>%
    select(Gene_ID, Reference_Gene_ID) %>% drop_na()
targeted.data <- read.table(file.path(workdir, "PMD/E1877/RNA/stripe/quality_control/stringtie/E1877_TEQUILA.gene_abundance.tsv"),
    sep = "\t", header = TRUE) %>% filter(!(Gene.ID %in% untargeted.map$Gene_ID)) %>% separate(Gene.ID, c("Gene_ID", NA), sep = "\\.") %>% 
    filter(Gene_ID %in% pmd.genes) %>% select(Gene_ID, TPM) %>% mutate(Library = "TEQUILA")

# Filter for genes with TPM >= 1 in untargeted long-read RNA-seq
outDF <- bind_rows(untargeted.data, targeted.data) %>% mutate(Library = factor(Library, levels = c("Untargeted", "TEQUILA"))) %>%
    spread(Library, TPM) %>% filter(Untargeted >= 1)
    
# Compute Spearman's correlation in target gene abundances between untargeted long-read RNA-seq and TEQUILA-seq data
corrDF <- outDF %>% summarise(Corr = cor(Untargeted, TEQUILA, method = "spearman")) %>% mutate(Label = paste("Spearman's r = ", 
    format(round(Corr, digits = 2), nsmall = 2), sep = ""))

# Plot target gene abundances between untargeted long-read RNA-seq and TEQUILA-seq data
p3 <- ggplot(outDF, aes(x = log2(Untargeted + 1), y = log2(TEQUILA + 1))) + theme_classic() + geom_segment(aes(x = 0, xend = 15, y = 0, yend = 15), 
    linetype = "dashed", linewidth = 0.25) + geom_point(stroke = NA, size = 1, alpha = 0.5, color = "#CA0020") + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    plot.title = element_text(color = "black", size = 7, hjust = 0.5)) + xlab("Untargeted, log2(TPM+1)") + ylab("TEQUILA-seq, log2(TPM+1)") + 
    ggtitle("Gene abundance (PMD-359)") + geom_text(data = corrDF, mapping = aes(x = 5.25, y = 3.5, label = Label), size = 5/14*6, hjust = 0)

# =====================================================================================================================
#                                                       PANEL D
# =====================================================================================================================

# Calculate usage frequencies for splice junctions mapping to target genes in untargeted long-read RNA-seq and TEQUILA-seq data on fibroblast cell line E1877
untargeted.junctions <- read.table(file.path(workdir, "PMD/E1877/RNA/stripe_qc/quality_control/E1877_ONT.junctions.tsv"),
    sep = "\t", header = TRUE) %>% rename(Untargeted = Read_Counts)
targeted.junctions <- read.table(file.path(workdir, "PMD/E1877/RNA/stripe/quality_control/E1877_TEQUILA.junctions.tsv"),
    sep = "\t", header = TRUE) %>% rename(TEQUILA = Read_Counts)
countDF <- full_join(untargeted.junctions, targeted.junctions, by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)
psiDF <- ComputePSI(countDF)

# Only keep junctions supported by at least 5 reads and have a usage frequency >= 5% in untargeted long-read RNA-seq
outDF <- psiDF[countDF$Untargeted >= 5, ] %>% filter(Untargeted >= 0.05) %>% drop_na %>% tibble

# Compute Spearman's correlation in junction usage frequencies between untargeted long-read RNA-seq and TEQUILA-seq
corrDF <- outDF %>% summarise(Corr = cor(Untargeted, TEQUILA, method = "spearman")) %>%
    mutate(Label = paste("Spearman's r = ", format(round(Corr, digits = 2), nsmall = 2), sep = ""))

# Plot usage frequencies for splice junctions mapping to target genes between untargeted long-read RNA-seq and TEQUILA-seq
p4 <- ggplot(outDF, aes(x = Untargeted * 100, y = TEQUILA * 100)) + theme_classic() + geom_point(stroke = NA, size = 1, alpha = 0.5, color = "#CA0020") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.text = element_text(color = "black", size = 6), axis.title = element_text(color = "black", size = 7), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), plot.title = element_text(color = "black", size = 7, hjust = 0.5)) + 
    xlab("Untargeted (%)") + ylab("TEQUILA-seq (%)") + ggtitle("Splice junction usage (PMD-359)") + coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    geom_text(data = corrDF, mapping = aes(x = 35, y = 10, label = Label), size = 5/14*6, hjust = 0)

# Assemble plots into the same grid and save to outfile
p <- plot_grid(p1, p2, p3, p4, nrow = 2, labels = c("a", "b", "c", "d"), label_size = 8)
ggsave(outfile, plot = p, width = 4, height = 4)
