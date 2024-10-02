#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.01
# Supplementary Figure 2

# Gene abundances based on untargeted long-read RNA-seq and TEQUILA-seq (PMD-359 gene panel) data for fibroblast cell line E1877. 
# Each bar represents one gene, and only the 2000 most abundant genes are shown. Bars are color coded as genes that were targeted (red)
# or not targeted (gray). On-target rates represent the fraction of transcriptional abundance corresponding to target genes. Fold enrichment 
# is calculated by dividing the on-target rate in the capturant by the on-target rate in the unenriched input. 

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
outfile <- file.path(workdir, "manuscript/Supplementary_Figures/Figure_S2/Figure_S2.pdf")

# Read in gene abundances for untargeted long-read RNA-seq data on fibroblast cell line E1877
cdg.genes <- read.table(file.path(workdir, "PMD/references/target_genes.txt"))$V1
untargeted.map <- read.table(file.path(workdir, "PMD/E1877/RNA/stripe_qc/quality_control/stringtie/E1877_ONT.gtf"),
    sep = "\t", header = FALSE) %>% filter(V3 == "transcript") %>% mutate(Gene_ID = unlist(lapply(V9, function(x) 
    PullFeature(x, "gene_id"))), Reference_Gene_ID = unlist(lapply(V9, function(x) PullFeature(x, "ref_gene_id")))) %>%
    select(Gene_ID, Reference_Gene_ID) %>% drop_na()
untargeted.data <- read.table(file.path(workdir, "PMD/E1877/RNA/stripe_qc/quality_control/stringtie/E1877_ONT.gene_abundance.tsv"),
    sep = "\t", header = TRUE) %>% filter(!(Gene.ID %in% untargeted.map$Gene_ID)) %>% separate(Gene.ID, c("Gene_ID", NA), sep = "\\.") %>% 
    mutate(Group = case_when(Gene_ID %in% cdg.genes ~ "PMD-359", TRUE ~ "Other")) %>% select(Gene_ID, Group, TPM) %>% 
    mutate(Rank = rank(-TPM, ties.method = "first"), Library = "Untargeted")

# Read in gene abundances for TEQUILA-seq data on fibroblast cell line E1877
targeted.map <- read.table(file.path(workdir, "PMD/E1877/RNA/stripe/quality_control/stringtie/E1877_TEQUILA.gtf"),
    sep = "\t", header = FALSE) %>% filter(V3 == "transcript") %>% mutate(Gene_ID = unlist(lapply(V9, function(x) 
    PullFeature(x, "gene_id"))), Reference_Gene_ID = unlist(lapply(V9, function(x) PullFeature(x, "ref_gene_id")))) %>%
    select(Gene_ID, Reference_Gene_ID) %>% drop_na()
targeted.data <- read.table(file.path(workdir, "PMD/E1877/RNA/stripe/quality_control/stringtie/E1877_TEQUILA.gene_abundance.tsv"),
    sep = "\t", header = TRUE) %>% filter(!(Gene.ID %in% untargeted.map$Gene_ID)) %>% separate(Gene.ID, c("Gene_ID", NA), sep = "\\.") %>% 
    mutate(Group = case_when(Gene_ID %in% cdg.genes ~ "PMD-359", TRUE ~ "Other")) %>% select(Gene_ID, Group, TPM) %>% 
    mutate(Rank = rank(-TPM, ties.method = "first"), Library = "TEQUILA")

# Compute on-target rates for untargeted long-read RNA-seq and TEQUILA-seq
outDF <- bind_rows(untargeted.data, targeted.data) %>% mutate(Library = factor(Library, levels = c("Untargeted", "TEQUILA")), 
    Group = factor(Group, levels = c("PMD-359", "Other")))
on.target.rates <- outDF %>% filter(Group == "PMD-359") %>% group_by(Library) %>% reframe(Rate = sum(TPM)/1e4) %>% 
    mutate(Label = paste("On-target: ", formatC(round(Rate, digits = 1), format = "f", flag = "0", digits = 1), "%", sep = ""))

# Compute the enrichment fold attained by TEQUILA-seq
enrichment.fold <- on.target.rates %>% select(!(Label)) %>% spread(Library, Rate) %>% mutate(Enrichment = TEQUILA/Untargeted) %>% 
    mutate(Label = paste("Enrichment: ", formatC(round(Enrichment, digits = 1), format = "f", flag = "0", digits = 1), "x", sep = "")) %>% 
    mutate(Library = factor("TEQUILA", levels = c("Untargeted", "TEQUILA")))

# Plot gene rank versus gene abundance for the 2000 most abundant genes in untargeted long-read RNA-seq and TEQUILA-seq
palette <- setNames(c("#CA0020", "gray80"), c("PMD-359", "Other"))
p <- ggplot(outDF %>% filter(Rank <= 2000), aes(x = Rank, y = log2(TPM + 1))) + facet_wrap(vars(Library), nrow = 1) + theme_classic() +
    geom_bar(aes(fill = Group), linewidth = 0, stat = "identity") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), axis.title = element_text(color = "black", size = 7), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), strip.text = element_text(color = "black", size = 6), legend.title = element_blank(), 
    legend.text = element_text(color = "black", size = 6), legend.key.size = unit(0.15, 'cm'), legend.position = "bottom", 
    legend.background = element_rect(color = "black", linewidth = 0.25, fill = NA)) + xlab("Genes ranked by abundance") + ylab("Gene abundance\nlog2(TPM+1)") +
    geom_text(data = on.target.rates, aes(x = 600, y = 12, label = Label), size = 5/14*6, hjust = 0) + scale_fill_manual(values = palette) +
    geom_text(data = enrichment.fold, aes(x = 600, y = 10, label = Label), size = 5/14*6, hjust = 0) 

# Save p to outfile
ggsave(outfile, plot = p, width = 3, height = 2)
