#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.01
# Figure 1c

# Gene abundances based on untargeted long-read RNA-seq and TEQUILA-seq (CDG-466 gene panel) data for fibroblast cell line CDG-152-1. 
# Each bar represents one gene, and only the 2000 most abundant genes are shown. Bars are color coded as genes that were targeted (blue)
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
outfile <- file.path(workdir, "manuscript/Main_Figures/Figure_1/Figure_1c.pdf")

# Read in gene abundances for untargeted long-read RNA-seq data on fibroblast cell line CDG-152-1
cdg.genes <- read.table(file.path(workdir, "CDG/references/target_genes.txt"))$V1
untargeted.map <- read.table(file.path(workdir, "CDG/CDG-152-1/RNA/stripe_qc/quality_control/stringtie/CDG-152-1_ONT.gtf"),
    sep = "\t", header = FALSE) %>% filter(V3 == "transcript") %>% mutate(Gene_ID = unlist(lapply(V9, function(x) 
    PullFeature(x, "gene_id"))), Reference_Gene_ID = unlist(lapply(V9, function(x) PullFeature(x, "ref_gene_id")))) %>%
    select(Gene_ID, Reference_Gene_ID) %>% drop_na()
untargeted.data <- read.table(file.path(workdir, "CDG/CDG-152-1/RNA/stripe_qc/quality_control/stringtie/CDG-152-1_ONT.gene_abundance.tsv"),
    sep = "\t", header = TRUE) %>% filter(!(Gene.ID %in% untargeted.map$Gene_ID)) %>% separate(Gene.ID, c("Gene_ID", NA), sep = "\\.") %>% 
    mutate(Group = case_when(Gene_ID %in% cdg.genes ~ "CDG-466", TRUE ~ "Other")) %>% select(Gene_ID, Group, TPM) %>% 
    mutate(Rank = rank(-TPM, ties.method = "first"), Library = "Untargeted")

# Read in gene abundances for TEQUILA-seq data on fibroblast cell line CDG-152-1
targeted.map <- read.table(file.path(workdir, "CDG/CDG-152-1/RNA/stripe/quality_control/stringtie/CDG-152-1_TEQUILA.gtf"),
    sep = "\t", header = FALSE) %>% filter(V3 == "transcript") %>% mutate(Gene_ID = unlist(lapply(V9, function(x) 
    PullFeature(x, "gene_id"))), Reference_Gene_ID = unlist(lapply(V9, function(x) PullFeature(x, "ref_gene_id")))) %>%
    select(Gene_ID, Reference_Gene_ID) %>% drop_na()
targeted.data <- read.table(file.path(workdir, "CDG/CDG-152-1/RNA/stripe/quality_control/stringtie/CDG-152-1_TEQUILA.gene_abundance.tsv"),
    sep = "\t", header = TRUE) %>% filter(!(Gene.ID %in% untargeted.map$Gene_ID)) %>% separate(Gene.ID, c("Gene_ID", NA), sep = "\\.") %>% 
    mutate(Group = case_when(Gene_ID %in% cdg.genes ~ "CDG-466", TRUE ~ "Other")) %>% select(Gene_ID, Group, TPM) %>% 
    mutate(Rank = rank(-TPM, ties.method = "first"), Library = "TEQUILA")

# Compute on-target rates for untargeted long-read RNA-seq and TEQUILA-seq
outDF <- bind_rows(untargeted.data, targeted.data) %>% mutate(Library = factor(Library, levels = c("Untargeted", "TEQUILA")), 
    Group = factor(Group, levels = c("CDG-466", "Other")))
on.target.rates <- outDF %>% filter(Group == "CDG-466") %>% group_by(Library) %>% reframe(Rate = sum(TPM)/1e4) %>% 
    mutate(Label = paste("On-target: ", formatC(round(Rate, digits = 1), format = "f", flag = "0", digits = 1), "%", sep = ""))

# Compute the enrichment fold attained by TEQUILA-seq
enrichment.fold <- on.target.rates %>% select(!(Label)) %>% spread(Library, Rate) %>% mutate(Enrichment = TEQUILA/Untargeted) %>% 
    mutate(Label = paste("Enrichment: ", formatC(round(Enrichment, digits = 1), format = "f", flag = "0", digits = 1), "x", sep = "")) %>% 
    mutate(Library = factor("TEQUILA", levels = c("Untargeted", "TEQUILA")))

# Plot gene rank versus gene abundance for the 2000 most abundant genes in untargeted long-read RNA-seq and TEQUILA-seq
palette <- setNames(c("#0571B0", "gray80"), c("CDG-466", "Other"))
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
