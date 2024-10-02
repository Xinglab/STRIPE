#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.01
# Figure 1e

# Principal component analysis (PCA) based on abundance ranks of target genes (CDG-466 gene panel) in fibroblasts from 10 
# healthy controls, 29 CDG patients, and GTEx samples for tissues that can be collected from skin punch biopsies. 

# =====================================================================================================================
#                                                      LIBRARIES 
# =====================================================================================================================

# Load required libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(pcaMethods))

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
outfile <- file.path(workdir, "manuscript/Main_Figures/Figure_1/Figure_1e.pdf")

# Read in target gene abundances for TEQUILA-seq data (CDG-466 gene panel) on 39 individuals in our cohort
cdg.genes <- read.table(file.path(workdir, "CDG/references/target_genes.txt"))$V1
cohort.samples <- read.table(file.path(workdir, "CDG/samples.txt"), sep = "\t", header = TRUE) %>%
    filter(Provider != "Lan Lin") %>% pull(ID)
tissue.assign <- setNames(rep("Fibroblasts", length(cohort.samples)), cohort.samples)
outDF <- tibble(Gene_ID = cdg.genes)

for (sample.id in cohort.samples) {
    sample.map <- read.table(file.path(workdir, "CDG", sample.id, "RNA/stripe/quality_control/stringtie", paste(sample.id, "TEQUILA.gtf", sep = "_")),
        sep = "\t", header = FALSE) %>% filter(V3 == "transcript") %>% mutate(Gene_ID = unlist(lapply(V9, function(x) 
        PullFeature(x, "gene_id"))), Reference_Gene_ID = unlist(lapply(V9, function(x) PullFeature(x, "ref_gene_id")))) %>%
        select(Gene_ID, Reference_Gene_ID) %>% drop_na()
    sample.data <- read.table(file.path(workdir, "CDG", sample.id, "RNA/stripe/quality_control/stringtie", paste(sample.id, "TEQUILA.gene_abundance.tsv", sep = "_")),
        sep = "\t", header = TRUE) %>% filter(!(Gene.ID %in% sample.map$Gene_ID)) %>% separate(Gene.ID, c("Gene_ID", NA), sep = "\\.") %>% 
        filter(Gene_ID %in% cdg.genes) %>% select(Gene_ID, TPM) %>% setNames(c("Gene_ID", sample.id))
    outDF <- full_join(outDF, sample.data, by = join_by(Gene_ID)) %>% replace(is.na(.), 0)
}

# Read in target gene abundances for GTEx fibroblasts, skin, adipose, and whole blood
for (tissue in c("Fibroblast", "Skin", "Adipose", "Blood")) {
    tissueDF <- read.table(file.path(workdir, "CDG/references/GTEx_v8", tissue, "target_genes.gene_tpm.txt"),
        sep = "\t", header = TRUE, check.names = FALSE) %>% separate(Name, c("Gene_ID", NA), sep = "\\.")
    tissue.assign <- c(tissue.assign, setNames(rep(tissue, ncol(tissueDF) - 1), colnames(tissueDF)[-1]))
    outDF <- full_join(outDF, tissueDF, by = join_by(Gene_ID)) %>% replace(is.na(.), 0)
}

# Perform PCA clustering of samples based on target gene abundance ranks
pcaResults <- pca(outDF %>% mutate_if(is.numeric, function(x) rank(-x, ties.method = "max")) %>%
    tibble::column_to_rownames(var = "Gene_ID") %>% t, method = "nipals", nPcs = 2, scale = "none", center = TRUE)
scoreDF <- as.data.frame(scores(pcaResults)) %>% tibble::rownames_to_column(var = "Sample_ID") %>%
    mutate(Tissue = tissue.assign[Sample_ID]) %>% mutate(Tissue = factor(case_when(Tissue == "Blood" ~ "Whole blood (GTEx v8)",
    Tissue == "Adipose" ~ "Subcutaneous adipose (GTEx v8)", Tissue == "Skin" ~ "Skin - lower leg (GTEx v8)", 
    Tissue == "Fibroblast" ~ "Fibroblasts (GTEx v8)", TRUE ~ "Fibroblasts"), levels = c("Whole blood (GTEx v8)",
    "Subcutaneous adipose (GTEx v8)", "Skin - lower leg (GTEx v8)", "Fibroblasts (GTEx v8)", "Fibroblasts")))
palette <- setNames(c("#E43321", "#F9BD96", "#46B658", "#1A78AC", "#6FC7CE"), levels(scoreDF$Tissue))
p <- ggplot(scoreDF[order(scoreDF$Tissue), ], aes(x = PC1, y = PC2, color = Tissue)) + geom_point(size = 1, stroke = NA) + 
    theme_classic() + xlab(paste("PC1 (", format(round(pcaResults@R2[1]*100, digits = 1), nsmall = 1), "%)", sep = "")) +
    ylab(paste("PC2 (", format(round(pcaResults@R2[2]*100, digits = 1), nsmall = 1), "%)", sep = "")) + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), legend.position = "none", 
    plot.title = element_text(color = "black", size = 7, hjust = 0.5)) + scale_color_manual(values = palette) + coord_cartesian(xlim = c(-1500, 1500), 
    ylim = c(-1500, 1500)) + ggtitle("Gene abundance (CDG-466)")
    
# Save p to outfile
ggsave(outfile, plot = p, width = 2, height = 2)
