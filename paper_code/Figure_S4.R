#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.02
# Supplementary Figure 4

# (a, b) Principal component analysis (PCA) based on (a) abundance ranks of target genes and (b) usage frequencies for splice 
# junctions in target genes (CDG-466 gene panel) for fibroblasts from 10 healthy controls, 29 CDG patients, and GTEx samples 
# for tissues that can be collected from skin punch biopsies. (c, d) Same as in (a, b) but for target genes (PMD-359 gene panel) 
# in fibroblasts from 1 healthy control, 12 PMD patients, and GTEx samples for tissues that can be collected from skin punch 
# biopsies. In (b, d), splice junctions meeting the following criteria were used for PCA: (i) no more than one-third of all 
# samples have missing usage frequencies, (ii) maximum usage frequency across all samples is at least 5%, and (iii) range of 
# usage frequencies across all samples is at least 5%. 

# =====================================================================================================================
#                                                      LIBRARIES 
# =====================================================================================================================

# Load required libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(cowplot))
suppressMessages(library(pcaMethods))

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
outfile <- file.path(workdir, "manuscript/Supplementary_Figures/Figure_S4/Figure_S4.pdf")

# =====================================================================================================================
#                                                       PANEL A
# =====================================================================================================================

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
p1 <- ggplot(scoreDF[order(scoreDF$Tissue), ], aes(x = PC1, y = PC2, color = Tissue)) + geom_point(size = 1, stroke = NA) + 
    theme_classic() + xlab(paste("PC1 (", format(round(pcaResults@R2[1]*100, digits = 1), nsmall = 1), "%)", sep = "")) +
    ylab(paste("PC2 (", format(round(pcaResults@R2[2]*100, digits = 1), nsmall = 1), "%)", sep = "")) + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), legend.position = "none", 
    plot.title = element_text(color = "black", size = 7, hjust = 0.5)) + scale_color_manual(values = palette) + coord_cartesian(xlim = c(-1500, 1500), 
    ylim = c(-1500, 1500)) + ggtitle("Gene abundance (CDG-466)")

# =====================================================================================================================
#                                                       PANEL B
# =====================================================================================================================

outDF <- tibble(Chrom = character(0), Start = integer(0), End = integer(0))

# Compute usage frequencies for splice junctions mapping to target genes for GTEx fibroblasts, skin, adipose, and whole blood
for (tissue in c("Fibroblast", "Skin", "Adipose", "Blood")) {
    tissueDF <- ComputePSI(read.table(file.path(workdir, "CDG/references/GTEx_v8", tissue, "target_genes.junction_counts.txt"), sep = "\t",
        header = TRUE, check.names = FALSE) %>% distinct() %>% separate(Name, c("Chrom", "Start", "End"), sep = "_") %>% mutate_at(c("Start", "End"), as.integer))
    outDF <- full_join(outDF, tissueDF, by = join_by(Chrom, Start, End))
}

# Compute usage frequencies for splice junctions mapping to target genes for TEQUILA-seq data (CDG-466 gene panel) on 39 individuals in our cohort
for (sample.id in cohort.samples) {
    sampleDF <- ComputePSI(read.table(file.path(workdir, "CDG", sample.id, "RNA/stripe/quality_control", paste(sample.id, "TEQUILA.junctions.tsv", 
        sep = "_")), sep = "\t", header = TRUE) %>% rename("{sample.id}" := Read_Counts))
    outDF <- left_join(outDF, sampleDF, by = join_by(Chrom, Start, End))
}

# Filter outDF for junctions in which (i) no more than one-third of the samples have missing usage frequencies, (ii) maximum
# usage frequency across samples is at least 5%, (iii) range of usage frequencies across samples is at least 5%. Then remove
# samples with missing usage frequencies for more than one-third of the remaining junctions
outDF <- outDF %>% unite("Name", Chrom:End) %>% tibble::column_to_rownames(var = "Name") %>% filter(rowSums(is.na(.))/ncol(.) <= 1/3) %>% 
    filter(apply(., 1, max, na.rm = TRUE) >= 0.05) %>% filter(apply(., 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) >= 0.05) %>%
    t %>% as.data.frame %>% filter(rowSums(is.na(.))/ncol(.) <= 1/3) %>% t %>% as.data.frame

# Perform PCA clustering of samples based on usage frequencies of splice junctions represented in outDF
pcaResults <- pca(t(outDF), method = "nipals", nPcs = 2, scale = "none", center = TRUE)
scoreDF <- as.data.frame(scores(pcaResults)) %>% tibble::rownames_to_column(var = "Sample_ID") %>%
    mutate(Tissue = tissue.assign[Sample_ID]) %>% mutate(Tissue = factor(case_when(Tissue == "Blood" ~ "Whole blood (GTEx v8)",
    Tissue == "Adipose" ~ "Subcutaneous adipose (GTEx v8)", Tissue == "Skin" ~ "Skin - lower leg (GTEx v8)", 
    Tissue == "Fibroblast" ~ "Fibroblasts (GTEx v8)", TRUE ~ "Fibroblasts"), levels = c("Whole blood (GTEx v8)",
    "Subcutaneous adipose (GTEx v8)", "Skin - lower leg (GTEx v8)", "Fibroblasts (GTEx v8)", "Fibroblasts")))
p2 <- ggplot(scoreDF[order(scoreDF$Tissue), ], aes(x = PC1, y = PC2, color = Tissue)) + geom_point(size = 1, stroke = NA) + 
    theme_classic() + xlab(paste("PC1 (", format(round(pcaResults@R2[1]*100, digits = 1), nsmall = 1), "%)", sep = "")) +
    ylab(paste("PC2 (", format(round(pcaResults@R2[2]*100, digits = 1), nsmall = 1), "%)", sep = "")) + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), legend.position = "none", 
    plot.title = element_text(color = "black", size = 7, hjust = 0.5)) + scale_color_manual(values = palette) + coord_cartesian(xlim = c(-3, 3), 
    ylim = c(-3, 3)) + ggtitle("Splice junction usage (CDG-466)")

# =====================================================================================================================
#                                                       PANEL C
# =====================================================================================================================

# Read in target gene abundances for TEQUILA-seq data (PMD-359 gene panel) on 13 individuals in our cohort
pmd.genes <- read.table(file.path(workdir, "PMD/references/target_genes.txt"))$V1
cohort.samples <- read.table(file.path(workdir, "PMD/samples.txt"), sep = "\t", header = TRUE) %>%
    filter(Provider == "Rebecca Ganetzky") %>% pull(ID)
tissue.assign <- setNames(rep("Fibroblasts", length(cohort.samples)), cohort.samples)
outDF <- tibble(Gene_ID = pmd.genes)

for (sample.id in cohort.samples) {
    sample.map <- read.table(file.path(workdir, "PMD", sample.id, "RNA/stripe/quality_control/stringtie", paste(sample.id, "TEQUILA.gtf", sep = "_")),
        sep = "\t", header = FALSE) %>% filter(V3 == "transcript") %>% mutate(Gene_ID = unlist(lapply(V9, function(x) 
        PullFeature(x, "gene_id"))), Reference_Gene_ID = unlist(lapply(V9, function(x) PullFeature(x, "ref_gene_id")))) %>%
        select(Gene_ID, Reference_Gene_ID) %>% drop_na()
    sample.data <- read.table(file.path(workdir, "PMD", sample.id, "RNA/stripe/quality_control/stringtie", paste(sample.id, "TEQUILA.gene_abundance.tsv", sep = "_")),
        sep = "\t", header = TRUE) %>% filter(!(Gene.ID %in% sample.map$Gene_ID)) %>% separate(Gene.ID, c("Gene_ID", NA), sep = "\\.") %>% 
        filter(Gene_ID %in% pmd.genes) %>% select(Gene_ID, TPM) %>% setNames(c("Gene_ID", sample.id))
    outDF <- full_join(outDF, sample.data, by = join_by(Gene_ID)) %>% replace(is.na(.), 0)
}

# Read in target gene abundances for GTEx fibroblasts, skin, adipose, and whole blood
for (tissue in c("Fibroblast", "Skin", "Adipose", "Blood")) {
    tissueDF <- read.table(file.path(workdir, "PMD/references/GTEx_v8", tissue, "target_genes.gene_tpm.txt"),
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
p3 <- ggplot(scoreDF[order(scoreDF$Tissue), ], aes(x = PC1, y = PC2, color = Tissue)) + geom_point(size = 1, stroke = NA) + 
    theme_classic() + xlab(paste("PC1 (", format(round(pcaResults@R2[1]*100, digits = 1), nsmall = 1), "%)", sep = "")) +
    ylab(paste("PC2 (", format(round(pcaResults@R2[2]*100, digits = 1), nsmall = 1), "%)", sep = "")) + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), legend.position = "none", 
    plot.title = element_text(color = "black", size = 7, hjust = 0.5)) + scale_color_manual(values = palette) + coord_cartesian(xlim = c(-1500, 1500), 
    ylim = c(-1500, 1500)) + ggtitle("Gene abundance (PMD-359)")

# =====================================================================================================================
#                                                       PANEL D
# =====================================================================================================================

outDF <- tibble(Chrom = character(0), Start = integer(0), End = integer(0))

# Compute usage frequencies for splice junctions mapping to target genes for GTEx fibroblasts, skin, adipose, and whole blood
for (tissue in c("Fibroblast", "Skin", "Adipose", "Blood")) {
    tissueDF <- ComputePSI(read.table(file.path(workdir, "PMD/references/GTEx_v8", tissue, "target_genes.junction_counts.txt"), sep = "\t",
        header = TRUE, check.names = FALSE) %>% distinct() %>% separate(Name, c("Chrom", "Start", "End"), sep = "_") %>% mutate_at(c("Start", "End"), as.integer))
    outDF <- full_join(outDF, tissueDF, by = join_by(Chrom, Start, End))
}

# Compute usage frequencies for splice junctions mapping to target genes for TEQUILA-seq data (PMD-359 gene panel) on 13 individuals in our cohort
for (sample.id in cohort.samples) {
    sampleDF <- ComputePSI(read.table(file.path(workdir, "PMD", sample.id, "RNA/stripe/quality_control", paste(sample.id, "TEQUILA.junctions.tsv", 
        sep = "_")), sep = "\t", header = TRUE) %>% rename("{sample.id}" := Read_Counts))
    outDF <- left_join(outDF, sampleDF, by = join_by(Chrom, Start, End))
}

# Filter outDF for junctions in which (i) no more than one-third of the samples have missing usage frequencies, (ii) maximum
# usage frequency across samples is at least 5%, (iii) range of usage frequencies across samples is at least 5%. Then remove
# samples with missing usage frequencies for more than one-third of the remaining junctions
outDF <- outDF %>% unite("Name", Chrom:End) %>% tibble::column_to_rownames(var = "Name") %>% filter(rowSums(is.na(.))/ncol(.) <= 1/3) %>% 
    filter(apply(., 1, max, na.rm = TRUE) >= 0.05) %>% filter(apply(., 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) >= 0.05) %>%
    t %>% as.data.frame %>% filter(rowSums(is.na(.))/ncol(.) <= 1/3) %>% t %>% as.data.frame

# Perform PCA clustering of samples based on usage frequencies of splice junctions represented in outDF
pcaResults <- pca(t(outDF), method = "nipals", nPcs = 2, scale = "none", center = TRUE)
scoreDF <- as.data.frame(scores(pcaResults)) %>% tibble::rownames_to_column(var = "Sample_ID") %>%
    mutate(Tissue = tissue.assign[Sample_ID]) %>% mutate(Tissue = factor(case_when(Tissue == "Blood" ~ "Whole blood (GTEx v8)",
    Tissue == "Adipose" ~ "Subcutaneous adipose (GTEx v8)", Tissue == "Skin" ~ "Skin - lower leg (GTEx v8)", 
    Tissue == "Fibroblast" ~ "Fibroblasts (GTEx v8)", TRUE ~ "Fibroblasts"), levels = c("Whole blood (GTEx v8)",
    "Subcutaneous adipose (GTEx v8)", "Skin - lower leg (GTEx v8)", "Fibroblasts (GTEx v8)", "Fibroblasts")))
p4 <- ggplot(scoreDF[order(scoreDF$Tissue), ], aes(x = PC1, y = PC2, color = Tissue)) + geom_point(size = 1, stroke = NA) + 
    theme_classic() + xlab(paste("PC1 (", format(round(pcaResults@R2[1]*100, digits = 1), nsmall = 1), "%)", sep = "")) +
    ylab(paste("PC2 (", format(round(pcaResults@R2[2]*100, digits = 1), nsmall = 1), "%)", sep = "")) + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), legend.position = "none", 
    plot.title = element_text(color = "black", size = 7, hjust = 0.5)) + scale_color_manual(values = palette) + coord_cartesian(xlim = c(-3, 3), 
    ylim = c(-3, 3)) + ggtitle("Splice junction usage (PMD-359)")

# Assemble plots into the same grid and save to outfile
p <- plot_grid(p1, p2, p3, p4, nrow = 2, labels = c("a", "b", "c", "d"), label_size = 8)
ggsave(outfile, plot = p, width = 4, height = 4)
