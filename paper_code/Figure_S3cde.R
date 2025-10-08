#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.10.07
# Supplementary Figure 3c,d,e

# (c, d) Principal component analysis (PCA) based on abundance ranks of target genes from the CDG-466 gene panel (left) and the 
# PMD-359 gene panel (right) in fibroblasts from 88 cohort individuals and GTEx samples for clinically accessible tissues
# focusing on (c) PC3 and PC4 and (d) PC5 and PC6. (e) Percentage of variance explained by the first ten principal components 
# for abundance ranks of target genes from the CDG-466 gene panel (left) and the PMD-359 gene panel (right). 

# =====================================================================================================================
#                                                      LIBRARIES 
# =====================================================================================================================

# Load required libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(pcaMethods))

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir <- "/mnt/isilon/lin_lab_share/STRIPE"
outdir <- file.path(workdir, "manuscript/Revisions/20250228/Supplementary_Figures/Figure_S3")

# Read in mapping file describing new IDs for cohort samples
sample.map <- read.table(file.path(workdir, "manuscript/Revisions/20250228/sample.map"), sep = "\t", header = FALSE) %>% 
    select(V2, V1) %>% tibble::deframe()

plotDF <- tibble(Panel = character(), Tissue = character(), PC1 = numeric(), PC2 = numeric(), PC3 = numeric(),
    PC4 = numeric(), PC5 = numeric(), PC6 = numeric(), PC7 = numeric(), PC8 = numeric(), PC9 = numeric(), PC10 = numeric())
varDF <- tibble(Panel = character(), PC = integer(), Variance = numeric())

for (disease in c("CDG", "PMD")) {
    target.genes <- read.table(file.path(workdir, disease, "references/target_genes.txt"))$V1
    outDF <- tibble(Gene_ID = target.genes)
    if (disease == "CDG") {
        cohort.samples <- read.table(file.path(workdir, disease, "samples.txt"), sep = "\t", header = TRUE) %>%
            filter(Provider != "Lan Lin" & ID != sample.map["PMD-C01"]) %>% pull(ID)
    } else {
        cohort.samples <- read.table(file.path(workdir, disease, "samples.txt"), sep = "\t", header = TRUE) %>%
            filter(Provider == "Rebecca Ganetzky" | (Provider == "Marni Falk" & Status != "Diagnosed" & 
            !grepl("-MF", ID))) %>% pull(ID)
    }
    tissue.assign <- setNames(rep("Fibroblasts (cohort)", length(cohort.samples)), cohort.samples)

    # Read in target gene abundances for TEQUILA-seq data on fibroblast cell lines from 88 cohort individuals
    for (sample.id in cohort.samples) {
        sample.data <- read.table(file.path(workdir, disease, sample.id, "RNA/stripe/quality_control/stringtie", 
            paste(sample.id, "TEQUILA.gene_abundance.tsv", sep = "_")), sep = "\t", header = TRUE) %>% 
            filter(grepl("ENSG", Gene.ID)) %>% separate(Gene.ID, c("Gene_ID", NA), sep = "\\.") %>% filter(
            Gene_ID %in% target.genes) %>% select(Gene_ID, TPM) %>% setNames(c("Gene_ID", sample.id))
        outDF <- left_join(outDF, sample.data, by = join_by(Gene_ID)) %>% replace(is.na(.), 0)
    }

    # Read in target gene abundances for GTEx fibroblasts, skin, adipose, and whole blood
    for (tissue in c("Fibroblast", "Skin", "Adipose", "Blood")) {
        tissueDF <- read.table(file.path(workdir, disease, "references/GTEx_v8", tissue, "target_genes.gene_tpm.txt"),
            sep = "\t", header = TRUE, check.names = FALSE) %>% separate(Name, c("Gene_ID", NA), sep = "\\.")
        tissue.assign <- c(tissue.assign, setNames(rep(ifelse(tissue == "Fibroblast", "Fibroblasts (GTEx v8)",
            ifelse(tissue == "Blood", "Whole blood (GTEx v8)", paste(tissue, "(GTEx v8)"))), ncol(tissueDF) - 1), 
            colnames(tissueDF)[-1]))
        outDF <- full_join(outDF, tissueDF, by = join_by(Gene_ID)) %>% replace(is.na(.), 0)
    }

    # Perform PCA clustering of samples based on target gene abundance ranks
    outDF <- outDF %>% mutate_if(is.numeric, function(x) rank(-x, ties.method = "max")) %>%
        tibble::column_to_rownames(var = "Gene_ID")
    pcaResults <- pca(t(outDF), method = "nipals", nPcs = 10, scale = "none", center = TRUE)
    scoreDF <- as.data.frame(scores(pcaResults)) %>% tibble::rownames_to_column(var = "Sample_ID") %>%
        mutate(Tissue = tissue.assign[Sample_ID], Panel = disease, Tissue = factor(Tissue, levels = 
        c("Whole blood (GTEx v8)", "Adipose (GTEx v8)", "Skin (GTEx v8)", "Fibroblasts (GTEx v8)", 
        "Fibroblasts (cohort)")))
    plotDF <- bind_rows(plotDF, scoreDF %>% select(Panel, Tissue, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))
    varDF <- bind_rows(varDF, tibble(PC = 1:10, Variance = attr(pcaResults, "R2")) %>% 
        mutate(Panel = disease) %>% select(Panel, PC, Variance))
}
plotDF$Panel <- factor(plotDF$Panel, levels = c("CDG", "PMD"))
levels(plotDF$Panel) <- c("CDG-466 gene panel", "PMD-359 gene panel")
plotDF$Tissue <- factor(plotDF$Tissue, levels = c("Whole blood (GTEx v8)", "Adipose (GTEx v8)", 
    "Skin (GTEx v8)", "Fibroblasts (GTEx v8)", "Fibroblasts (cohort)"))
varDF$Panel <- factor(varDF$Panel, levels = c("CDG", "PMD"))
levels(varDF$Panel) <- c("CDG-466 gene panel", "PMD-359 gene panel")

# PC3 vs. PC4
p <- ggplot(plotDF[order(plotDF$Tissue), ], aes(x = PC3, y = PC4, color = Tissue)) + facet_wrap(. ~ Panel, nrow = 1) + 
    geom_point(size = 1, stroke = NA) + theme_classic() + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    axis.line = element_line(color = "black", linewidth = 0.25), strip.text = element_text(color = "black", size = 6),
    legend.text = element_text(color = "black", size = 6), legend.title = element_blank(), legend.key.size = unit(0.3, 'cm')) + 
    scale_color_manual(values = c("#D6604D", "#E0C27C", "#689E45", "#5A8BA8", "#94C1B6"))

ggsave(file.path(outdir, "Figure_S3c.pdf"), plot = p, width = 4.5, height = 1.75)

# PC5 vs. PC6
p <- ggplot(plotDF[order(plotDF$Tissue), ], aes(x = PC5, y = PC6, color = Tissue)) + facet_wrap(. ~ Panel, nrow = 1) + 
    geom_point(size = 1, stroke = NA) + theme_classic() + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    axis.line = element_line(color = "black", linewidth = 0.25), strip.text = element_text(color = "black", size = 6),
    legend.text = element_text(color = "black", size = 6), legend.title = element_blank(), legend.key.size = unit(0.3, 'cm')) + 
    scale_color_manual(values = c("#D6604D", "#E0C27C", "#689E45", "#5A8BA8", "#94C1B6"))

ggsave(file.path(outdir, "Figure_S3d.pdf"), plot = p, width = 4.5, height = 1.75)

# Plot percentage of variance explained by individual PCs:
varDF$Label <- sprintf("%.1f", varDF$Variance * 100)

p <- ggplot(varDF, aes(x = factor(PC), y = Variance * 100)) + facet_wrap(. ~ Panel, nrow = 1) + 
    geom_bar(stat = "identity", fill = "#bbd6ff") + theme_classic() + scale_x_discrete(breaks = c(1, 4, 7, 10)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.text = element_text(color = "black", size = 6), axis.title = element_text(color = "black", size = 7), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), axis.line = element_line(color = "black", linewidth = 0.25), 
    strip.text = element_text(color = "black", size = 6)) + ylim(0, 50) + xlab("Principal component index") +
    ylab("% of variance explained") + geom_text(aes(label = Label), vjust = -0.5, hjust = 0.5, color = "black",
    size = 5/14*6)

ggsave(file.path(outdir, "Figure_S3e.pdf"), plot = p, width = 4.5, height = 1.75)