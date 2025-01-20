#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.12.28
# Supplementary Figure 3a

# Comparison of target gene abundances between untargeted long-read RNA-seq and TEQUILA-seq data for fibroblast cell lines CDG-152-1 
# (CDG-466 gene panel, left) and E1877 (PMD-359 gene panel, right). Only target genes with TPM (transcripts per million) > 1 in 
# untargeted data are displayed. 

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
outfile <- file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S3/Figure_S3a.pdf")

# Read in target gene abundances for untargeted long-read RNA-seq and TEQUILA-seq data on fibroblast cell lines 
# CDG-152-1 (CDG-466 gene panel) and E1877 (PMD-359 gene panel)

outDF <- tibble(Gene_ID = character(), TPM = numeric(), Library = character(), Panel = character())
for (disease in c("CDG", "PMD")) {
    sample.id <- ifelse(disease == "CDG", "CDG-152-1", "E1877")
    target.genes <- read.table(file.path(workdir, disease, "references/target_genes.txt"))$V1
    for (library in c("ONT", "TEQUILA")) {
        gene.map <- read.table(file.path(workdir, disease, sample.id, "RNA", ifelse(library == "ONT", "stripe_qc", "stripe"), 
            "quality_control/stringtie", paste(sample.id, "_", library, ".gtf", sep = "")), sep = "\t", header = FALSE) %>% 
            filter(V3 == "transcript") %>% mutate(Gene_ID = unlist(lapply(V9, function(x) PullFeature(x, "gene_id"))), 
            Reference_Gene_ID = unlist(lapply(V9, function(x) PullFeature(x, "ref_gene_id")))) %>% select(Gene_ID, Reference_Gene_ID) %>% drop_na()
        input.data <- read.table(file.path(workdir, disease, sample.id, "RNA", ifelse(library == "ONT", "stripe_qc", "stripe"),
            "quality_control/stringtie", paste(sample.id, "_", library, ".gene_abundance.tsv", sep = "")), sep = "\t", header = TRUE) %>% 
            filter(!(Gene.ID %in% gene.map$Gene_ID)) %>% separate(Gene.ID, c("Gene_ID", NA), sep = "\\.") %>%
            filter(Gene_ID %in% target.genes) %>% select(Gene_ID, TPM) %>% mutate(Library = library, Panel = disease) 
        outDF <- bind_rows(outDF, input.data)
    }
}
outDF$Panel <- factor(outDF$Panel, levels = c("CDG", "PMD"))
levels(outDF$Panel) <- c("CDG-466 gene panel", "PMD-359 gene panel")

# Filter for genes with TPM > 1 in untargeted long-read RNA-seq
outDF <- outDF %>% spread(Library, TPM) %>% filter(ONT > 1)
    
# Compute Spearman's correlation in target gene abundances between untargeted long-read RNA-seq and TEQUILA-seq data
corrDF <- outDF %>% group_by(Panel) %>% summarise(Corr = cor(ONT, TEQUILA, method = "spearman")) %>% mutate(Label = paste("Spearman's r = ", 
    format(round(Corr, digits = 2), nsmall = 2), sep = ""))

# Plot target gene abundances between untargeted long-read RNA-seq and TEQUILA-seq data
p <- ggplot(outDF, aes(x = log2(ONT + 1), y = log2(TEQUILA + 1))) + facet_wrap(. ~ Panel, nrow = 1) + theme_classic() + 
    geom_point(aes(color = Panel), stroke = NA, size = 1, alpha = 0.5) + geom_segment(data = tibble(Panel = factor(c("CDG-466 gene panel", 
    "PMD-359 gene panel"), levels = c("CDG-466 gene panel", "PMD-359 gene panel")), XStart = c(0, 0), XEnd = c(15, 15), YStart = c(0, 0), 
    YEnd = c(15, 15)), aes(x = XStart, xend = XEnd, y = YStart, yend = YEnd), linetype = "dashed", linewidth = 0.25) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.text = element_text(color = "black", size = 6), axis.title = element_text(color = "black", size = 7), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), strip.text = element_text(color = "black", size = 6),
    axis.line = element_line(color = "black", linewidth = 0.25), legend.position = "none") + 
    xlab("Untargeted (cDNA-PCR), log2(TPM+1)") + ylab("TEQUILA-seq, log2(TPM+1)") + geom_text(data = corrDF, mapping = 
    aes(x = 4, y = 2, label = Label), size = 5/14*6, hjust = 0) + scale_color_manual(values = c("#3188BD", "#8AC28E"))

# Save p to outfile
ggsave(outfile, plot = p, width = 3, height = 1.75)