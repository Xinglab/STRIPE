#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.12.27
# Supplementary Figure 1

# Cumulative distribution of normalized abundances (transcripts per million, TPM) for 183 genes known to be implicated in congenital disorders of glycosylation (CDG), 
# and 292 genes known to be implicated in primary mitochondrial diseases (PMDs) among 504 fibroblast (blue) and 755 whole blood (red) GTEx samples. 

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
outfile <- file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S1/Figure_S1.pdf")

# Gather median TPMs for known CDG genes and known PMD genes among GTEx v8 whole blood and fibroblast samples
outDF <- tibble(Gene_ID = character(), Median = numeric(), Tissue = character(), Disease = character())
for (disease in c("CDG", "PMD")) {
    input.genes <- read.table(file.path(workdir, disease, "references/known_disease_genes.txt"))$V1
    for (tissue in c("Fibroblast", "Blood")) {
        outDF <- bind_rows(outDF, read.table(file.path(workdir, disease, "references/GTEx_v8", tissue, "target_genes.gene_tpm.txt"), sep = "\t", 
            header = TRUE, stringsAsFactors = FALSE, check.names = FALSE) %>% separate(Name, c("Gene_ID", NA), sep = "\\.") %>% 
            filter(Gene_ID %in% input.genes) %>% rowwise(Gene_ID) %>% mutate(Median = median(c_across(where(is.numeric)))) %>% 
            select(Gene_ID, Median) %>% mutate(Tissue = tissue, Disease = disease))
        
    }
}
outDF$Tissue <- factor(outDF$Tissue, levels = c("Fibroblast", "Blood"))
levels(outDF$Tissue) <- c("Fibroblasts (n = 504)", "Whole blood (n = 755)")
outDF$Disease <- factor(outDF$Disease, levels = c("CDG", "PMD"))
levels(outDF$Disease) <- c("183 known CDG genes", "292 known PMD genes")

p <- ggplot(outDF, aes(x = Median, color = Tissue, group = Tissue)) + facet_wrap(. ~ Disease, nrow = 1) + stat_ecdf(geom = "step", linewidth = 0.5) + 
    theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 6), axis.title = element_text(color = "black", size = 7), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), axis.line = element_line(color = "black", linewidth = 0.25),
    strip.text = element_text(color = "black", size = 6), legend.position = "right", legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 6), legend.key.height = unit(3, "mm")) + xlab("Median abundance (TPM)") + 
    ylab("Cumulative fraction") + scale_x_continuous(trans = scales::pseudo_log_trans(base = 10), limits = c(0, 1000), breaks = c(0, 10, 100, 1000)) + 
    guides(color = guide_legend(ncol = 1)) + geom_vline(xintercept = 10, linewidth = 0.25, linetype = "dashed") +
    scale_color_manual(values = c("#5A8CA8", "#D6604D"))

# Save p to outfile
ggsave(outfile, plot = p, width = 4.5, height = 1.5)