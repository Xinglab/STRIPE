#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.01
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
suppressMessages(library(cowplot))

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir <- "/mnt/isilon/lin_lab_share/STRIPE"
palette <- setNames(c("#0571B0", "#CA0020"), c("Fibroblasts (n = 504)", "Whole blood (n = 755)"))
outfile <- file.path(workdir, "manuscript/Supplementary_Figures/Figure_S1/Figure_S1.pdf")

# Plot median TPM for known disease genes (CDG) among GTEx v8 whole blood and fibroblast samples
cdg.genes <- read.table(file.path(workdir, "CDG/references/known_disease_genes.txt"))$V1
cdg.gtex.blood <- read.table(file.path(workdir, "CDG/references/GTEx_v8/Blood/target_genes.gene_tpm.txt"), sep = "\t", 
    header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
cdg.gtex.fibroblast <- read.table(file.path(workdir, "CDG/references/GTEx_v8/Fibroblast/target_genes.gene_tpm.txt"), sep = "\t", 
    header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
cdg.data <- bind_rows(
    cdg.gtex.blood %>% separate(Name, c("Gene_ID", NA), sep = "\\.") %>% filter(Gene_ID %in% cdg.genes) %>%
    rowwise(Gene_ID) %>% mutate(Median = median(c_across(where(is.numeric)))) %>% select(Gene_ID, Median) %>% 
    mutate(Tissue = "Whole blood (n = 755)"),
    cdg.gtex.fibroblast %>% separate(Name, c("Gene_ID", NA), sep = "\\.") %>% filter(Gene_ID %in% cdg.genes) %>%
    rowwise(Gene_ID) %>% mutate(Median = median(c_across(where(is.numeric)))) %>% select(Gene_ID, Median) %>% 
    mutate(Tissue = "Fibroblasts (n = 504)")
)
p1 <- ggplot(cdg.data, aes(x = Median, color = Tissue, group = Tissue)) + stat_ecdf(geom = "step", linewidth = 0.5) + 
    theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 6), axis.title = element_text(color = "black", size = 7), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), plot.title = element_text(color = "black", size = 7, hjust = 0.5), 
    legend.position = "none") + xlab("Median abundance (TPM)") + ylab("Cumulative fraction") + scale_x_continuous(trans = 
    scales::pseudo_log_trans(base = 10), limits = c(0, 1000), breaks = c(0, 10, 100, 1000)) + scale_color_manual(values = palette) + 
    ggtitle("183 known genes in congenital disorders of glycosylation") + geom_vline(xintercept = 10, linewidth = 0.25, linetype = "dashed")

# Plot median TPM for known disease genes (PMD) among GTEx v8 whole blood and fibroblast samples
pmd.genes <- read.table(file.path(workdir, "PMD/references/known_disease_genes.txt"))$V1
pmd.gtex.blood <- read.table(file.path(workdir, "PMD/references/GTEx_v8/Blood/target_genes.gene_tpm.txt"), sep = "\t", 
    header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
pmd.gtex.fibroblast <- read.table(file.path(workdir, "PMD/references/GTEx_v8/Fibroblast/target_genes.gene_tpm.txt"), sep = "\t", 
    header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
pmd.data <- bind_rows(
    pmd.gtex.blood %>% separate(Name, c("Gene_ID", NA), sep = "\\.") %>% filter(Gene_ID %in% pmd.genes) %>%
    rowwise(Gene_ID) %>% mutate(Median = median(c_across(where(is.numeric)))) %>% select(Gene_ID, Median) %>% 
    mutate(Tissue = "Whole blood (n = 755)"),
    pmd.gtex.fibroblast %>% separate(Name, c("Gene_ID", NA), sep = "\\.") %>% filter(Gene_ID %in% pmd.genes) %>%
    rowwise(Gene_ID) %>% mutate(Median = median(c_across(where(is.numeric)))) %>% select(Gene_ID, Median) %>% 
    mutate(Tissue = "Fibroblasts (n = 504)")
)
p2 <- ggplot(pmd.data, aes(x = Median, color = Tissue, group = Tissue)) + stat_ecdf(geom = "step", linewidth = 0.5) + 
    theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 6), axis.title = element_text(color = "black", size = 7), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), legend.text = element_text(color = "black", size = 6), 
    legend.title = element_blank(), plot.title = element_text(color = "black", size = 7, hjust = 0.5), legend.key.height = unit(3, "mm"), 
    legend.position = "bottom") + xlab("Median abundance (TPM)") + ylab("Cumulative fraction") + scale_x_continuous(trans = 
    scales::pseudo_log_trans(base = 10), limits = c(0, 1000), breaks = c(0, 10, 100, 1000)) + scale_color_manual(values = palette) + 
    ggtitle("292 known genes in primary mitochondrial diseases") + geom_vline(xintercept = 10, linewidth = 0.25, linetype = "dashed")

# Assemble p1 and p2 into the same grid and save to outfile
p <- plot_grid(p1, p2, align = "v", ncol = 1, rel_heights = c(0.8, 1))
ggsave(outfile, plot = p, width = 3, height = 4)
