#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.04
# Figure 2d

# (d) Read support and total coverage for a splice junction involving skipping of EARS2 exons 3 to 7 across 
# cohort samples and tissue-matched GTEx controls.

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
outfile <- file.path(workdir, "manuscript/Main_Figures/Figure_2/Figure_2d.pdf")

# Retrieve read support and total coverage for splice junction chr16:23525380-23552148 in GTEx controls
gtex.splice <- full_join(read.table(file.path(workdir, "PMD/references/GTEx_v8/Fibroblast/target_genes.junction_counts.txt"), 
    sep = "\t", header = TRUE, check.names = FALSE) %>% separate(Name, c("Chrom", "Start", "End"), sep = "_") %>% 
    filter(Chrom == "chr16" & (Start == "23525380" | End == "23552148")) %>% distinct(), tibble(Chrom = "chr16", 
    Start = "23525380", End = "23552148"), by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)
junction.counts <- as.numeric(gtex.splice %>% filter(Chrom == "chr16" & Start == "23525380" & End == "23552148") %>% 
    select(-c(Chrom, Start, End)))
junction.coverage <- as.numeric(gtex.splice %>% select(-c(Chrom, Start, End)) %>% colSums)
outDF <- tibble(Count = junction.counts, Coverage = junction.coverage) %>% mutate(Group = "Fibroblasts (GTEx)")

# Retrieve read support and total coverage for splice junction chr16:23525380-23552148 in cohort samples
cohort.samples <- read.table(file.path(workdir, "PMD/samples.txt"), sep = "\t", header = TRUE) %>% filter(Provider == "Rebecca Ganetzky") %>% pull(ID)
sample.labels <- c()
junction.counts <- c()
junction.coverage <- c()

for (sampid in cohort.samples) {
    if (sampid == "BS2-1") {
        inDF <- read.table(file.path(workdir, "PMD", sampid, "RNA/stripe/target_genes/EARS2/aberrant_splicing/output_hap1.txt"), 
            header = TRUE, sep = "\t") %>% filter(Junction == "chr16:23525380-23552148")
        junction.counts <- c(junction.counts, inDF$Count)
        junction.coverage <- c(junction.coverage, inDF$Coverage)
        sample.labels <- c(sample.labels, "BS2-1 (haplotype 1)")
    } else {
        inDF <- read.table(file.path(workdir, "PMD", sampid, "RNA/stripe/quality_control", paste(sampid, "TEQUILA.junctions.tsv", 
            sep = "_")), header = TRUE, sep = "\t")
        inDF <- full_join(inDF %>% filter(Chrom == "chr16" & (Start == 23525380 | End == 23552148)), tibble(Chrom = "chr16", Start = 23525380, 
            End = 23552148), by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)
        junction.counts <- c(junction.counts, inDF %>% filter(Chrom == "chr16" & Start == 23525380 & End == 23552148) %>% pull(Read_Counts))
        junction.coverage <- c(junction.coverage, sum(inDF$Read_Counts))
        sample.labels <- c(sample.labels, "Fibroblasts (cohort)")
    }
}
outDF <- bind_rows(outDF, tibble(Count = junction.counts, Coverage = junction.coverage, Group = sample.labels)) %>% 
    mutate(Group = factor(as.character(Group), levels = c("Fibroblasts (GTEx)", "Fibroblasts (cohort)", "BS2-1 (haplotype 1)")))

# Plot read support versus total coverage for splice junction chr16:23525380-23552148 in cohort samples and GTEx controls
p <- ggplot(outDF, aes(x = Coverage, y = Count, color = Group)) + geom_point(stroke = NA, alpha = 0.8) + theme_bw() + 
    scale_x_continuous(trans = "pseudo_log", breaks = c(0, 10, 100, 1000), limits = c(0, 2000)) + 
    scale_y_continuous(trans = "pseudo_log", breaks = c(0, 10, 100, 1000), limits = c(0, 2000)) + 
    theme(panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    legend.position = "none", plot.title = element_text(color = "black", size = 7, hjust = 0.5)) + ylab("Junction read count") + 
    xlab("Total junction coverage") + scale_color_manual(values = c("#1A78AC", "#6FC7CE", "#E43321")) + ggtitle("EARS2 exons 3-7 skipping")

ggsave(outfile, plot = p, width = 2.5, height = 2)
