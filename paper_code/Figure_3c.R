#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.05
# Figure 3c

# (c) Read support and total coverage for splice junctions involving skipping of ATP5MK exons 2 and 3 (top) as well as
# activation of a cryptic donor splice site within exon 3 (bottom) across cohort samples and tissue-matched GTEx controls.

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
outfile <- file.path(workdir, "manuscript/Main_Figures/Figure_3/Figure_3c.pdf")

# Retrieve read support and total coverage for splice junction chr10:103392284-103392381 in GTEx controls
gtex.splice <- full_join(read.table(file.path(workdir, "PMD/references/GTEx_v8/Fibroblast/target_genes.junction_counts.txt"), 
    sep = "\t", header = TRUE, check.names = FALSE) %>% separate(Name, c("Chrom", "Start", "End"), sep = "_") %>% 
    filter(Chrom == "chr10" & (Start == "103392284" | End == "103392381")) %>% distinct(), tibble(Chrom = "chr16", 
    Start = "103392284", End = "103392381"), by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)
junction.counts <- as.numeric(gtex.splice %>% filter(Chrom == "chr16" & Start == "103392284" & End == "103392381") %>% 
    select(-c(Chrom, Start, End)))
junction.coverage <- as.numeric(gtex.splice %>% select(-c(Chrom, Start, End)) %>% colSums)
outDF <- tibble(Count = junction.counts, Coverage = junction.coverage) %>% mutate(Group = "Fibroblasts (GTEx)", 
    Junction = "ATP5MK cryptic donor splice site (exon 3)")

# Retrieve read support and total coverage for splice junction chr10:103392284-103396408 in GTEx controls
gtex.splice <- full_join(read.table(file.path(workdir, "PMD/references/GTEx_v8/Fibroblast/target_genes.junction_counts.txt"), 
    sep = "\t", header = TRUE, check.names = FALSE) %>% separate(Name, c("Chrom", "Start", "End"), sep = "_") %>% 
    filter(Chrom == "chr10" & (Start == "103392284" | End == "103396408")) %>% distinct(), tibble(Chrom = "chr16", 
    Start = "103392284", End = "103396408"), by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)
junction.counts <- as.numeric(gtex.splice %>% filter(Chrom == "chr16" & Start == "103392284" & End == "103396408") %>% 
    select(-c(Chrom, Start, End)))
junction.coverage <- as.numeric(gtex.splice %>% select(-c(Chrom, Start, End)) %>% colSums)
outDF <- bind_rows(outDF, tibble(Count = junction.counts, Coverage = junction.coverage) %>% mutate(Group = "Fibroblasts (GTEx)", 
    Junction = "ATP5MK exons 2-3 skipping"))

# Retrieve read support and total coverage for splice junction chr10:103392284-103392381 in cohort samples
cohort.samples <- read.table(file.path(workdir, "PMD/samples.txt"), sep = "\t", header = TRUE) %>% filter(Provider == "Rebecca Ganetzky") %>% pull(ID)
sample.labels <- c()
junction.counts <- c()
junction.coverage <- c()

for (sampid in cohort.samples) {
    if (sampid == "Q1663") {
        inDF <- read.table(file.path(workdir, "PMD", sampid, "RNA/stripe/target_genes/ATP5MK/aberrant_splicing/output_all.txt"), 
            header = TRUE, sep = "\t") %>% filter(Junction == "chr10:103392284-103392381")
        junction.counts <- c(junction.counts, inDF$Count)
        junction.coverage <- c(junction.coverage, inDF$Coverage)
        sample.labels <- c(sample.labels, "Q1663")
    } else {
        inDF <- read.table(file.path(workdir, "PMD", sampid, "RNA/stripe/quality_control", paste(sampid, "TEQUILA.junctions.tsv", 
            sep = "_")), header = TRUE, sep = "\t")
        inDF <- full_join(inDF %>% filter(Chrom == "chr10" & (Start == 103392284 | End == 103392381)), tibble(Chrom = "chr16", Start = 103392284, 
            End = 103392381), by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)
        junction.counts <- c(junction.counts, inDF %>% filter(Chrom == "chr16" & Start == 103392284 & End == 103392381) %>% pull(Read_Counts))
        junction.coverage <- c(junction.coverage, sum(inDF$Read_Counts))
        sample.labels <- c(sample.labels, "Fibroblasts (cohort)")
    }
}
outDF <- bind_rows(outDF, tibble(Count = junction.counts, Coverage = junction.coverage, Group = sample.labels) %>% 
    mutate(Junction = "ATP5MK cryptic donor splice site (exon 3)")) 

# Retrieve read support and total coverage for splice junction chr10:103392284-103396408 in cohort samples
sample.labels <- c()
junction.counts <- c()
junction.coverage <- c()

for (sampid in cohort.samples) {
    if (sampid == "Q1663") {
        inDF <- read.table(file.path(workdir, "PMD", sampid, "RNA/stripe/target_genes/ATP5MK/aberrant_splicing/output_all.txt"), 
            header = TRUE, sep = "\t") %>% filter(Junction == "chr10:103392284-103396408")
        junction.counts <- c(junction.counts, inDF$Count)
        junction.coverage <- c(junction.coverage, inDF$Coverage)
        sample.labels <- c(sample.labels, "Q1663")
    } else {
        inDF <- read.table(file.path(workdir, "PMD", sampid, "RNA/stripe/quality_control", paste(sampid, "TEQUILA.junctions.tsv", 
            sep = "_")), header = TRUE, sep = "\t")
        inDF <- full_join(inDF %>% filter(Chrom == "chr10" & (Start == 103392284 | End == 103396408)), tibble(Chrom = "chr16", Start = 103392284, 
            End = 103396408), by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)
        junction.counts <- c(junction.counts, inDF %>% filter(Chrom == "chr16" & Start == 103392284 & End == 103396408) %>% pull(Read_Counts))
        junction.coverage <- c(junction.coverage, sum(inDF$Read_Counts))
        sample.labels <- c(sample.labels, "Fibroblasts (cohort)")
    }
}
outDF <- bind_rows(outDF, tibble(Count = junction.counts, Coverage = junction.coverage, Group = sample.labels) %>% 
    mutate(Junction = "ATP5MK exons 2-3 skipping")) %>% mutate(Group = factor(as.character(Group), levels = c("Fibroblasts (GTEx)", 
    "Fibroblasts (cohort)", "Q1663")), Junction = factor(as.character(Junction), levels = c("ATP5MK exons 2-3 skipping",
    "ATP5MK cryptic donor splice site (exon 3)")))

# Plot read support versus total coverage for splice junctions chr10:103392284-103392381 and chr10:103392284-103396408 in cohort samples and GTEx controls
p <- ggplot(outDF, aes(x = Coverage, y = Count, color = Group)) + facet_wrap(vars(Junction), ncol = 1) + geom_point(stroke = NA, alpha = 0.8) + theme_bw() + 
    scale_x_continuous(trans = "pseudo_log", breaks = c(0, 10, 100, 1000, 10000), limits = c(0, 50000)) + 
    scale_y_continuous(trans = "pseudo_log", breaks = c(0, 10, 100, 1000, 10000), limits = c(0, 50000)) + 
    theme(panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), axis.title = element_text(color = "black", size = 7), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), legend.position = "none", strip.text = element_text(color = "black", size = 6)) + 
    ylab("Junction read count") + xlab("Total junction coverage") + scale_color_manual(values = c("#1A78AC", "#6FC7CE", "#E43321"))

ggsave(outfile, plot = p, width = 2.25, height = 2.5)
