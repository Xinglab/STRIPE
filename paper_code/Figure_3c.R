#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.05
# Figure 3c

# Usage frequencies for splice junctions involving skipping of ATP5MK exons 2 and 3 (top) as well as activation of a 
# cryptic donor splice site within exon 3 (bottom) across cohort samples and tissue-matched GTEx controls.

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
outDF <- tibble(Count = junction.counts, Coverage = junction.coverage) %>% mutate(Group = "Fibroblasts\n(GTEx v8)", Label = "",
    Junction = "ATP5MK cryptic donor splice site (exon 3)")

# Retrieve read support and total coverage for splice junction chr10:103392284-103396408 in GTEx controls
gtex.splice <- full_join(read.table(file.path(workdir, "PMD/references/GTEx_v8/Fibroblast/target_genes.junction_counts.txt"), 
    sep = "\t", header = TRUE, check.names = FALSE) %>% separate(Name, c("Chrom", "Start", "End"), sep = "_") %>% 
    filter(Chrom == "chr10" & (Start == "103392284" | End == "103396408")) %>% distinct(), tibble(Chrom = "chr16", 
    Start = "103392284", End = "103396408"), by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)
junction.counts <- as.numeric(gtex.splice %>% filter(Chrom == "chr16" & Start == "103392284" & End == "103396408") %>% 
    select(-c(Chrom, Start, End)))
junction.coverage <- as.numeric(gtex.splice %>% select(-c(Chrom, Start, End)) %>% colSums)
outDF <- bind_rows(outDF, tibble(Count = junction.counts, Coverage = junction.coverage) %>% mutate(Group = "Fibroblasts\n(GTEx v8)", Label = "", 
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
        sample.labels <- c(sample.labels, "")
    }
}
outDF <- bind_rows(outDF, tibble(Count = junction.counts, Coverage = junction.coverage, Group = "Fibroblasts\n(cohort)", Label = sample.labels, 
    Junction = "ATP5MK cryptic donor splice site (exon 3)"))

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
        sample.labels <- c(sample.labels, "")
    }
}
outDF <- bind_rows(outDF, tibble(Count = junction.counts, Coverage = junction.coverage, Group = "Fibroblasts\n(cohort)", Label = sample.labels, 
    Junction = "ATP5MK exons 2-3 skipping")) %>% mutate(Group = factor(as.character(Group), levels = c("Fibroblasts\n(GTEx v8)", "Fibroblasts\n(cohort)"))) %>% 
    filter(Coverage >= 20) %>% mutate(Usage = Count/Coverage, Junction = factor(as.character(Junction), levels = c("ATP5MK exons 2-3 skipping",
    "ATP5MK cryptic donor splice site (exon 3)")))

p <- ggplot() + geom_jitter(data = outDF %>% filter(Label == ""), aes(x = Group, y = Usage * 100, color = Group), stroke = NA, alpha = 0.5, height = 0, width = 0.2) + 
    geom_point(data = outDF %>% filter(Label != ""), aes(x = Group, y = Usage * 100), stroke = NA, color = "#E43321") + theme_bw() + ylim(0, 100) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", 
    size = 6), axis.title.y = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    legend.position = "none", plot.title = element_text(color = "black", size = 7, hjust = 0.5), axis.title.x = element_blank(), 
    strip.text = element_text(color = "black", size = 6)) + ylab("Splice junction usage (%)") + scale_color_manual(values = c("#1A78AC", "#6FC7CE")) + 
    facet_wrap(vars(Junction), ncol = 1)

ggsave(outfile, plot = p, width = 2, height = 2.75)
