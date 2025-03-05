#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.05
# Figure 7b

# (b) Multiple sequence alignment between the coding sequence of ALG1 (partitioned by exon) and its 8 known pseudogenes 
# (TSSC2, ALG1L6P, ALG1L3P, ALG1L5P, ALG1L9P, ALG1L12P, ALG1L7P, and ALG1L8P). Matches, mismatches, and gaps are 
# represented in blue, red, and gray respectively. 

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
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_7/Figure_7b.pdf")

# Read in multiple sequence alignment file as a dataframe
gene.list <- setNames(c("ALG1", "TSSC2", "ALG1L6P", "ALG1L3P", "ALG1L5P", "ALG1L9P", "ALG1L12P", 
    "ALG1L7P", "ALG1L8P"), c("ENST00000262374", "ENST00000450217", "ENST00000515248", "ENST00000503331",
    "ENST00000482043", "ENST00000532875", "ENST00000515046", "ENST00000502317", "ENST00000533887"))
tx.info <- read.table(file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_7/ALG1/sequences.fa"), sep = "\t", 
    header = FALSE) %>% filter(row_number() %% 2 == 1)
inDF <- bind_cols(data.frame(do.call(rbind, strsplit(read.table(file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_7/ALG1/sequences.fa"), 
    sep = "\t", header = FALSE) %>% filter(row_number() %% 2 == 0) %>% pull(V1), ""))), suppressWarnings(separate(tx.info, V1, "TXID", 
    sep = " ")) %>% mutate(TXID = gsub(">", "", TXID)) %>% separate(TXID, c("TXID", NA), sep = "\\.") %>% 
    mutate(TXID = recode(TXID, !!!gene.list, .default = TXID))) %>% tibble::column_to_rownames(var = "TXID")

# Extract CDS and exon boundary information for ALG1 from tx.info
canonical.tx <- unlist(strsplit(filter(tx.info, grepl("ENST00000262374", V1)) %>% pull(V1), " "))
exon.intervals <- data.frame(Exons = unlist(strsplit(gsub("segs:", "", grep("segs", canonical.tx, value = TRUE)), ","))) %>%
    separate(Exons, c("Start", "End"), sep = "-") %>% mutate_if(is.character, as.integer) %>% tibble::rownames_to_column(var = "Exon")
cds.interval <- as.integer(unlist(strsplit(gsub("CDS=", "", grep("CDS=", canonical.tx, value = TRUE)), "-")))

# Annotate each column of inDF with the exon it corresponds to relative to ALG1 and whether it resides in the CDS
exon.assign <- c()
within.cds <- c()
ctr <- 0
for (align.col in colnames(inDF)) {
    if (inDF["ALG1", align.col] != "-") {
        ctr <- ctr + 1
    }
    exon.assign <- c(exon.assign, filter(exon.intervals, Start <= ctr & End >= ctr) %>% pull(Exon))
    within.cds <- c(within.cds, cds.interval[1] <= ctr & cds.interval[2] >= ctr)
}

# Filter inDF for alignment positions corresponding to the CDS of ALG1
exon.assign <- exon.assign[within.cds]
inDF <- inDF[, within.cds]
inDF <- bind_rows(mutate(inDF["ALG1", ], across(everything(), ~ case_when(. == "-" ~ 0, TRUE ~ 1))),
    inDF[rownames(inDF) != "ALG1", ] %>% mutate(across(everything(), ~ case_when(. == "-" ~ 0,
    . == inDF["ALG1", cur_column()] ~ 1, TRUE ~ -1))))

p <- ggplot(inDF %>% tibble::rownames_to_column(var = "Gene") %>% gather(Position, Status, -Gene) %>%
    mutate(Position = recode(Position, !!!setNames(1:ncol(inDF), colnames(inDF)))) %>%
    mutate(Exon = recode(Position, !!!setNames(exon.assign, 1:ncol(inDF)))) %>% 
    mutate(Gene = factor(Gene, levels = rev(as.character(gene.list))), 
    Status = factor(Status, levels = c(-1, 0, 1)), Exon = factor(Exon, levels = unique(exon.assign))) %>%
    group_by(Exon) %>% mutate(Position = rank(Position)), aes(x = Position,
    y = Gene, fill = Status)) + geom_tile() + theme_classic() + facet_grid(. ~ Exon, scales = "free_x", space = "free_x") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), 
    axis.text.y = element_text(color = "black", size = 6, face = "italic"), legend.position = "none",
    strip.text = element_text(color = "black", size = 6), axis.line = element_blank(), panel.spacing = unit(0.01, "cm"),
    panel.border = element_blank(), strip.background = element_blank()) + scale_fill_manual(values = c("#CC6677", "#DDDDDD", "#88CCEE"))

ggsave(outfile, plot = p, width = 3, height = 1.5)
