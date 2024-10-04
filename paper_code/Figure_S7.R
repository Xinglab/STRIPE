#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.04
# Supplementary Figure 7

# (a) Haplotype-resolved expression of TRMU in individual Q1642, who is heterozygous for a 1.8-kb structural deletion
# that removes the promoter region and first exon. (b) Major and minor haplotype read counts for TRMU across cohort samples
# and tissue-matched GTEx controls. 

# =====================================================================================================================
#                                                      LIBRARIES 
# =====================================================================================================================

# Load required libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(cowplot))

# =====================================================================================================================
#                                                   HELPER FUNCTIONS
# =====================================================================================================================

FlatternCoverage <- function(fp) {
    # Function designed to report discrete genomic intervals with the same read coverage
    input.file <- file(fp, "r")
    ref.line <- unlist(strsplit(readLines(input.file, 1), "\t"))
    outDF <- tibble(Chrom = character(), Start = integer(), End = integer(), Coverage = integer())

    while (TRUE) {
        curr.line <- unlist(strsplit(readLines(input.file, 1), "\t"))
        if (length(curr.line) == 0) {
            break
        }
        if (curr.line[3] != ref.line[3]) {
            outDF <- bind_rows(outDF, tibble(Chrom = ref.line[1], Start = as.integer(ref.line[2]) - 1,
                End = as.integer(curr.line[2]) - 1, Coverage = as.integer(ref.line[3])))
            ref.line <- curr.line
        } 
    }
    close(input.file)

    return(outDF)
}

RescaleFeature <- function(x) {
    return(ifelse(x < 50 & x > 0, 50, ifelse(x > 300, 300 + round(sqrt(x)), x)))
}

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir <- "/mnt/isilon/lin_lab_share/STRIPE"
target.gene <- read.table(file.path(workdir, "/PMD/references/target_genes.bed"), sep = "\t", header = FALSE) %>% filter(V5 == "TRMU")
gencode.gtf <- "/scr1/users/wangr5/references/gencode.v45.annotation.gtf"
samtools <- "/scr1/users/wangr5/tools/samtools-1.21/samtools"
outfile <- file.path(workdir, "manuscript/Supplementary_Figures/Figure_S7/Figure_S7.pdf")

# =====================================================================================================================
#                                                      PANEL A
# =====================================================================================================================

# Create a temporary directory in the folder for Supplementary Figure 7
dir.create(file.path(dirname(outfile), "tmp"), showWarnings = FALSE)

# Extract canonical transcript annotations for target.gene from gencode.gtf
system(paste("grep \"", target.gene$V4, "\" ", gencode.gtf, " | grep \"Ensembl_canonical\" > ", file.path(dirname(outfile), 
    "tmp/gene.gtf"), sep = ""))

# Compute per-base coverage in target gene from patient haplotype-specific BAM files
system(paste(samtools, "depth -r", paste(target.gene$V1, ":", target.gene$V2+1, "-", target.gene$V3, sep = ""),
    file.path(workdir, "PMD/Q1642/RNA/stripe/target_genes/TRMU/hap1_reads.bam"), ">",
    file.path(dirname(outfile), "tmp/hap1_coverage.txt")))
system(paste(samtools, "depth -r", paste(target.gene$V1, ":", target.gene$V2+1, "-", target.gene$V3, sep = ""),
    file.path(workdir, "PMD/Q1642/RNA/stripe/target_genes/TRMU/hap2_reads.bam"), ">",
    file.path(dirname(outfile), "tmp/hap2_coverage.txt")))

hap1DF <- FlatternCoverage(file.path(dirname(outfile), "tmp/hap1_coverage.txt"))
hap2DF <- FlatternCoverage(file.path(dirname(outfile), "tmp/hap2_coverage.txt"))
gtfDF <- read.table(file.path(dirname(outfile), "tmp/gene.gtf"), sep = "\t", header = FALSE, quote = "") %>% 
    filter(V3 %in% c("exon", "CDS", "UTR"))
gtfDF <- mutate(gtfDF, V3 = case_when(!("CDS" %in% gtfDF$V3) & V3 == "exon" ~ "UTR", TRUE ~ V3), V4 = V4 - 1) %>% 
    filter(V3 %in% c("UTR", "CDS")) %>% arrange(V4)

# Shrink genomic features in gtfDF, hap1DF, and hap2DF
oldCoord <- sort(unique(c(gtfDF$V4, gtfDF$V5, hap1DF$Start, hap1DF$End, hap2DF$Start, hap2DF$End)))
newCoord <- cumsum(c(0, RescaleFeature(tail(oldCoord, -1) - head(oldCoord, -1))))
gtfDF <- mutate(gtfDF, V4 = recode(V4, !!!setNames(newCoord, oldCoord))/max(newCoord),
    V5 = recode(V5, !!!setNames(newCoord, oldCoord))/max(newCoord))
hap1DF <- mutate(hap1DF, Start = recode(Start, !!!setNames(newCoord, oldCoord))/max(newCoord),
    End = recode(End, !!!setNames(newCoord, oldCoord))/max(newCoord))
hap2DF <- mutate(hap2DF, Start = recode(Start, !!!setNames(newCoord, oldCoord))/max(newCoord),
    End = recode(End, !!!setNames(newCoord, oldCoord))/max(newCoord))
intronDF <- data.frame(V3 = rep("intron", nrow(gtfDF) - 1), V4 = head(gtfDF$V5, -1), V5 = tail(gtfDF$V4, -1)) %>% filter(V5 > V4)
utrDF <- filter(gtfDF, V3 == "UTR")
cdsDF <- filter(gtfDF, V3 == "CDS")

gene.coverage <- ggplot(bind_rows(tibble(Coord = c(head(hap1DF$Start, 1), rep(tail(hap1DF$Start, -1), each = 2), tail(hap1DF$End, 1)), 
    Coverage = rep(hap1DF$Coverage, each = 2)) %>% mutate(Group = "Q1642 (haplotype 1)"), tibble(Coord = c(head(hap2DF$Start, 1), 
    rep(tail(hap2DF$Start, -1), each = 2), tail(hap2DF$End, 1)), Coverage = rep(hap2DF$Coverage, each = 2)) %>% mutate(Group = "Q1642 (haplotype 2)")) %>%
    mutate(Group = factor(Group, levels = c("Q1642 (haplotype 1)", "Q1642 (haplotype 2)"))), aes(x = Coord, y = Coverage, color = Group)) + 
    geom_path(alpha = 0.75, linewidth = 0.5) + theme_classic() + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
    axis.text.x = element_blank(), axis.line.x = element_blank(), axis.title.y = element_text(color = "black", size = 7), 
    axis.ticks.y = element_line(color = "black", linewidth = 0.5), axis.text.y = element_text(color = "black", size = 6),
    axis.line.y = element_line(color = "black", linewidth = 0.5), legend.title = element_blank(), legend.position = c(0.9, 0.9),
    legend.text = element_text(color = "black", size = 6)) + ylab("TEQUILA-seq read coverage") + scale_color_manual(values = c("#59AAD1", "#DF5C79")) + ylim(0, 1200)
gene.structure <- ggplot() + geom_rect(data = utrDF, fill = "white", xmin = 0.75, xmax = 1.25, ymin = utrDF$V4, ymax = utrDF$V5, 
    color = "black", linewidth = 0.5) + geom_rect(data = cdsDF, fill = "#BFBFBF", xmin = 0.75, xmax = 1.25, ymin = cdsDF$V4,
    ymax = cdsDF$V5, color = "black", linewidth = 0.5) + geom_segment(data = intronDF, x = 1, xend = 1, y = intronDF$V4, yend = intronDF$V5,
    linewidth = 0.5, color = "black") + theme_classic() + xlim(0.5, 1.5) + theme(axis.title.y = element_blank(), axis.ticks = element_blank(),
    axis.text = element_blank(), axis.line = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.title.x = element_text(color = "black", size = 7)) + coord_flip() + ylab("TRMU (NM_018006.5)") 
p1 <- plot_grid(gene.coverage, gene.structure, ncol = 1, align = "v", rel_heights = c(2.25, 1))

# Remove intermediate files
system(paste("rm -rf", file.path(dirname(outfile), "tmp")))

# =====================================================================================================================
#                                                      PANEL B
# =====================================================================================================================

# Retrieve minor haplotype expression ratios for target gene across GTEx controls
gtex.haplotype <- tail(as.character((read.table(file.path(workdir, "PMD/references/GTEx_v8/Fibroblast/target_genes.haplotype_expression.txt"), 
    sep = "\t", header = TRUE, check.names = FALSE) %>% filter(grepl(target.gene$V4, name)))[1,]), -1)
minor.counts <- unlist(lapply(strsplit(gtex.haplotype, "\\|"), function(x) min(as.integer(x[1]), as.integer(x[2]))))
major.counts <- unlist(lapply(strsplit(gtex.haplotype, "\\|"), function(x) max(as.integer(x[1]), as.integer(x[2]))))
outDF <- tibble(Minor = minor.counts, Major = major.counts) %>% mutate(Group = "Fibroblasts (GTEx)")

# Retrieve minor haplotype expression ratios for target gene across cohort samples
cohort.samples <- read.table(file.path(workdir, "PMD/samples.txt"), sep = "\t", header = TRUE) %>% filter(Provider == "Rebecca Ganetzky") %>% pull(ID)
sample.labels <- c()
minor.counts <- c()
major.counts <- c()
for (sampid in cohort.samples) {
    infile <- file.path(workdir, "PMD", sampid, "RNA/stripe/target_genes/TRMU/phasing_stats.txt")
    if (file.exists(infile)) {
        if (sampid == "Q1642") {
            sample.labels <- c(sample.labels, sampid)
        } else {
            sample.labels <- c(sample.labels, "Fibroblasts (cohort)")
        }
        inDF <- read.table(infile, header = FALSE, sep = "\t")
        minor.counts <- c(minor.counts, min(inDF[2,2], inDF[3,2]))
        major.counts <- c(major.counts, max(inDF[2,2], inDF[3,2]))
    }
}
outDF <- bind_rows(outDF, tibble(Minor = minor.counts, Major = major.counts, Group = sample.labels)) %>% mutate(Group = factor(as.character(Group), 
    levels = c("Fibroblasts (GTEx)", "Fibroblasts (cohort)", "Q1642")))

p2 <- ggplot(outDF, aes(x = Major, y = Minor, color = Group)) + geom_point(stroke = NA, alpha = 0.8) + theme_bw() + 
    scale_x_continuous(trans = "pseudo_log", breaks = c(0, 10, 100, 1000)) + scale_y_continuous(trans = "pseudo_log", 
    breaks = c(0, 10, 100, 1000)) + theme(panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    legend.position = "none", plot.title = element_text(color = "black", size = 7, hjust = 0.5)) + ylab("Minor haplotype read count") + 
    xlab("Major haplotype read count") + scale_color_manual(values = c("#1A78AC", "#6FC7CE", "#E43321")) + ggtitle("TRMU")

# Assemble p1 and p2 onto the same plotting grid
p <- plot_grid(p1, p2, nrow = 1, rel_widths = c(2, 1), labels = c("a", "b"), label_size = 8)
ggsave(outfile, plot = p, width = 6.5, height = 2.5)
