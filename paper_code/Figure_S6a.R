#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.01.01
# Supplementary Figure 6a

# Haplotype-resolved read coverage over MPDU1 in patient JaWe, who is heterozygous for a stop-gain variant 
# (NM_004870.4:c.19G>T, p.Gly7Ter). 

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
target.gene <- read.table(file.path(workdir, "/CDG/references/target_genes.bed"), sep = "\t", header = FALSE) %>% filter(V5 == "MPDU1")
gencode.gtf <- "/scr1/users/wangr5/references/gencode.v45.annotation.gtf"
samtools <- "/scr1/users/wangr5/tools/samtools-1.21/samtools"
outfile <- file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S6/Figure_S6a.pdf")

# Create a temporary directory
dir.create(file.path(dirname(outfile), "tmp"), showWarnings = FALSE)

# Extract canonical transcript annotations for target.gene from gencode.gtf
system(paste("grep \"", target.gene$V4, "\" ", gencode.gtf, " | grep \"Ensembl_canonical\" > ", file.path(dirname(outfile), 
    "tmp/gene.gtf"), sep = ""))

# Compute per-base coverage in target gene from patient haplotype-specific BAM files
system(paste(samtools, "depth -r", paste(target.gene$V1, ":", target.gene$V2+1, "-", target.gene$V3, sep = ""),
    file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/JaWe/MPDU1/hap1_reads.bam"),
    ">", file.path(dirname(outfile), "tmp/hap1_coverage.txt")))
system(paste(samtools, "depth -r", paste(target.gene$V1, ":", target.gene$V2+1, "-", target.gene$V3, sep = ""),
    file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/JaWe/MPDU1/hap2_reads.bam"),
    ">", file.path(dirname(outfile), "tmp/hap2_coverage.txt")))

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

plotDF <- bind_rows(tibble(Coord = c(head(hap1DF$Start, 1), rep(tail(hap1DF$Start, -1), each = 2), tail(hap1DF$End, 1)), 
    Coverage = rep(hap1DF$Coverage, each = 2)) %>% mutate(Group = "Haplotype 1"), tibble(Coord = c(head(hap2DF$Start, 1), 
    rep(tail(hap2DF$Start, -1), each = 2), tail(hap2DF$End, 1)), Coverage = rep(hap2DF$Coverage, each = 2)) %>% 
    mutate(Group = "Haplotype 2")) %>% mutate(Group = factor(Group, levels = c("Haplotype 1", "Haplotype 2")))
gene.coverage <- ggplot(plotDF, aes(x = Coord, y = Coverage, color = Group)) + 
    geom_path(alpha = 0.75, linewidth = 0.5) + theme_classic() + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
    axis.text.x = element_blank(), axis.line.x = element_blank(), axis.title.y = element_text(color = "black", size = 7), 
    axis.ticks.y = element_line(color = "black", linewidth = 0.5), axis.text.y = element_text(color = "black", size = 6),
    axis.line.y = element_line(color = "black", linewidth = 0.5), legend.title = element_blank(), legend.position = "top",
    legend.text = element_text(color = "black", size = 6), legend.background = element_blank()) + 
    ylab("Read coverage") + scale_color_manual(values = c("#EC7B02", "#90AFB1")) + ylim(0, 25000) + 
    guides(color = guide_legend(nrow = 1))
gene.structure <- ggplot() + geom_rect(data = utrDF, fill = "white", xmin = 0.75, xmax = 1.25, ymin = utrDF$V4, ymax = utrDF$V5, 
    color = "black", linewidth = 0.5) + geom_rect(data = cdsDF, fill = "#BFBFBF", xmin = 0.75, xmax = 1.25, ymin = cdsDF$V4,
    ymax = cdsDF$V5, color = "black", linewidth = 0.5) + geom_segment(data = intronDF, x = 1, xend = 1, y = intronDF$V4, yend = intronDF$V5,
    linewidth = 0.5, color = "black") + theme_classic() + xlim(0.5, 1.5) + theme(axis.title.y = element_blank(), axis.ticks = element_blank(),
    axis.text = element_blank(), axis.line = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.title.x = element_text(color = "black", size = 7)) + coord_flip() + ylab("MPDU1 (NM_004870.4)") 
p <- cowplot::plot_grid(gene.coverage, gene.structure, ncol = 1, align = "v", rel_heights = c(2.25, 1))

ggsave(outfile, plot = p, width = 3, height = 2)

# Remove intermediate files
system(paste("rm -rf", file.path(dirname(outfile), "tmp")))
