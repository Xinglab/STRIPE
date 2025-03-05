#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.05
# Figure 7d

# (d) Number of uniquely mapped and multi-mapped TEQUILA-seq reads from patient CDG-P22 that are aligned to ALG1 and its 
# known pseudogenes.

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
samtools <- "/scr1/users/wangr5/tools/samtools-1.21/samtools"
gencode <- "/scr1/users/wangr5/references/gencode.v45.annotation.gtf"
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_7/Figure_7d.pdf")

# Pull out all reads supporting ALG1 and its pseudogenes
dir.create(file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_7/tmp"), showWarnings = FALSE)
gene.list <- setNames(c("ALG1", "TSSC2", "ALG1L6P", "ALG1L3P", "ALG1L5P", "ALG1L9P", "ALG1L12P", 
    "ALG1L7P", "ALG1L8P"), c("ENSG00000033011", "ENSG00000223756", "ENSG00000238278", 
    "ENSG00000251087", "ENSG00000226943", "ENSG00000254978", "ENSG00000250794", 
    "ENSG00000251271", "ENSG00000227620"))

# Read in mapping file describing new IDs for cohort samples
sample.map <- read.table(file.path(workdir, "manuscript/Revisions/20250228/sample.map"), sep = "\t", header = FALSE) %>% 
    select(V2, V1) %>% tibble::deframe()

outDF <- tibble(Gene = character(), UniqueCount = integer(), MultiCount = integer())
for (i in 1:length(gene.list)) {
    # Pull out feature annotations and genomic region for gene
    system(paste("grep", names(gene.list[i]), gencode, ">", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_7/tmp/input.gtf")))
    gene.coord <- read.table(file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_7/tmp/input.gtf"), sep = "\t",
        header = FALSE, quote = "") %>% filter(V3 == "gene") %>% reframe(Coord = paste(V1, ":", V4, 
        "-", V5, sep = "")) %>% pull(Coord)
    
    # Pull out TEQUILA-seq reads mapping to gene region
    system(paste(samtools, "view -hb", file.path(workdir, "CDG", sample.map["CDG-P22"], "RNA", paste(sample.map["CDG-P22"], "TEQUILA.bam", sep = "_")), gene.coord,
        ">", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_7/tmp/input.bam")))
    system(paste(samtools, "flagstat -O tsv", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_7/tmp/input.bam"), ">",
        file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_7/tmp/flagstat.txt")))
    
    # Retrieve number of uniquely mapped reads and multi-mapped reads from flagstat.txt
    # Assume that primary mapped reads are uniquely mapped whereas secondary/supplementary reads are multi-mapped
    flagDF <- read.table(file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_7/tmp/flagstat.txt"), sep = "\t", header = FALSE)
    uniqueCount <- as.integer(flagDF[9,1])
    multiCount <- as.integer(flagDF[7,1]) - uniqueCount
    outDF <- bind_rows(outDF, tibble(Gene = gene.list[i], UniqueCount = uniqueCount, MultiCount = multiCount))
}

outDF <- outDF %>% gather("Category", "Count", -Gene) %>% mutate(Category = factor(as.character(Category), levels = c("MultiCount", "UniqueCount")))
levels(outDF$Category) <- c("Mapped to multiple loci", "Uniquely mapped")
geneOrder <- outDF %>% select(-Category) %>% group_by(Gene) %>% summarise(Count = sum(Count)) %>% ungroup %>% arrange(Count) %>% pull(Gene)
outDF$Gene <- factor(as.character(outDF$Gene), levels = geneOrder)

p <- ggplot(outDF, aes(x = Count, y = Gene, fill = Category, group = Category)) + geom_bar(position = "stack", width = 0.75, stat = "identity") + 
    theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.text.x = element_text(color = "black", size = 6), axis.ticks = element_line(color = "black", linewidth = 0.5), 
    axis.text.y = element_text(color = "black", size = 6, face = "italic"), legend.position = c(0.6, 0.4),
    axis.title.x = element_text(color = "black", size = 7), axis.title.y = element_blank(), legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 6), legend.background = element_rect(color = "black", linewidth = 0.5, fill = NA),
    legend.key.size = unit(0.3, "cm")) + xlab("Number of TEQUILA-seq reads (CDG-P22)") + scale_x_continuous(limits = c(0, 20000)) +
    geom_text(data = outDF %>% select(-Category) %>% group_by(Gene) %>% summarise(Count = sum(Count)) %>% ungroup, aes(x = Count, y = Gene, label = Count), hjust = -0.3, 
    size = 5/14*6, inherit.aes = FALSE) + scale_fill_manual(values = c("#E69F00", "#0072B2"))

ggsave(outfile, plot = p, width = 2.25, height = 1.5)

# Remove intermediate files and folders
system(paste("rm", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_7/tmp/*")))
system(paste("rmdir", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_7/tmp")))
