#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.02.28
# Figure 1e

# Comparison of usage frequencies for annotated splice junctions in target genes between untargeted long-read RNA-seq and 
# TEQUILA-seq data for fibroblast cell lines CDG-P12 (CDG-466 gene panel, left) and PMD-C01 (PMD-359 gene panel, right). 

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

ComputePSI <- function(countDF) {
    # Function designed to compute usage frequencies for splice junctions in countDF
    startCoverage <- select(countDF, !(End)) %>% group_by(Chrom, Start) %>% mutate_if(is.numeric, sum) %>% ungroup %>% select(-c(Chrom, Start))
    endCoverage <- select(countDF, !(Start)) %>% group_by(Chrom, End) %>% mutate_if(is.numeric, sum) %>% ungroup %>% select(-c(Chrom, End))
    totalCoverage <- startCoverage + endCoverage - select(countDF, -c(Chrom, Start, End))
    return(bind_cols(countDF %>% select(c(Chrom, Start, End)), select(countDF, -c(Chrom, Start, End)) / mutate(totalCoverage, across(everything(), 
        ~ ifelse(. < 20, NA, .)))))
}

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir <- "/mnt/isilon/lin_lab_share/STRIPE"
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_1/Figure_1e.pdf")

# Read in mapping file describing new IDs for cohort samples
sample.map <- read.table(file.path(workdir, "manuscript/Revisions/20250228/sample.map"), sep = "\t", header = FALSE) %>% 
    select(V2, V1) %>% tibble::deframe()

outDF <- tibble(Chrom = character(), Start = integer(), End = integer(), ONT = numeric(), TEQUILA = numeric(),
    Panel = character())
for (disease in c("CDG", "PMD")) {
    # Read in GENCODE v45 annotations as a dataframe and only retain exon-level annotations
    system(paste("grep -f", file.path(workdir, disease, "references/target_genes.txt"), 
        "/scr1/users/wangr5/references/gencode.v45.annotation.gtf", "| awk '$3 == \"exon\"' >",
        file.path(workdir, "tmp/input.gtf")))
    gencode <- read.table(file.path(workdir, "tmp/input.gtf"), sep = "\t", header = FALSE) %>%
        mutate(Transcript_ID = unlist(lapply(V9, function(x) PullFeature(x, "transcript_id")))) %>%
        select(V1, V4, V5, Transcript_ID) %>% group_by(V1, Transcript_ID) %>% 
        summarise(V4 = paste(sort(V4), collapse = ","), V5 = paste(sort(V5), collapse = ","))
    
    # Retrieve annotated junctions
    junctions <- c()
    for (i in 1:nrow(gencode)) {
        exonStarts <- as.numeric(unlist(strsplit(gencode[[i,3]], ",")))
        exonEnds <- as.numeric(unlist(strsplit(gencode[[i,4]], ",")))
        if (length(exonStarts) > 0) {
            for (j in 1:(length(exonStarts) - 1)) {
                junctions <- c(junctions, paste(gencode[[i,1]], exonEnds[j]+1, exonStarts[j+1]-1, sep = "_"))
            }
        }
    }
    junctions <- unique(junctions)

    # Calculate usage frequencies for annotated junctions from target genes in untargeted long-read RNA-seq
    # and TEQUILA-seq data
    sample.id <- ifelse(disease == "CDG", sample.map["CDG-P12"], sample.map["PMD-C01"])
    currDF <- tibble(Chrom = character(), Start = integer(), End = integer())
    for (library in c("ONT", "TEQUILA")) {
        currDF <- full_join(currDF, read.table(file.path(workdir, disease, sample.id, "RNA", ifelse(library == "ONT", 
            "stripe_qc", "stripe"), "quality_control", paste(sample.id, "_", library, ".junctions.tsv", sep = "")), sep = "\t",
            header = TRUE) %>% filter(paste(Chrom, Start, End, sep = "_") %in% junctions) %>%
            setNames(c("Chrom", "Start", "End", library)), by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)
    }
    outDF <- bind_rows(outDF, ComputePSI(currDF) %>% drop_na %>% mutate(Panel = disease))

    # Remove intermediate file
    system(paste("rm", file.path(workdir, "tmp/input.gtf")))
}
outDF$Panel <- factor(outDF$Panel, levels = c("CDG", "PMD"))
levels(outDF$Panel) <- c("CDG-466 gene panel", "PMD-359 gene panel")

# Compute Spearman's correlation in junction usage frequencies between untargeted long-read RNA-seq and TEQUILA-seq
corrDF <- outDF %>% group_by(Panel) %>% summarise(Corr = cor(ONT, TEQUILA, method = "spearman")) %>%
    mutate(Label = paste("Spearman's r = ", format(round(Corr, digits = 2), nsmall = 2), sep = ""))

# Plot usage frequencies for annotated splice junctions mapping to target genes between untargeted long-read RNA-seq and TEQUILA-seq
p <- ggplot(outDF, aes(x = ONT * 100, y = TEQUILA * 100)) + facet_wrap(. ~ Panel, nrow = 1) + theme_classic() + 
    geom_point(aes(color = Panel), stroke = NA, size = 1, alpha = 0.5) + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    strip.text = element_text(color = "black", size = 6), axis.line = element_line(color = "black", linewidth = 0.25), legend.position = "none") + 
    xlab("Untargeted (cDNA-PCR) (%)") + ylab("TEQUILA-seq (%)") + geom_text(data = corrDF, mapping = 
    aes(x = 35, y = 10, label = Label), size = 5/14*6, hjust = 0) + scale_color_manual(values = c("#3188BD", "#8AC28E"))

# Save p to outfile
ggsave(outfile, plot = p, width = 3, height = 1.75)
