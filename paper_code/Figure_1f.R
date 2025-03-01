#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.02.28
# Figure 1f

# Principal component analysis (PCA) based on usage frequencies for annotated splice junctions in target genes from the 
# CDG-466 gene panel (left) and the PMD-359 gene panel (right) in fibroblasts from 88 cohort individuals and GTEx samples 
# for clinically accessible tissues.  

# =====================================================================================================================
#                                                      LIBRARIES 
# =====================================================================================================================

# Load required libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(pcaMethods))

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
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_1/Figure_1f.pdf")

# Read in mapping file describing new IDs for cohort samples
sample.map <- read.table(file.path(workdir, "manuscript/Revisions/20250228/sample.map"), sep = "\t", header = FALSE) %>% 
    select(V2, V1) %>% tibble::deframe()

plotDF <- tibble(Panel = character(), Tissue = character(), PC1 = numeric(), PC2 = numeric())
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

    # Calculate usage frequencies for annotated junctions from target genes in TEQUILA-seq data on fibroblast cell
    # lines from 88 cohort individuals
    outDF <- tibble(Chrom = character(), Start = integer(), End = integer())
    if (disease == "CDG") {
        cohort.samples <- read.table(file.path(workdir, disease, "samples.txt"), sep = "\t", header = TRUE) %>%
            filter(Provider != "Lan Lin" & ID != sample.map["PMD-C01"]) %>% pull(ID)
    } else {
        cohort.samples <- read.table(file.path(workdir, disease, "samples.txt"), sep = "\t", header = TRUE) %>%
            filter(Provider == "Rebecca Ganetzky" | (Provider == "Marni Falk" & Status != "Diagnosed" & 
            !grepl("-MF", ID))) %>% pull(ID)
    }
    tissue.assign <- setNames(rep("Fibroblasts (cohort)", length(cohort.samples)), cohort.samples)
    for (sample.id in cohort.samples) {
        outDF <- full_join(outDF, read.table(file.path(workdir, disease, sample.id, "RNA/stripe/quality_control", 
            paste(sample.id, "_TEQUILA.junctions.tsv", sep = "")), sep = "\t", header = TRUE) %>% 
            filter(paste(Chrom, Start, End, sep = "_") %in% junctions) %>% setNames(c("Chrom", "Start", "End", sample.id)), 
            by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)
    }

    # Remove intermediate file
    system(paste("rm", file.path(workdir, "tmp/input.gtf")))

    # Calculate usage frequencies for annotated junctions from target genes for GTEx fibroblasts, skin, adipose, and whole blood
    for (tissue in c("Fibroblast", "Skin", "Adipose", "Blood")) {
        tissueDF <- read.table(file.path(workdir, disease, "references/GTEx_v8", tissue, "target_genes.junction_counts.txt"), sep = "\t",
            header = TRUE, check.names = FALSE) %>% distinct() %>% filter(Name %in% junctions) %>% separate(Name, c("Chrom", "Start", "End"), 
            sep = "_") %>% mutate_at(c("Start", "End"), as.integer)
        tissue.assign <- c(tissue.assign, setNames(rep(ifelse(tissue == "Fibroblast", "Fibroblasts (GTEx v8)",
            ifelse(tissue == "Blood", "Whole blood (GTEx v8)", paste(tissue, "(GTEx v8)"))), ncol(tissueDF) - 3), 
            colnames(tissueDF)[-(1:3)]))
        outDF <- full_join(outDF, tissueDF, by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)
    }

    outDF <- ComputePSI(outDF)

    # Perform PCA clustering of samples based on usage frequencies of annotated splice junctions represented in outDF
    # Only keep junctions with no more than one-third of the samples have missing usage frequencies
    outDF <- outDF %>% unite("Name", Chrom:End) %>% tibble::column_to_rownames(var = "Name") %>% filter(rowSums(is.na(.))/ncol(.) <= 1/3)
    pcaResults <- pca(t(outDF), method = "nipals", nPcs = 2, scale = "none", center = TRUE)
    scoreDF <- as.data.frame(scores(pcaResults)) %>% tibble::rownames_to_column(var = "Sample_ID") %>%
        mutate(Tissue = tissue.assign[Sample_ID], Panel = disease, Tissue = factor(Tissue, levels = 
        c("Whole blood (GTEx v8)", "Adipose (GTEx v8)", "Skin (GTEx v8)", "Fibroblasts (GTEx v8)", 
        "Fibroblasts (cohort)")))
    plotDF <- bind_rows(plotDF, scoreDF %>% select(Panel, Tissue, PC1, PC2))
}
plotDF$Panel <- factor(plotDF$Panel, levels = c("CDG", "PMD"))
levels(plotDF$Panel) <- c("CDG-466 gene panel", "PMD-359 gene panel")
plotDF$Tissue <- factor(plotDF$Tissue, levels = c("Whole blood (GTEx v8)", "Adipose (GTEx v8)", 
    "Skin (GTEx v8)", "Fibroblasts (GTEx v8)", "Fibroblasts (cohort)"))

p <- ggplot(plotDF[order(plotDF$Tissue), ], aes(x = PC1, y = PC2, color = Tissue)) + facet_wrap(. ~ Panel, nrow = 1) + 
    geom_point(size = 1, stroke = NA) + theme_classic() + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    axis.line = element_line(color = "black", linewidth = 0.25), strip.text = element_text(color = "black", size = 6),
    legend.text = element_text(color = "black", size = 6), legend.title = element_blank(), legend.key.size = unit(0.3, 'cm')) + 
    scale_color_manual(values = c("#D6604D", "#E0C27C", "#689E45", "#5A8BA8", "#94C1B6"))

# Save p to outfile
ggsave(outfile, plot = p, width = 4, height = 1.75)