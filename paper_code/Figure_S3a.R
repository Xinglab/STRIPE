#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.10.07
# Supplementary Figure 3a

# Comparison of target and non-target gene abundances between untargeted long-read RNA-seq and TEQUILA-seq data for fibroblast 
# cell lines CDG-P12 (CDG-466 gene panel, left) and PMD-C01 (PMD-359 gene panel, right). Only target genes with TPM (transcripts 
# per million) > 1 in untargeted data are displayed. 

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

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir <- "/mnt/isilon/lin_lab_share/STRIPE"
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Supplementary_Figures/Figure_S3/Figure_S3a.pdf")

# Read in mapping file describing new IDs for cohort samples
sample.map <- read.table(file.path(workdir, "manuscript/Revisions/20250228/sample.map"), sep = "\t", header = FALSE) %>% 
    select(V2, V1) %>% tibble::deframe()

# Read in target and non-target gene abundances for untargeted long-read RNA-seq and TEQUILA-seq data on fibroblast cell lines 
# CDG-P12 (CDG-466 gene panel) and PMD-C01 (PMD-359 gene panel)

outDF <- tibble(Gene_ID = character(), Category = character(), TPM = numeric(), Library = character(), Panel = character())
for (disease in c("CDG", "PMD")) {
    sample.id <- ifelse(disease == "CDG", sample.map["CDG-P12"], sample.map["PMD-C01"])
    target.genes <- read.table(file.path(workdir, disease, "references/target_genes.txt"))$V1
    for (lib in c("ONT", "TEQUILA")) {
        input.data <- read.table(file.path(workdir, disease, sample.id, "RNA", ifelse(lib == "ONT", "stripe_qc", "stripe"),
            "quality_control/stringtie", paste(sample.id, "_", lib, ".gene_abundance.tsv", sep = "")), sep = "\t", header = TRUE) %>%
            filter(grepl("ENSG", Gene.ID)) %>% separate(Gene.ID, c("Gene_ID", NA), sep = "\\.") %>%
            mutate(Category = ifelse(Gene_ID %in% target.genes, "Target", "Non-target")) %>% select(Gene_ID, Category, TPM) %>% 
            mutate(Library = lib, Panel = disease) 
        outDF <- bind_rows(outDF, input.data)
    }
}
outDF$Panel <- factor(outDF$Panel, levels = c("CDG", "PMD"))
levels(outDF$Panel) <- c("CDG-466 gene panel", "PMD-359 gene panel")
outDF$Category <- factor(outDF$Category, levels = c("Non-target", "Target"))

# Filter for genes with TPM > 1 in untargeted long-read RNA-seq
outDF <- outDF %>% spread(Library, TPM) %>% filter(ONT > 1)
    
# Compute Spearman's correlation in target gene and non-target gene abundances between untargeted long-read RNA-seq and TEQUILA-seq data
corrDF <- outDF %>% group_by(Panel, Category) %>% summarise(Corr = cor(ONT, TEQUILA, method = "spearman")) %>% mutate(Label = paste("Spearman's r (",
    tolower(Category), ") = ", format(round(Corr, digits = 2), nsmall = 2), sep = "")) %>% group_by(Panel) %>%
    summarize(Label = paste(Label, collapse = "\n"), .groups = "drop")

# Plot target and non-target gene abundances between untargeted long-read RNA-seq and TEQUILA-seq data
p <- ggplot(outDF, aes(x = log2(ONT + 1), y = log2(TEQUILA + 1))) + facet_wrap(. ~ Panel, nrow = 1) + theme_classic() + 
    geom_point(aes(color = Category), stroke = NA, size = 1, alpha = 0.5) + geom_segment(data = tibble(Panel = factor(c("CDG-466 gene panel", 
    "PMD-359 gene panel"), levels = c("CDG-466 gene panel", "PMD-359 gene panel")), XStart = c(0, 0), XEnd = c(15, 15), YStart = c(0, 0), 
    YEnd = c(15, 15)), aes(x = XStart, xend = XEnd, y = YStart, yend = YEnd), linetype = "dashed", linewidth = 0.25) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.text = element_text(color = "black", size = 6), axis.title = element_text(color = "black", size = 7), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), strip.text = element_text(color = "black", size = 6),
    axis.line = element_line(color = "black", linewidth = 0.25), plot.title = element_text(color = "black",
    size = 7, hjust = 0.5), legend.text = element_text(color = "black", size =6), legend.title = element_blank(), 
    legend.key.size = unit(0.3, "cm"), legend.position = "bottom") + ggtitle("Normalized gene abundance") + xlab("Untargeted (cDNA-PCR), log2(TPM+1)") + ylab("TEQUILA-seq, log2(TPM+1)") + 
    geom_text(data = corrDF, mapping = aes(x = 4, y = 2, label = Label), size = 5/14*6, hjust = 0) + scale_color_manual(values = c("#D3D3D3", "#D83F50"))

# Save p to outfile
ggsave(outfile, plot = p, width = 4, height = 2.75)
