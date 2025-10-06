#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.10.06
# Figure 1d

# Gene abundances based on untargeted long-read RNA-seq and TEQUILA-seq data for fibroblast cell lines CDG-P12
# (CDG-466 gene panel, left) and PMD-C01 (PMD-359 gene panel, right). Each bar represents one gene, and only the 1000 
# most abundant genes are shown. On-target rates represent the fraction of transcriptional abundance corresponding 
# to target genes. Fold enrichment is calculated by dividing the on-target rate in the capturant by the on-target 
# rate in the unenriched input. 

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
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_1/Figure_1d.pdf")

# Read in mapping file describing new IDs for cohort samples
sample.map <- read.table(file.path(workdir, "manuscript/Revisions/20250228/sample.map"), sep = "\t", header = FALSE) %>% 
    select(V2, V1) %>% tibble::deframe()

# Read in gene abundances for untargeted long-read RNA-seq and TEQUILA-seq data on fibroblast cell lines
# CDG-P12 (CDG-466 gene panel) and PMD-C01 (PMD-359 gene panel)
outDF <- tibble(Gene_ID = character(), Group = character(), TPM = numeric(), Rank = integer(), Library = character(), Panel = character())
for (disease in c("CDG", "PMD")) {
    sample.id <- ifelse(disease == "CDG", sample.map["CDG-P12"], sample.map["PMD-C01"])
    target.genes <- read.table(file.path(workdir, disease, "references/target_genes.txt"))$V1
    for (library in c("ONT", "TEQUILA")) {
        input.data <- read.table(file.path(workdir, disease, sample.id, "RNA", ifelse(library == "ONT", "stripe_qc", "stripe"),
            "quality_control/stringtie", paste(sample.id, "_", library, ".gene_abundance.tsv", sep = "")), sep = "\t", header = TRUE) %>%
            filter(grepl("ENSG", Gene.ID)) %>% separate(Gene.ID, c("Gene_ID", NA), sep = "\\.") %>% 
            mutate(Group = case_when(Gene_ID %in% target.genes ~ "Target", TRUE ~ "Non-target")) %>% select(Gene_ID, Group, TPM) %>% 
            mutate(Rank = rank(-TPM, ties.method = "first"), Library = library, Panel = disease)
        outDF <- bind_rows(outDF, input.data)
    }
}
outDF$Group <- factor(outDF$Group, levels = c("Non-target", "Target"))
outDF$Library <- factor(outDF$Library, levels = c("ONT", "TEQUILA"))
levels(outDF$Library) <- c("Untargeted (cDNA-PCR)", "TEQUILA-seq")
outDF$Panel <- factor(outDF$Panel, levels = c("CDG", "PMD"))
levels(outDF$Panel) <- c("CDG-466 gene panel", "PMD-359 gene panel")

# Compute on-target rates for untargeted long-read RNA-seq and TEQUILA-seq
on.target.rates <- outDF %>% filter(Group == "Target") %>% group_by(Panel, Library) %>% reframe(Rate = sum(TPM)/1e4) %>% 
    mutate(Label = paste("On-target: ", formatC(round(Rate, digits = 1), format = "f", flag = "0", digits = 1), "%", sep = ""))

# Compute the enrichment fold attained by TEQUILA-seq
enrichment.fold <- on.target.rates %>% select(-Label) %>% spread(Library, Rate) %>% mutate(Enrichment = `TEQUILA-seq`/`Untargeted (cDNA-PCR)`) %>% 
    mutate(Label = paste("Enrichment: ", formatC(round(Enrichment, digits = 1), format = "f", flag = "0", digits = 1), "x", sep = "")) %>% 
    mutate(Library = factor("TEQUILA-seq", levels = c("Untargeted (cDNA-PCR)", "TEQUILA-seq")))

# Plot gene rank versus gene abundance for the 1000 most abundant genes in untargeted long-read RNA-seq and TEQUILA-seq
p <- ggplot(outDF %>% filter(Rank <= 1000), aes(x = Rank, y = log2(TPM + 1))) + facet_grid(rows = vars(Panel), cols = vars(Library) ) + 
    theme_classic() + geom_bar(aes(fill = Group), linewidth = 0, stat = "identity") + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    strip.text = element_text(color = "black", size = 6), legend.title = element_blank(), legend.text = element_text(color = "black", size = 6), 
    legend.key.size = unit(0.25, 'cm'), legend.position = "bottom", axis.line = element_line(color = "black", linewidth = 0.25)) + 
    xlab("Genes ranked by abundance") + ylab("Gene abundance, log2(TPM+1)") + geom_text(data = on.target.rates, aes(x = 300, y = 13, label = Label), 
    size = 5/14*6, hjust = 0) + geom_text(data = enrichment.fold, aes(x = 300, y = 11, label = Label), size = 5/14*6, hjust = 0) +
    scale_fill_manual(values = c("#D3D3D3", "#D83F50"))

# Save p to outfile
ggsave(outfile, plot = p, width = 3, height = 3)
