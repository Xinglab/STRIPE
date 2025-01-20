#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.12.27
# Supplementary Figure 2b

# On-target rates for CDG-466 and PMD-359 gene panels across TEQUILA-seq samples generated on fibroblasts from 88 cohort individuals.

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
outfile <- file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S2/Figure_S2b.pdf")

# Construct dataframe summarizing on-target rates for CDG-466 and PMD-359 gene panels
# across TEQUILA-seq samples generated on fibroblasts from 88 cohort individuals

outDF = tibble(Label = character(), Sample_ID = character(), Rate = numeric())
for (disease in c("CDG", "PMD")) {
    if (disease == "CDG") {
        input.samples <- read.table(file.path(workdir, disease, "samples.txt"), sep = "\t", header = TRUE) %>%
            filter(Provider != "Lan Lin" & ID != "E1877") %>% pull(ID)
        label <- paste("CDG-466", paste("(n = ", length(input.samples), ")", sep = ""), sep = "\n")
        target.genes <- read.table(file.path(workdir, disease, "references/target_genes.bed"), sep = "\t", header = FALSE) %>% pull(V4) 
    } else {
        input.samples <- read.table(file.path(workdir, disease, "samples.txt"), sep = "\t", header = TRUE) %>%
            filter(Provider == "Rebecca Ganetzky" | (Provider == "Marni Falk" & Status != "Diagnosed" & 
            !grepl("-MF", ID))) %>% pull(ID)
        label <- paste("PMD-359", paste("(n = ", length(input.samples), ")", sep = ""), sep = "\n")
        target.genes <- read.table(file.path(workdir, disease, "references/target_genes.bed"), sep = "\t", header = FALSE) %>% pull(V4) 
    }
    for (sample.id in input.samples) {
        outDF <- bind_rows(outDF, tibble(Label = label, Sample_ID = sample.id, 
            Rate = read.table(file.path(workdir, disease, sample.id, "RNA/stripe/quality_control",
            paste(sample.id, "TEQUILA.mapping_stats.txt", sep = "_")), sep = "\t", header = FALSE) %>% 
            pull(V1) %>% tail(1) %>% gsub("%", "", .) %>% as.numeric()))
    }    
}

# Plot on-target rates for each gene panel
p <- ggplot(outDF, aes(x = Label, y = Rate, color = Label)) + geom_boxplot() + theme_classic() + ylim(0, 100) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.text = element_text(color = "black", size = 6), axis.title.y = element_text(color = "black", size = 7), 
    axis.title.x = element_blank(), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    axis.line = element_line(color = "black", linewidth = 0.25), legend.position = "none") + 
    ylab("On-target rate (%)") + scale_color_manual(values = c("#3188BD", "#AADDA3"))

# Save p to outfile
ggsave(outfile, plot = p, width = 1.75, height = 1.75)
