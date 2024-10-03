#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.02
# Supplementary Figure 5a

# Number of known pathogenic small variants in exonic regions that were detected and phased using information from TEQUILA-seq 
# data for previously diagnosed patients. 

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
outfile <- file.path(workdir, "manuscript/Supplementary_Figures/Figure_S5/Figure_S5a.pdf")

# Read in Table_S5.txt and annotate each variant based on whether it was detected or not. Also check whether reads for
# the corresponding gene were phased or not.
inDF <- read.table(file.path(workdir, "manuscript/Supplementary_Tables/Table_S5/Table_S5.txt"), sep = "\t", 
    header = TRUE) %>% mutate(Detected = ifelse(longcallR_gt == "." & clair3_gt == "." & deepvariant_gt == "." &
    phasing_gt == ".", "No", "Yes"))

phased <- c()
for(i in 1:nrow(inDF)){
    phased <- c(phased, ifelse(
        file.exists(file.path(workdir, "manuscript/Supplementary_Tables/Table_S5", inDF[i, "patient_id"], 
        inDF[i, "gene_name"], "phasing_stats.txt")), "Yes", "No"
    ))
}

inDF$Phased <- phased
inDF$Category <- ifelse(inDF$Detected == "No", "Not detected", ifelse(inDF$Phased == "No", "Detected, not phased", "Detected, phased"))
plotDF <- data.frame(table(inDF[c("variant_class", "Category")])) %>% mutate(variant_class = factor(variant_class,
    levels = names(sort(table(inDF$variant_class)))), Category = factor(Category,
    levels = c("Detected, phased", "Detected, not phased", "Not detected")))

p <- ggplot(plotDF, aes(x = Freq, y = variant_class, fill = Category)) + geom_bar(position = "stack", stat = "identity", width = 0.75) +
    theme_classic() + xlab("Number of known pathogenic small variants in exonic regions") + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title.x = element_text(color = "black", size = 7), axis.title.y = element_blank(), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    legend.title = element_blank(), legend.text = element_text(color = "black", size = 6), legend.key.size = unit(0.4, "cm"),
    legend.position = c(0.7, 0.4)) + scale_fill_manual(values = c("#00BA38", "#FFE699", "#EEA3D6"))

ggsave(outfile, plot = p, width = 4, height = 1.75)
