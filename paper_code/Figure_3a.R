#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.03
# Figure 3a

# (a) Number of unique RNA processing defects observed for splice site region variants that were identified from clinical 
# testing of previously diagnosed and undiagnosed patients. Variants located in the -3 to +6 region of the donor splice 
# site or the -20 to +3 region of the acceptor splice site were classified as "splice site region" variants. 

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
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/Figure_3a.pdf")
outDF <- tibble(ID = character(), Gene = character(), Variant = character(), Variant_Class = character(), Position = character(),
    Index = integer(), Exon_Skipping = integer(), Cryptic_SS = integer(), Intron_Retention = integer(), Intron_Polyadenylation = integer())

# =====================================================================================================================
#                                                      PMD-P06
# =====================================================================================================================

# MRPS34:
#   * chr16:1772656 C>T (NM_023936.2:c.322-10G>A)
#   * Splice acceptor region (Position: -10)
#   * Intron 1 retention, cryptic acceptor splice site (intron 1)

outDF <- bind_rows(outDF, tibble(ID = "PMD-P06", Gene = "MRPS34", Variant = "chr16:1772656 C>T", Variant_Class = "A", 
    Position = "-10", Index = 1, Exon_Skipping = 0, Cryptic_SS = 1, Intron_Retention = 1, Intron_Polyadenylation = 0))

# =====================================================================================================================
#                                                      PMD-P11
# =====================================================================================================================

# SUCLG1:
#   * chr2:84449761 T>C (NM_003849.4:c.98-9A>G)
#   * Splice acceptor region (Position: -9)
#   * Cryptic acceptor splice site (intron 1)

outDF <- bind_rows(outDF, tibble(ID = "PMD-P11", Gene = "SUCLG1", Variant = "chr2:84449761 T>C", Variant_Class = "A", 
    Position = "-9", Index = 1, Exon_Skipping = 0, Cryptic_SS = 1, Intron_Retention = 0, Intron_Polyadenylation = 0))

# =====================================================================================================================
#                                                      PMD-P04
# =====================================================================================================================

# ELAC2:
#   * chr17:13016933 CT>A (NM_018127.7:c.297-2_297-1delinsT)
#   * Splice acceptor (Position: -2)
#   * Cryptic acceptor splice site (exon 3)

outDF <- bind_rows(outDF, tibble(ID = "PMD-P04", Gene = "ELAC2", Variant = "chr17:13016933 CT>A", Variant_Class = "A", 
    Position = "-2", Index = 1, Exon_Skipping = 0, Cryptic_SS = 1, Intron_Retention = 0, Intron_Polyadenylation = 0))

# =====================================================================================================================
#                                                      PMD-P12
# =====================================================================================================================

# TRIT1:
#   * chr1:39850120 T>C (NM_017646.6:c.702A>G, p.Ala234=)
#   * Splice donor region (Position: -2)
#   * Exons 4 and 5 skipping, exons 2-5 skipping

outDF <- bind_rows(outDF, tibble(ID = "PMD-P12", Gene = "TRIT1", Variant = "chr1:39850120 T>C", Variant_Class = "D", 
    Position = "-2", Index = 1, Exon_Skipping = 2, Cryptic_SS = 0, Intron_Retention = 0, Intron_Polyadenylation = 0))

# =====================================================================================================================
#                                                      CDG-P05
# =====================================================================================================================

# PIGN:
#   * chr18:62143306 C>T (NM_176787.5:c.963G>A, p.Gln321=)
#   * Splice donor region (Position: -1)
#   * Exon 11 skipping, cryptic donor splice site (intron 11), intron 11 polyadenylation

outDF <- bind_rows(outDF, tibble(ID = "CDG-P05", Gene = "PIGN", Variant = "chr18:62143306 C>T", Variant_Class = "D", 
    Position = "-1", Index = 1, Exon_Skipping = 1, Cryptic_SS = 1, Intron_Retention = 0, Intron_Polyadenylation = 1))

# =====================================================================================================================
#                                                      CDG-P02
# =====================================================================================================================

# PIGN:
#   * chr18:62106985 C>G (NM_176787.5:c.1674+1G>C)
#   * Splice donor (Position: +1)
#   * Exon 18 skipping

outDF <- bind_rows(outDF, tibble(ID = "CDG-P02", Gene = "PIGN", Variant = "chr18:62106985 C>G", Variant_Class = "D", 
    Position = "+1", Index = 1, Exon_Skipping = 1, Cryptic_SS = 0, Intron_Retention = 0, Intron_Polyadenylation = 0))

# =====================================================================================================================
#                                                      CDG-P04
# =====================================================================================================================

# PIGN:
#   * chr18:62106985 C>G (NM_176787.5:c.1674+1G>C)
#   * Splice donor (Position: +1)
#   * Exon 18 skipping

outDF <- bind_rows(outDF, tibble(ID = "CDG-P04", Gene = "PIGN", Variant = "chr18:62106985 C>G", Variant_Class = "D", 
    Position = "+1", Index = 2, Exon_Skipping = 1, Cryptic_SS = 0, Intron_Retention = 0, Intron_Polyadenylation = 0))

# =====================================================================================================================
#                                                      PMD-P01
# =====================================================================================================================

# ATP5MK:
#   * chr10:103392370 C>G (NM_001206427.2:c.87+1G>C)
#   * Splice donor (Position: +1)
#   * Exon 3 skipping, exons 2 and 3 skipping, cryptic donor splice site (exon 3)

outDF <- bind_rows(outDF, tibble(ID = "PMD-P01", Gene = "ATP5MK", Variant = "chr10:103392370 C>G", Variant_Class = "D", 
    Position = "+1", Index = 3, Exon_Skipping = 2, Cryptic_SS = 1, Intron_Retention = 0, Intron_Polyadenylation = 0))

# =====================================================================================================================
#                                                      PMD-P07
# =====================================================================================================================

# NUBPL:
#   * chr14:31826715 G>A (NM_025152.3:c.693+1G>A)
#   * Splice donor (Position: +1)
#   * Intron 8 polyadenylation

outDF <- bind_rows(outDF, tibble(ID = "PMD-P07", Gene = "NUBPL", Variant = "chr14:31826715 G>A", Variant_Class = "D", 
    Position = "+1", Index = 4, Exon_Skipping = 0, Cryptic_SS = 0, Intron_Retention = 0, Intron_Polyadenylation = 1))

# =====================================================================================================================
#                                                      CDG-P09
# =====================================================================================================================

# PIGQ:
#   * chr16:576255 G>A (NM_004204.5:c.942+1G>A)
#   * Splice donor (Position: +1)
#   * Exon 4 skipping, intron 4 polyadenylation

outDF <- bind_rows(outDF, tibble(ID = "CDG-P09", Gene = "PIGQ", Variant = "chr16:576255 G>A", Variant_Class = "D", 
    Position = "+1", Index = 5, Exon_Skipping = 1, Cryptic_SS = 0, Intron_Retention = 0, Intron_Polyadenylation = 1))

# =====================================================================================================================
#                                                      PMD-P02
# =====================================================================================================================

# ATP5PO:
#   * chr21:33914447 T>C (NM_001697.3:c.87+3A>G)
#   * Splice donor region (Position: +3)
#   * Exon 2 skipping

outDF <- bind_rows(outDF, tibble(ID = "PMD-P02", Gene = "ATP5PO", Variant = "chr21:33914447 T>C", Variant_Class = "D", 
    Position = "+3", Index = 1, Exon_Skipping = 1, Cryptic_SS = 0, Intron_Retention = 0, Intron_Polyadenylation = 0))

# =====================================================================================================================
#                                                      PMD-P05
# =====================================================================================================================

# FBXL4:
#   * chr6:98880546 CACTT>C (NM_001278716.2:c.1389+3_1389+6del)
#   * Splice donor region (Position: +3)
#   * Exon 8 skipping

outDF <- bind_rows(outDF, tibble(ID = "PMD-P05", Gene = "FBXL4", Variant = "chr6:98880546 CACTT>C", Variant_Class = "D", 
    Position = "+3", Index = 2, Exon_Skipping = 1, Cryptic_SS = 0, Intron_Retention = 0, Intron_Polyadenylation = 0))

# =====================================================================================================================
#                                                      PMD-P08
# =====================================================================================================================

# OPA1:
#   * chr3:193643448 AG>A (NM_130837.3:c.1377+5del)
#   * Splice donor region (Position: +5)
#   * Exon 14 skipping

outDF <- bind_rows(outDF, tibble(ID = "PMD-P08", Gene = "OPA1", Variant = "chr3:193643448 AG>A", Variant_Class = "D", 
    Position = "+5", Index = 1, Exon_Skipping = 1, Cryptic_SS = 0, Intron_Retention = 0, Intron_Polyadenylation = 0))

# =====================================================================================================================
#                                                      CDG-P02
# =====================================================================================================================

# PIGN:
#   * chr18:62113129 C>T (NM_176787.5:c.1434+5G>A)
#   * Splice donor (Position: +5)
#   * Exon 16 skipping, intron 16 polyadenylation

outDF <- bind_rows(outDF, tibble(ID = "CDG-P02", Gene = "PIGN", Variant = "chr18:62113129 C>T", Variant_Class = "D", 
    Position = "+5", Index = 2, Exon_Skipping = 1, Cryptic_SS = 0, Intron_Retention = 0, Intron_Polyadenylation = 1))

# =====================================================================================================================
#                                                      CDG-P06
# =====================================================================================================================

# PIGN:
#   * chr18:62113129 C>T (NM_176787.5:c.1434+5G>A)
#   * Splice donor (Position: +5)
#   * Exon 16 skipping, intron 16 polyadenylation

outDF <- bind_rows(outDF, tibble(ID = "CDG-P06", Gene = "PIGN", Variant = "chr18:62113129 C>T", Variant_Class = "D", 
    Position = "+5", Index = 3, Exon_Skipping = 1, Cryptic_SS = 0, Intron_Retention = 0, Intron_Polyadenylation = 1))

outDF$Position <- paste(outDF$Variant_Class, outDF$Position, sep = "")
outDF$Position <- factor(outDF$Position, levels = unique(outDF$Position))
outDF <- outDF %>% gather("Event", "Count", -c(ID, Gene, Variant, Variant_Class, Position, Index)) %>%
    mutate(Event = factor(Event, levels = c("Intron_Polyadenylation", "Intron_Retention", "Cryptic_SS", "Exon_Skipping")))
levels(outDF$Event) <- c("Intron polyadenylation", "Intron retention", "Cryptic splice site activation", "Exon skipping")

p <- ggplot(outDF %>% mutate(Count = factor(Count), Index = factor(Index)), aes(x = Index, y = Event, fill = Count)) + 
    facet_grid(cols = vars(Position), scales = "free_x", space = "free_x") + geom_tile(color = "black", linewidth = 0.25) + 
    geom_text(aes(label = Count), size = 6*5/14) + theme_classic() + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.y = element_text(color = "black", size = 6), 
    axis.text.x = element_text(color = "black", size = 6, angle = -30, hjust = 0), axis.title = element_blank(), axis.line = element_blank(), 
    legend.position = "right", legend.key.size = unit(0.4, "cm"), legend.text = element_text(color = "black", size = 6), 
    legend.title = element_text(color = "black", size = 6), axis.ticks = element_line(color = "black", linewidth = 0.25)) + 
    scale_fill_manual(values = c("#FFFFFF", "#C5E7C2", "#66C2A5")) + guides(fill = guide_legend(title = "Number of unique events"))

ggsave(outfile, plot = p, width = 6.5, height = 1.5)