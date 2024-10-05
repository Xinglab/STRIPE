#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.05
# Figure 3a

# (a) Number of aberrant RNA processing events observed around known pathogenic splice site variants using TEQUILA-seq 
# data of previously diagnosed patients. 

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
outfile <- file.path(workdir, "manuscript/Main_Figures/Figure_3/Figure_3a.pdf")
outDF <- tibble(ID = character(), Panel = character(), Gene = character(), Variant = character(), Variant_Class = character(), 
    Exon_Skipping = integer(), Cryptic_SS = integer(), Intron_Retention = integer(), Intron_Polyadenylation = integer())

# =====================================================================================================================
#                                                     AnJa (CDG)
# =====================================================================================================================

# PIGN:
#   * chr18:62106985 C>G (splice donor): exon 18 skipping
#   * chr18:62113129 C>T (splice donor): exon 16 skipping, activation of a cryptic polyadenylation site within intron 16

outDF <- bind_rows(outDF, tibble(ID = rep("AnJa", 2), Panel = rep("CDG-466", 2), Gene = rep("PIGN", 2),
    Variant = c("chr18:62106985 C>G", "chr18:62113129 C>T"), Variant_Class = rep("Splice donor", 2),
    Exon_Skipping = c(1, 1), Cryptic_SS = c(0, 0), Intron_Retention = c(0, 0), Intron_Polyadenylation = c(0, 1)))

# =====================================================================================================================
#                                                  CDG-132-1 (CDG)
# =====================================================================================================================

# PIGQ:
#   * chr16:576255 G>A (splice donor): exon 4 skipping, intron 3 retention, activation of a cryptic polyadenylation site 
#     within intron 4

outDF <- bind_rows(outDF, tibble(ID = rep("CDG-132-1", 1), Panel = rep("CDG-466", 1), Gene = rep("PIGQ", 1),
    Variant = c("chr16:576255 G>A"), Variant_Class = rep("Splice donor", 1), Exon_Skipping = c(1), Cryptic_SS = c(0), 
    Intron_Retention = c(0), Intron_Polyadenylation = c(1)))

# =====================================================================================================================
#                                                  CDG-137-1 (CDG)
# =====================================================================================================================

# PIGN:
#   * chr18:62106985 C>G (splice donor): exon 18 skipping

outDF <- bind_rows(outDF, tibble(ID = rep("CDG-137-1", 1), Panel = rep("CDG-466", 1), Gene = rep("PIGN", 1),
    Variant = c("chr18:62106985 C>G"), Variant_Class = rep("Splice donor", 1), Exon_Skipping = c(1), Cryptic_SS = c(0), 
    Intron_Retention = c(0), Intron_Polyadenylation = c(0)))

# =====================================================================================================================
#                                                  CDG-147-1 (CDG)
# =====================================================================================================================

# PIGN:
#   * chr18:62143306 C>T (splice donor): exon 11 skipping, activation of a cryptic donor splice site in intron 11, 
#     activation of a cryptic polyadenylation site within intron 11

outDF <- bind_rows(outDF, tibble(ID = rep("CDG-147-1", 1), Panel = rep("CDG-466", 1), Gene = rep("PIGN", 1),
    Variant = c("chr18:62143306 C>T"), Variant_Class = rep("Splice donor", 1), Exon_Skipping = c(1), Cryptic_SS = c(1), 
    Intron_Retention = c(0), Intron_Polyadenylation = c(1)))

# =====================================================================================================================
#                                                  CDG-161-1 (CDG)
# =====================================================================================================================

# PIGN:
#   * chr18:62113129 C>T (splice donor): exon 16 skipping, activation of a cryptic polyadenylation site within intron 16

outDF <- bind_rows(outDF, tibble(ID = rep("CDG-161-1", 1), Panel = rep("CDG-466", 1), Gene = rep("PIGN", 1),
    Variant = c("chr18:62113129 C>T"), Variant_Class = rep("Splice donor", 1), Exon_Skipping = c(1), Cryptic_SS = c(0), 
    Intron_Retention = c(0), Intron_Polyadenylation = c(1)))

# =====================================================================================================================
#                                                    Q1363 (PMD)
# =====================================================================================================================

# MRPS34:
#   * chr16:1772656 C>T (splice acceptor): intron 1 retention, creation of a cryptic acceptor splice site within intron 1

outDF <- bind_rows(outDF, tibble(ID = rep("Q1363", 1), Panel = rep("PMD-359", 1), Gene = rep("MRPS34", 1),
    Variant = c("chr16:1772656 C>T"), Variant_Class = rep("Splice acceptor", 1), Exon_Skipping = c(0), Cryptic_SS = c(1), 
    Intron_Retention = c(1), Intron_Polyadenylation = c(0)))

# =====================================================================================================================
#                                                    Q1513 (PMD)
# =====================================================================================================================

# FBXL4:
#   * chr6:98880546 CACTT>C (splice donor): exon 8 skipping

outDF <- bind_rows(outDF, tibble(ID = rep("Q1513", 1), Panel = rep("PMD-359", 1), Gene = rep("FBXL4", 1),
    Variant = c("chr6:98880546 CACTT>C"), Variant_Class = rep("Splice donor", 1), Exon_Skipping = c(1), Cryptic_SS = c(0), 
    Intron_Retention = c(0), Intron_Polyadenylation = c(0)))

# =====================================================================================================================
#                                                    Q1663 (PMD)
# =====================================================================================================================

# ATP5MK:
#   * chr10:103392370 C>G (splice donor): exon 3 skipping, skipping of exons 2 and 3, activation of a cryptic donor splice 
#     site within exon 3

outDF <- bind_rows(outDF, tibble(ID = rep("Q1663", 1), Panel = rep("PMD-359", 1), Gene = rep("ATP5MK", 1),
    Variant = c("chr10:103392370 C>G"), Variant_Class = rep("Splice donor", 1), Exon_Skipping = c(2), Cryptic_SS = c(1), 
    Intron_Retention = c(0), Intron_Polyadenylation = c(0)))

# =====================================================================================================================
#                                                    Q1687 (PMD)
# =====================================================================================================================

# NUBPL:
#   * chr14:31826715 G>A (splice donor): activation of a cryptic polyadenylation site within intron 8

outDF <- bind_rows(outDF, tibble(ID = rep("Q1687", 1), Panel = rep("PMD-359", 1), Gene = rep("NUBPL", 1),
    Variant = c("chr14:31826715 G>A"), Variant_Class = rep("Splice donor", 1), Exon_Skipping = c(0), Cryptic_SS = c(0), 
    Intron_Retention = c(0), Intron_Polyadenylation = c(1)))

# =====================================================================================================================
#                                                    Q1695 (PMD)
# =====================================================================================================================

# ELAC2:
#   * chr17:13016933 CT>A (splice acceptor): activation of a cryptic acceptor splice site within exon 3

outDF <- bind_rows(outDF, tibble(ID = rep("Q1695", 1), Panel = rep("PMD-359", 1), Gene = rep("ELAC2", 1),
    Variant = c("chr17:13016933 CT>A"), Variant_Class = rep("Splice acceptor", 1), Exon_Skipping = c(0), Cryptic_SS = c(1), 
    Intron_Retention = c(0), Intron_Polyadenylation = c(0)))

# =====================================================================================================================
#                                                    Q1695 (PMD)
# =====================================================================================================================

# ELAC2:
#   * chr17:13016933 CT>A (splice acceptor): activation of a cryptic acceptor splice site within exon 3

outDF <- bind_rows(outDF, tibble(ID = rep("Q1695", 1), Panel = rep("PMD-359", 1), Gene = rep("ELAC2", 1),
    Variant = c("chr17:13016933 CT>A"), Variant_Class = rep("Splice acceptor", 1), Exon_Skipping = c(0), Cryptic_SS = c(1), 
    Intron_Retention = c(0), Intron_Polyadenylation = c(0)))

# =====================================================================================================================
#                                                    Q1819 (PMD)
# =====================================================================================================================

# TRIT1:
#   * chr1:39850120 T>C (splice donor): skipping of exons 4 and 5, skipping of exons 2 to 5

outDF <- bind_rows(outDF, tibble(ID = rep("Q1819", 1), Panel = rep("PMD-359", 1), Gene = rep("TRIT1", 1),
    Variant = c("chr1:39850120 T>C"), Variant_Class = rep("Splice donor", 1), Exon_Skipping = c(2), Cryptic_SS = c(0), 
    Intron_Retention = c(0), Intron_Polyadenylation = c(0)))

# =====================================================================================================================
#                                                    Q1964 (PMD)
# =====================================================================================================================

# ATP5PO:
#   * chr21:33914447 T>C (splice donor): exon 2 skipping

outDF <- bind_rows(outDF, tibble(ID = rep("Q1964", 1), Panel = rep("PMD-359", 1), Gene = rep("ATP5PO", 1),
    Variant = c("chr21:33914447 T>C"), Variant_Class = rep("Splice donor", 1), Exon_Skipping = c(1), Cryptic_SS = c(0), 
    Intron_Retention = c(0), Intron_Polyadenylation = c(0)))

# =====================================================================================================================
#                                                   Q2032s1 (PMD)
# =====================================================================================================================

# OPA1:
#   * chr3:193643448 AG>A (splice donor): exon 14 skipping

outDF <- bind_rows(outDF, tibble(ID = rep("Q2032s1", 1), Panel = rep("PMD-359", 1), Gene = rep("OPA1", 1),
    Variant = c("chr3:193643448 AG>A"), Variant_Class = rep("Splice donor", 1), Exon_Skipping = c(1), Cryptic_SS = c(0), 
    Intron_Retention = c(0), Intron_Polyadenylation = c(0)))

# =====================================================================================================================
#                                                   Q2319 (PMD)
# =====================================================================================================================

# SUCLG1:
#   * chr2:84449761 T>C (splice acceptor): creation of a cryptic acceptor splice site within intron 1

outDF <- bind_rows(outDF, tibble(ID = rep("Q2319", 1), Panel = rep("PMD-359", 1), Gene = rep("SUCLG1", 1),
    Variant = c("chr2:84449761 T>C"), Variant_Class = rep("Splice acceptor", 1), Exon_Skipping = c(0), Cryptic_SS = c(1), 
    Intron_Retention = c(0), Intron_Polyadenylation = c(0)))

outDF <- outDF %>% unite("Index", ID:Variant) %>% gather("Event", "Count", -c(Index, Variant_Class)) %>% 
    mutate(Variant_Class = factor(Variant_Class, levels = c("Splice donor", "Splice acceptor")),
    Event = factor(Event, levels = c("Intron_Polyadenylation", "Intron_Retention", "Cryptic_SS", "Exon_Skipping")))
levels(outDF$Event) <- c("Intron polyadenylation", "Intron retention", "Cryptic splice site activation", "Exon skipping")
sampleOrder <- outDF %>% mutate(Weight = ifelse(Event == "Intron polyadenylation", 1.5, 1) * Count) %>% select(Index, Weight) %>% 
    group_by(Index) %>% summarise(Total = sum(Weight)) %>% arrange(desc(Total)) %>% pull(Index)

p <- ggplot(outDF %>% mutate(Count = factor(Count), Index = factor(Index, levels = sampleOrder)), aes(x = Index, y = Event, fill = Count)) + 
    geom_tile() + facet_grid(cols = vars(Variant_Class), scales = "free_x", space = "free_x") + geom_text(aes(label = Count), size = 6*5/14) +
    theme_bw() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), 
    axis.ticks.y = element_line(color = "black", linewidth = 0.25), axis.text.y = element_text(color = "black", size = 6), 
    axis.title.y = element_blank(), legend.position = "none", strip.text = element_text(color = "black", size = 6)) +
    scale_fill_manual(values = c("#FFFFFF", "#C5E7C2", "#66C2A5"))

ggsave(outfile, plot = p, width = 4.5, height = 1.25)
