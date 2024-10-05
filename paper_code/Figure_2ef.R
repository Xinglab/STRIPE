#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.04
# Figure 2e,f

# (e) Full-length structures and (f) isoform-level proportions of EARS2 transcripts detected from TEQUILA-seq data of 
# individuals in our study cohort. Shaded and unshaded regions within transcript structures represent putative coding 
# sequences and untranslated regions respectively. The arrow in (f) points to individual BS2-1 (haplotype 1). 

# =====================================================================================================================
#                                                      LIBRARIES 
# =====================================================================================================================

# Load required libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(cowplot))

# =====================================================================================================================
#                                                   HELPER FUNCTIONS
# =====================================================================================================================

PullFeature <- function(infoString, featureName) {
    # Function designed to pull out the value for featureName in infoString
    present <- unlist(lapply(strsplit(infoString, "; "), function(x) lapply(strsplit(x, " "), "[[", 1) == featureName))
    return(ifelse(sum(present) == 1, unlist(lapply(strsplit(unlist(strsplit(infoString, "; "))[present], " "), "[[", 2)), NA))
}

RescaleFeature <- function(x) {
    return(ifelse(x < 50 & x > 0, 50, ifelse(x > 300, 300 + round(sqrt(x)), x)))
}

GetIntrons <- function(x, transcripts) {
    intronDF <- tibble(Transcript_ID = character(), V3 = character(), V4 = numeric(), V5 = numeric(), Transcript_Number = integer())
    for (tx.id in transcripts) {
        intronDF <- bind_rows(intronDF, tibble(Transcript_ID = head(filter(x, Transcript_ID == tx.id)$Transcript_ID, -1),
            V4 = head(filter(x, Transcript_ID == tx.id)$V5, -1), V5 = tail(filter(x, Transcript_ID == tx.id)$V4, -1),
            Transcript_Number = head(filter(x, Transcript_ID == tx.id)$Transcript_Number, -1)) %>%
            mutate(V3 = "intron") %>% select(Transcript_ID, V3, V4, V5, Transcript_Number) %>% filter(V5 > V4))
    }
    return(intronDF)
}

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir <- "/mnt/isilon/lin_lab_share/STRIPE"
stringtie <- "/scr1/users/wangr5/tools/stringtie-2.2.3/stringtie"
samtools <- "/scr1/users/wangr5/tools/samtools-1.21/samtools"
target.gene <- read.table(file.path(workdir, "PMD/references/target_genes.bed"), sep = "\t", header = FALSE) %>% filter(V5 == "EARS2")
gencode.gtf <- "/scr1/users/wangr5/references/gencode.v45.annotation.gtf"
outfile <- file.path(workdir, "manuscript/Main_Figures/Figure_2/Figure_2ef.pdf")

# =====================================================================================================================
#                                                      PANEL E
# =====================================================================================================================

# Create a temporary directory in the folder for Figure 2
dir.create(file.path(dirname(outfile), "tmp"), showWarnings = FALSE)

# Extract transcript annotations for target.gene from gencode.gtf
system(paste("grep \"", target.gene$V4, "\" ", gencode.gtf, " > ", file.path(dirname(outfile), "tmp/gene.gtf"), sep = ""))

# Run stringtie on haplotype-specific BAM file for BS2-1
system(paste(stringtie, file.path(workdir, "PMD/BS2-1/RNA/stripe/target_genes/EARS2/hap1_reads.bam"), "-G", 
    file.path(dirname(outfile), "tmp/gene.gtf"), "-o", file.path(dirname(outfile), "tmp/output.gtf"), "-L -s 5 -c 5 -u -M 0"))
system(paste(stringtie, "--merge -G", file.path(dirname(outfile), "tmp/gene.gtf"), "-o", file.path(dirname(outfile), "tmp/merged.gtf"), 
    "-i", file.path(dirname(outfile), "tmp/output.gtf")))
system(paste("python /scr1/users/wangr5/tools/Annotate_ORF.py -i", file.path(dirname(outfile), "tmp/merged.gtf"),
    "-a", file.path(dirname(outfile), "tmp/gene.gtf"), "-f /scr1/users/wangr5/references/GRCh38.primary_assembly.genome.fa",
    "-o", file.path(dirname(outfile), "tmp/merged.updated.gtf")))

# Estimate abundances of transcripts in merged.updated.gtf from haplotype-specific BAM file for BS2-1
dir.create(file.path(dirname(outfile), "tmp/BS2-1"), showWarnings = FALSE)
system(paste(stringtie, file.path(workdir, "PMD/BS2-1/RNA/stripe/target_genes/EARS2/hap1_reads.bam"), "-G", 
    file.path(dirname(outfile), "tmp/merged.updated.gtf"), "-o", file.path(dirname(outfile), "tmp/BS2-1/output.gtf"), "-L -s 5 -c 5 -u -M 0"))

# Estimate abundances of transcripts in merged.updated.gtf from other cohort samples
cohort.samples <- read.table(file.path(workdir, "PMD/samples.txt"), sep = "\t", header = TRUE) %>% filter(Provider == "Rebecca Ganetzky" &
    ID != "BS2-1") %>% pull(ID)
for (sample.id in cohort.samples) {
    dir.create(file.path(dirname(outfile), "tmp", sample.id), showWarnings = FALSE)
    system(paste(samtools, "view -hb -F 256 -q 1", file.path(workdir, "PMD", sample.id, "RNA", paste(sample.id, "TEQUILA.bam", sep = "_")), 
        paste(target.gene$V1, ":", target.gene$V2+1, "-", target.gene$V3, sep = ""), " > ", file.path(dirname(outfile), "tmp", sample.id, "input.bam")))
    system(paste(samtools, "index", file.path(dirname(outfile), "tmp", sample.id, "input.bam")))
    system(paste(stringtie, file.path(dirname(outfile), "tmp", sample.id, "input.bam"), "-G", file.path(dirname(outfile), "tmp/merged.updated.gtf"), 
        "-o", file.path(dirname(outfile), "tmp", sample.id, "output.gtf"), "-L -u -M 0"))
}

# Assemble an TPM matrix for transcripts in merged.updated.gtf
outDF <- read.table(file.path(dirname(outfile), "tmp/merged.updated.gtf"), sep = "\t", header = FALSE, comment = "#") %>%
    filter(V3 == "transcript") %>% mutate(Transcript_ID = unlist(lapply(V9, function(x) gsub(";", "", PullFeature(x, "transcript_id"))))) %>%
    select(Transcript_ID) %>% tibble
outDF <- left_join(outDF, read.table(file.path(dirname(outfile), "tmp/BS2-1/output.gtf"), sep = "\t", header = FALSE, comment = "#") %>%
    filter(V3 == "transcript") %>% mutate(Transcript_ID = unlist(lapply(V9, function(x) PullFeature(x, "reference_id"))),
    TPM = as.numeric(unlist(lapply(V9, function(x) gsub(";", "", PullFeature(x, "TPM")))))) %>% select(Transcript_ID, TPM) %>% drop_na %>%
    setNames(c("Transcript_ID", "BS2-1")), by = join_by(Transcript_ID)) %>% replace(is.na(.), 0)

for (sample.id in cohort.samples) {
    outDF <- left_join(outDF, read.table(file.path(dirname(outfile), "tmp", sample.id, "output.gtf"), sep = "\t", header = FALSE, comment = "#") %>%
        filter(V3 == "transcript") %>% mutate(Transcript_ID = unlist(lapply(V9, function(x) PullFeature(x, "reference_id"))),
        TPM = as.numeric(unlist(lapply(V9, function(x) gsub(";", "", PullFeature(x, "TPM")))))) %>% select(Transcript_ID, TPM) %>% drop_na %>%
        setNames(c("Transcript_ID", sample.id)), by = join_by(Transcript_ID)) %>% replace(is.na(.), 0)
}

# Convert TPM matrix into a proportion matrix and identify transcripts with an isoform-level proportion of at least 5% in at least one sample
propMatrix <- bind_cols(outDF[,1], sweep(outDF[,-1], 2, colSums(outDF[,-1]), `/`))
keepTranscripts <- sort(propMatrix$Transcript_ID[rowSums(propMatrix[,-1] >= 0.05) >= 1])
propMatrix$Transcript_ID[!(propMatrix$Transcript_ID %in% keepTranscripts)] <- "Other"
propMatrix <- propMatrix %>% group_by(Transcript_ID) %>% summarise(across(everything(), sum)) %>% ungroup

# Parse individual genomic features from merged.updated.gtf and rescale them
gtfDF <- read.table(file.path(dirname(outfile), "tmp/merged.updated.gtf"), sep = "\t", header = FALSE, comment = "#") %>%
    mutate(Transcript_ID = unlist(lapply(V9, function(x) gsub(";", "", PullFeature(x, "transcript_id"))))) %>%
    filter(Transcript_ID %in% keepTranscripts) %>% select(Transcript_ID, V3, V4, V5)
gtfDF <- mutate(gtfDF, V3 = case_when(!(Transcript_ID %in% unique(filter(gtfDF, V3 == "CDS") %>% 
    pull(Transcript_ID))) & V3 == "exon" ~ "UTR", TRUE ~ V3)) %>% filter(V3 %in% c("UTR", "CDS")) %>%
    mutate(V4 = V4 - 1)
oldCoord <- sort(unique(c(pull(gtfDF, V4), pull(gtfDF, V5))))
newCoord <- cumsum(c(0, RescaleFeature(tail(oldCoord, -1) - head(oldCoord, -1))))
gtfDF <- mutate(gtfDF, V4 = recode(V4, !!!setNames(newCoord, oldCoord))/max(newCoord), V5 = recode(V5, !!!setNames(newCoord, oldCoord))/max(newCoord), 
    Transcript_Number = as.numeric(factor(Transcript_ID, levels = rev(keepTranscripts))))
intronDF <- GetIntrons(gtfDF, keepTranscripts)
utrDF <- filter(gtfDF, V3 == "UTR")
cdsDF <- filter(gtfDF, V3 == "CDS")
labelDF <- gtfDF %>% select(Transcript_ID, Transcript_Number) %>% distinct %>% mutate(Transcript_ID = gsub("MSTRG", "NovelTx", Transcript_ID))

palette <- c("#598CA8", "#94C1B6", "#689E45", "#D6604D", "#DFC17D", "#D7869D", "#CAB2D2")
p1 <- ggplot() + geom_rect(data = utrDF, fill = "white", xmin = 1 - utrDF$V4, xmax = 1 - utrDF$V5, ymin = utrDF$Transcript_Number - 0.25,
    ymax = utrDF$Transcript_Number + 0.25, color = "black", linewidth = 0.5) + geom_rect(data = cdsDF %>% mutate(Transcript_Number = factor(Transcript_Number)), 
    aes(fill = Transcript_Number), xmin = 1 - cdsDF$V4, xmax = 1 - cdsDF$V5, ymin = cdsDF$Transcript_Number - 0.25, ymax = cdsDF$Transcript_Number + 0.25, 
    color = "black", linewidth = 0.5) + geom_segment(data = intronDF, x = 1 - intronDF$V4, xend = 1 - intronDF$V5, y = intronDF$Transcript_Number, 
    yend = intronDF$Transcript_Number, linewidth = 0.5, color = "black") + theme_classic() + geom_text(data = labelDF, aes(x = -0.05, y = Transcript_Number, 
    label = Transcript_ID), size = 6*5/14, color = "black", hjust = 1) + coord_cartesian(xlim = c(-0.25, 1), ylim = c(0.5, length(keepTranscripts) + 0.5), 
    clip = "off") + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + scale_fill_manual(values = rev(palette))

# =====================================================================================================================
#                                                      PANEL F
# =====================================================================================================================

propDF <- gather(propMatrix, "Sample_ID", "Proportion", -Transcript_ID)
sampleOrder <- propDF %>% filter(Transcript_ID == "ENST00000449606.7") %>% arrange(Proportion) %>% pull(Sample_ID)
propDF$Sample_ID <- factor(propDF$Sample_ID, levels = sampleOrder)

p2 <- ggplot(propDF %>% mutate(Transcript_ID = gsub("MSTRG", "NovelTx", Transcript_ID)), aes(x = Sample_ID, y = Proportion, fill = Transcript_ID)) + 
    geom_bar(stat = "identity", position = "stack", color = NA) + theme_classic() + ylab("Isoform proportion (EARS2)") + 
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.y = element_line(color = "black", 
    linewidth = 0.25), axis.text.y = element_text(color = "black", size = 6), axis.title.y = element_text(color = "black", size = 7), legend.text = 
    element_text(color = "black", size = 6), legend.title = element_blank(), legend.key.size = unit(0.3, "cm"),
    plot.margin = margin(t = 20, l = 5, unit = "pt"), legend.position = "bottom") + scale_fill_manual(values = c(palette, "#BFBFBF")) +
    guides(fill = guide_legend(ncol = 2))

# Remove intermediate files
system(paste("rm -rf", file.path(dirname(outfile), "tmp")))

# Assemble p1 and p2 onto the same plotting grid
p <- plot_grid(p1, p2, nrow = 1, rel_widths = c(1.75, 1), labels = c("e", "f"), label_size = 8)
ggsave(outfile, plot = p, width = 6.5, height = 2.5)
