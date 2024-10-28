#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.05
# Figure 3d,e

# (e) Full-length structures and (f) isoform-level proportions of ATP5MK transcripts detected from TEQUILA-seq data of 
# individuals in our study cohort. Shaded and unshaded regions within transcript structures represent putative coding 
# sequences and untranslated regions respectively.

# =====================================================================================================================
#                                                      LIBRARIES 
# =====================================================================================================================

# Load required libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(cowplot))
suppressMessages(library(RColorBrewer))

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
target.gene <- read.table(file.path(workdir, "PMD/references/target_genes.bed"), sep = "\t", header = FALSE) %>% filter(V5 == "ATP5MK")
gencode.gtf <- "/scr1/users/wangr5/references/gencode.v45.annotation.gtf"
outfile <- file.path(workdir, "manuscript/Main_Figures/Figure_3/Figure_3de.pdf")

# =====================================================================================================================
#                                                      PANEL D
# =====================================================================================================================

# Create a temporary directory in the folder for Figure 3
dir.create(file.path(dirname(outfile), "tmp"), showWarnings = FALSE)

# Extract transcript annotations for target.gene from gencode.gtf
system(paste("grep \"", target.gene$V4, "\" ", gencode.gtf, " > ", file.path(dirname(outfile), "tmp/gene.gtf"), sep = ""))

# Run stringtie on TEQUILA-seq BAM files for cohort samples
cohort.samples <- read.table(file.path(workdir, "PMD/samples.txt"), sep = "\t", header = TRUE) %>% filter(Provider == "Rebecca Ganetzky") %>% pull(ID)
for (sample.id in cohort.samples) {
    dir.create(file.path(dirname(outfile), "tmp", sample.id), showWarnings = FALSE)
    system(paste(samtools, "view -hb -F 256 -q 1", file.path(workdir, "PMD", sample.id, "RNA", paste(sample.id, "TEQUILA.bam", sep = "_")), 
        paste(target.gene$V1, ":", target.gene$V2+1, "-", target.gene$V3, sep = ""), " > ", file.path(dirname(outfile), "tmp", sample.id, "input.bam")))
    system(paste(samtools, "index", file.path(dirname(outfile), "tmp", sample.id, "input.bam")))
    # Use -E 10 here since the default settings for stringtie "over-correct" the cryptic splice site observed in Q1663 
    system(paste(stringtie, file.path(dirname(outfile), "tmp", sample.id, "input.bam"), "-G", file.path(dirname(outfile), "tmp/gene.gtf"), "-o", 
        file.path(dirname(outfile), "tmp", sample.id, "output.gtf"), "-L -E 10 -s 5 -c 5 -u -M 0"))
}

# Merge sample-specific GTF files with existing GENCODE annotations (use -f 0 during merging to include the exon 3 skipping isoform)
data.table::fwrite(list(file.path(dirname(outfile), "tmp", cohort.samples, "output.gtf")), file = file.path(dirname(outfile), "tmp/mergelist.txt"))
system(paste(stringtie, "--merge", "-G", file.path(dirname(outfile), "tmp/gene.gtf"), "-o", file.path(dirname(outfile), "tmp/merged.gtf"), "-f 0",
    file.path(dirname(outfile), "tmp/mergelist.txt")))

# Re-run stringtie on TEQUILA-seq BAM files
for (sample.id in cohort.samples) {
    system(paste(stringtie, file.path(dirname(outfile), "tmp", sample.id, "input.bam"), "-G", file.path(dirname(outfile), "tmp/merged.gtf"), "-o", 
        file.path(dirname(outfile), "tmp", sample.id, "output.gtf"), "-s 5 -c 5 -u -M 0 -e -B"))
}

system(paste("python /scr1/users/wangr5/tools/Annotate_ORF.py -i", file.path(dirname(outfile), "tmp/merged.gtf"), "-a", 
    file.path(dirname(outfile), "tmp/gene.gtf"), "-f /scr1/users/wangr5/references/GRCh38.primary_assembly.genome.fa", 
    "-o", file.path(dirname(outfile), "tmp/merged.updated.gtf")))

# Construct FPKM matrix across samples
outDF <- tibble(Transcript_ID = character())
for (sample.id in cohort.samples) {
    outDF <- full_join(outDF, read.table(file.path(dirname(outfile), "tmp", sample.id, "t_data.ctab"), sep = "\t", header = TRUE) %>%
        select(t_name, FPKM) %>% setNames(c("Transcript_ID", sample.id)), by = join_by(Transcript_ID))
}

# Convert TPM matrix into a proportion matrix
propMatrix <- bind_cols(outDF[,1], sweep(outDF[,-1], 2, colSums(outDF[,-1]), `/`))
sampleTx <- select(propMatrix, Transcript_ID, Q1663) %>% arrange(desc(Q1663))
cohortTx <- select(propMatrix, -Q1663) %>% mutate(Total = rowSums(across(where(is.numeric)))) %>% select(Transcript_ID, Total) %>% arrange(desc(Total))
keepTranscripts <- c("ENST00000369815.6", "ENST00000309579.7", "MSTRG.1.4", "MSTRG.1.2", "MSTRG.1.3")
propMatrix$Transcript_ID[!(propMatrix$Transcript_ID %in% keepTranscripts)] <- "Other"
propMatrix <- propMatrix %>% group_by(Transcript_ID) %>% summarise(across(everything(), sum)) %>% ungroup

# Rescale features in gtfDF
gtfDF <- read.table(file.path(dirname(outfile), "tmp/merged.updated.gtf"), sep = "\t", header = FALSE, comment = "#") %>%
    mutate(Transcript_ID = unlist(lapply(V9, function(x) gsub(";", "", PullFeature(x, "transcript_id"))))) %>%
    filter(Transcript_ID %in% keepTranscripts) %>% select(Transcript_ID, V3, V4, V5) %>%
    arrange(Transcript_ID, V4, V5)
gtfDF <- mutate(gtfDF, V3 = case_when(!(Transcript_ID %in% unique(filter(gtfDF, V3 == "CDS") %>% 
    pull(Transcript_ID))) & V3 == "exon" ~ "UTR", TRUE ~ V3)) %>% filter(V3 %in% c("UTR", "CDS")) %>% mutate(V4 = V4 - 1)
oldCoord <- sort(unique(c(pull(gtfDF, V4), pull(gtfDF, V5))))
newCoord <- cumsum(c(0, RescaleFeature(tail(oldCoord, -1) - head(oldCoord, -1))))
gtfDF <- mutate(gtfDF, V4 = recode(V4, !!!setNames(newCoord, oldCoord))/max(newCoord), V5 = recode(V5, !!!setNames(newCoord, oldCoord))/max(newCoord), 
    Transcript_Number = as.numeric(factor(Transcript_ID, levels = rev(keepTranscripts))))
intronDF <- GetIntrons(gtfDF, keepTranscripts)
utrDF <- filter(gtfDF, V3 == "UTR")
cdsDF <- filter(gtfDF, V3 == "CDS")
labelDF <- gtfDF %>% select(Transcript_ID, Transcript_Number) %>% distinct %>% mutate(Transcript_ID = recode(Transcript_ID, !!!setNames( 
    paste("NovelTx", 1:length(unique(gtfDF$Transcript_ID[!grepl("ENST", gtfDF$Transcript_ID)])), sep = "."), 
    unique(gtfDF$Transcript_ID[!grepl("ENST", gtfDF$Transcript_ID)]))))

palette <- setNames(c(brewer.pal(9, "Blues")[c(6,4)], brewer.pal(9, "Reds")[3:5]), seq(length(keepTranscripts), 1))
p1 <- ggplot() + geom_rect(data = utrDF, fill = "white", xmin = 1 - utrDF$V4, xmax = 1 - utrDF$V5, ymin = utrDF$Transcript_Number - 0.25,
    ymax = utrDF$Transcript_Number + 0.25, color = "black", linewidth = 0.5) + geom_rect(data = cdsDF %>% mutate(Transcript_Number = factor(Transcript_Number)), 
    aes(fill = Transcript_Number), xmin = 1 - cdsDF$V4, xmax = 1 - cdsDF$V5, ymin = cdsDF$Transcript_Number - 0.25, ymax = cdsDF$Transcript_Number + 0.25, 
    color = "black", linewidth = 0.5) + geom_segment(data = intronDF, x = 1 - intronDF$V4, xend = 1 - intronDF$V5, y = intronDF$Transcript_Number, 
    yend = intronDF$Transcript_Number, linewidth = 0.5, color = "black") + theme_classic() + geom_text(data = labelDF, aes(x = -0.05, y = Transcript_Number, 
    label = Transcript_ID), size = 6*5/14, color = "black", hjust = 1) + coord_cartesian(xlim = c(-0.25, 1), ylim = c(0.5, length(keepTranscripts) + 0.5), 
    clip = "off") + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + scale_fill_manual(values = rev(palette))

# =====================================================================================================================
#                                                      PANEL E
# =====================================================================================================================

propDF <- gather(propMatrix, "Sample_ID", "Proportion", -Transcript_ID) %>% mutate(Transcript_ID = factor(Transcript_ID, levels = c(keepTranscripts, "Other")))
palette <- setNames(c(brewer.pal(9, "Blues")[c(6,4)], brewer.pal(9, "Reds")[3:5], "#BDBDBD"), c(recode(keepTranscripts, !!!setNames(paste("NovelTx", 
    1:length(unique(gtfDF$Transcript_ID[!grepl("ENST", gtfDF$Transcript_ID)])), sep = "."), unique(gtfDF$Transcript_ID[!grepl("ENST", 
    gtfDF$Transcript_ID)]))), "Other"))
propDF <- propDF %>% mutate(Transcript_ID = factor(recode(Transcript_ID, !!!setNames(paste("NovelTx", 1:length(unique(gtfDF$Transcript_ID[!grepl("ENST", 
    gtfDF$Transcript_ID)])), sep = "."), unique(gtfDF$Transcript_ID[!grepl("ENST", gtfDF$Transcript_ID)]))), levels = names(palette)), 
    Group = factor(ifelse(Sample_ID == "Q1663", "Q1663", "Fibroblasts\n(cohort)"), levels = c("Fibroblasts\n(cohort)", "Q1663")))
summaryDF <- propDF %>% select(-Sample_ID) %>% group_by(Group, Transcript_ID) %>% summarise(Prop_Mean = mean(Proportion))

p2 <- ggplot() + geom_bar(data = summaryDF, aes(x = Group, y = Prop_Mean * 100, color = Transcript_ID, group = Transcript_ID), fill = NA, linewidth = 0.5,
    stat = "identity", position = position_dodge(), width = 0.75) + geom_point(data = propDF, aes(x = Group, y = Proportion * 100, color = Transcript_ID,
    group = Transcript_ID), stroke = NA, size = 1, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.75)) + 
    theme_classic() + ylab("Isoform proportion (%)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.text = element_text(color = "black", size = 6), axis.title.y = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", 
    linewidth = 0.25), plot.title = element_text(color = "black", size = 7, hjust = 0.5, face = "italic"), axis.title.x = element_blank(),
    legend.position = "bottom", legend.key.size = unit(0.3, "cm"), legend.title = element_blank(), legend.text = element_text(color = "black", size = 6)) + 
    ggtitle("ATP5MK") + scale_color_manual(values = palette) + guides(color = guide_legend(ncol = 2)) + ylim(0, 100)
    
# Remove intermediate files
system(paste("rm -rf", file.path(dirname(outfile), "tmp")))

# Assemble p1 and p2 onto the same plotting grid
p <- plot_grid(p1, p2, nrow = 1, rel_widths = c(1.75, 1), labels = c("d", "e"), label_size = 8)
ggsave(outfile, plot = p, width = 6.5, height = 2)
