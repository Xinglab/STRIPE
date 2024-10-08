#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.07
# Supplementary Figure 11b,c

# (b) Full-length structures and (c) isoform-level proportions of PIGN transcripts detected from TEQUILA-seq data of 
# individuals in our study cohort. Shaded and unshaded regions within transcript structures represent putative coding 
# sequences and untranslated regions respectively.

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
gffcompare <- "/scr1/users/wangr5/tools/gffcompare-0.12.6/gffcompare"
target.gene <- read.table(file.path(workdir, "CDG/references/target_genes.bed"), sep = "\t", header = FALSE) %>% filter(V5 == "PIGN")
gencode.gtf <- "/scr1/users/wangr5/references/gencode.v45.annotation.gtf"
outfile1 <- file.path(workdir, "manuscript/Supplementary_Figures/Figure_S11/Figure_S11b.pdf")
outfile2 <- file.path(workdir, "manuscript/Supplementary_Figures/Figure_S11/Figure_S11c.pdf")

# =====================================================================================================================
#                                                      PANEL B
# =====================================================================================================================

# Create a temporary directory in the folder for Supplementary Figure 11
dir.create(file.path(dirname(outfile1), "tmp"), showWarnings = FALSE)

# Extract transcript annotations for target.gene from gencode.gtf
system(paste("grep \"", target.gene$V4, "\" ", gencode.gtf, " > ", file.path(dirname(outfile1), "tmp/gene.gtf"), sep = ""))

# Run stringtie on haplotype-specific BAM file for CDG-147-1
dir.create(file.path(dirname(outfile1), "tmp/CDG-147-1"), showWarnings = FALSE)
system(paste(stringtie, file.path(workdir, "CDG/CDG-147-1/RNA/stripe/target_genes/PIGN/hap1_reads.bam"), "-G", 
    file.path(dirname(outfile1), "tmp/gene.gtf"), "-o", file.path(dirname(outfile1), "tmp/CDG-147-1/output.gtf"), "-L -s 5 -c 5 -u -M 0 -l CDG-147-1"))

# Run stringtie on TEQUILA-seq BAM files for other cohort samples
cohort.samples <- read.table(file.path(workdir, "CDG/samples.txt"), sep = "\t", header = TRUE) %>% filter(Provider != "Lan Lin" &
    ID != "CDG-147-1") %>% pull(ID)
for (sample.id in cohort.samples) {
    dir.create(file.path(dirname(outfile1), "tmp", sample.id), showWarnings = FALSE)
    system(paste(samtools, "view -hb -F 256 -q 1", file.path(workdir, "CDG", sample.id, "RNA", paste(sample.id, "TEQUILA.bam", sep = "_")), 
        paste(target.gene$V1, ":", target.gene$V2+1, "-", target.gene$V3, sep = ""), " > ", file.path(dirname(outfile1), "tmp", sample.id, "input.bam")))
    system(paste(samtools, "index", file.path(dirname(outfile1), "tmp", sample.id, "input.bam")))
    system(paste(stringtie, file.path(dirname(outfile1), "tmp", sample.id, "input.bam"), "-G", file.path(dirname(outfile1), "tmp/gene.gtf"), "-o", 
        file.path(dirname(outfile1), "tmp", sample.id, "output.gtf"), "-L -s 5 -c 5 -u -M 0 -l", sample.id))
}

# Compare each sample-specific GTF file to existing GENCODE annotations
tpmDF <- tibble(Transcript_ID = character())

system(paste("cp", file.path(dirname(outfile1), "tmp/gene.gtf"), file.path(dirname(outfile1), "tmp/reference.gtf")))
current.ref <- file.path(dirname(outfile1), "tmp/reference.gtf")
for (sample.id in c("CDG-147-1", cohort.samples)) {
    system(paste(gffcompare, "-r", current.ref, "-o", file.path(dirname(outfile1), "tmp", sample.id, "gffcmp"), 
        file.path(dirname(outfile1), "tmp", sample.id, "output.gtf")), ignore.stdout = TRUE, ignore.stderr = TRUE)
    sampleDF <- read.table(file.path(dirname(outfile1), "tmp", sample.id, "output.gtf"), sep = "\t", header = FALSE, comment = "#") %>%
        filter(V3 == "transcript") %>% mutate(Transcript_ID = unlist(lapply(V9, function(x) gsub(";", "", PullFeature(x, "transcript_id")))),
        TPM = as.numeric(unlist(lapply(V9, function(x) gsub(";", "", PullFeature(x, "TPM")))))) %>% select(Transcript_ID, TPM) %>% tibble
    mappingDF <- read.table(file.path(dirname(outfile1), "tmp", sample.id, "gffcmp.tracking"), sep = "\t", header = FALSE) %>%
        separate(V3, c(NA, "Ref_Transcript_ID"), sep = "\\|") %>% separate(V5, c(NA, "Transcript_ID", NA, NA, NA, NA, NA), sep = "\\|") %>%
        filter(!(V4 %in% c("s", "x", "p", "e", "r", "u", "c"))) %>% mutate(New_ID = ifelse(V4 == "=", Ref_Transcript_ID, Transcript_ID))
    sampleDF <- mutate(sampleDF, New_ID = recode(Transcript_ID, !!!setNames(mappingDF$New_ID, mappingDF$Transcript_ID), .default = "Other")) %>%
        select(New_ID, TPM) %>% group_by(New_ID) %>% summarise(TPM = sum(TPM)) %>% ungroup %>% setNames(c("Transcript_ID", sample.id))
    tpmDF <- full_join(tpmDF, sampleDF, by = join_by(Transcript_ID)) %>% replace(is.na(.), 0)
    novel_transcripts <- mappingDF %>% filter(!(V4 %in% c("s", "x", "p", "e", "r", "u", "=", "c"))) %>% pull(New_ID)
    if(length(novel_transcripts) > 0){
        system(paste("grep -E", paste("\"", paste(paste("\"", novel_transcripts, "\"", sep = ""), collapse = "|"), "\"", sep = ""), file.path(dirname(outfile1), "tmp", sample.id, 
            "gffcmp.annotated.gtf"), ">>", current.ref))
    }
}

system(paste("tail -n +2", current.ref, ">", paste(current.ref, "tmp", sep = ".")))
system(paste("python /scr1/users/wangr5/tools/Annotate_ORF.py -i", paste(current.ref, "tmp", sep = "."), "-a", file.path(dirname(outfile1), "tmp/gene.gtf"), 
    "-f /scr1/users/wangr5/references/GRCh38.primary_assembly.genome.fa", "-o", file.path(dirname(outfile1), "tmp/reference.updated.gtf")))

# Convert TPM matrix into a proportion matrix and identify transcripts that:
#   * Have an isoform-level proportion of at least 5% in sample of interest
#   * Have an isoform-level proportion of at least 5% in at least five samples among remaining cohort individuals
propMatrix <- bind_cols(tpmDF[,1], sweep(tpmDF[,-1], 2, colSums(tpmDF[,-1]), `/`))
keepTranscripts <- unique(c(propMatrix$Transcript_ID[propMatrix[,2] >= 0.05],
    propMatrix$Transcript_ID[rowSums(propMatrix[,-c(1,2)] >= 0.05) >= 5]))
keepTranscripts <- keepTranscripts[!grepl("Other", keepTranscripts)]
propMatrix$Transcript_ID[!(propMatrix$Transcript_ID %in% keepTranscripts)] <- "Other"
propMatrix <- propMatrix %>% group_by(Transcript_ID) %>% summarise(across(everything(), sum)) %>% ungroup

# Rescale features in gtfDF
gtfDF <- read.table(file.path(dirname(outfile1), "tmp/reference.updated.gtf"), sep = "\t", header = FALSE, comment = "#") %>%
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
    paste("NovelTx", 1:length(unique(gtfDF$Transcript_ID[!grepl("ENST", gtfDF$Transcript_ID)])), sep = "."), unique(gtfDF$Transcript_ID[!grepl("ENST", gtfDF$Transcript_ID)]))))

newTxAssign <- setNames(seq(length(keepTranscripts),1), c(3,9,10,4,6,8,5,2,15,13,12,11,14,7,1))
intronDF <- mutate(intronDF, Transcript_Number = recode(Transcript_Number, !!!newTxAssign))
utrDF <- mutate(utrDF, Transcript_Number = recode(Transcript_Number, !!!newTxAssign))
cdsDF <- mutate(cdsDF, Transcript_Number = recode(Transcript_Number, !!!newTxAssign))
labelDF <- mutate(labelDF, Transcript_Number = recode(Transcript_Number, !!!newTxAssign))

palette <- setNames(c("#1C5BA6", "#357EB9", "#5A9ECC", "#8EBEDA", "#BAD2EB", "#6BC2B7", "#99D7A6", "#C2E8B9", "#FCB93F", "#FAAC90",
    "#F97D5F", "#F7523A", "#E82322", "#8C86BC", "#CFD1E6"), seq(length(keepTranscripts), 1))
p1 <- ggplot() + geom_rect(data = utrDF, fill = "white", xmin = 1 - utrDF$V4, xmax = 1 - utrDF$V5, ymin = utrDF$Transcript_Number - 0.25,
    ymax = utrDF$Transcript_Number + 0.25, color = "black", linewidth = 0.5) + geom_rect(data = cdsDF %>% mutate(Transcript_Number = factor(Transcript_Number)), 
    aes(fill = Transcript_Number), xmin = 1 - cdsDF$V4, xmax = 1 - cdsDF$V5, ymin = cdsDF$Transcript_Number - 0.25, ymax = cdsDF$Transcript_Number + 0.25, 
    color = "black", linewidth = 0.5) + geom_segment(data = intronDF, x = 1 - intronDF$V4, xend = 1 - intronDF$V5, y = intronDF$Transcript_Number, 
    yend = intronDF$Transcript_Number, linewidth = 0.5, color = "black") + theme_classic() + geom_text(data = labelDF, aes(x = -0.05, y = Transcript_Number, 
    label = Transcript_ID), size = 6*5/14, color = "black", hjust = 1) + coord_cartesian(xlim = c(-0.25, 1), ylim = c(0.5, length(keepTranscripts) + 0.5), 
    clip = "off") + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + scale_fill_manual(values = rev(palette))

ggsave(outfile1, plot = p1, width = 3, height = 2.5)

# =====================================================================================================================
#                                                      PANEL C
# =====================================================================================================================

propDF <- gather(propMatrix, "Sample_ID", "Proportion", -Transcript_ID)
sampleOrder <- propDF %>% filter(Transcript_ID == "CDG-147-1.1.3") %>% arrange(Proportion) %>% pull(Sample_ID)
propDF$Sample_ID <- factor(propDF$Sample_ID, levels = sampleOrder)

palette <- setNames(c("#1C5BA6", "#357EB9", "#5A9ECC", "#8EBEDA", "#BAD2EB", "#6BC2B7", "#99D7A6", "#C2E8B9", "#FCB93F", "#FAAC90",
    "#F97D5F", "#F7523A", "#E82322", "#8C86BC", "#CFD1E6", "#BDBDBD"), c(recode(rev(keepTranscripts)[as.integer(names(newTxAssign))], 
    !!!setNames(paste("NovelTx", 1:length(unique(gtfDF$Transcript_ID[!grepl("ENST", gtfDF$Transcript_ID)])), sep = "."), 
    unique(gtfDF$Transcript_ID[!grepl("ENST", gtfDF$Transcript_ID)]))), "Other"))
p2 <- ggplot(propDF %>% mutate(Transcript_ID = factor(recode(Transcript_ID, !!!setNames(paste("NovelTx", 1:length(unique(gtfDF$Transcript_ID[!grepl("ENST", 
    gtfDF$Transcript_ID)])), sep = "."), unique(gtfDF$Transcript_ID[!grepl("ENST", gtfDF$Transcript_ID)]))), levels = c(recode(rev(keepTranscripts)[as.integer(names(newTxAssign))], 
    !!!setNames(paste("NovelTx", 1:length(unique(gtfDF$Transcript_ID[!grepl("ENST", gtfDF$Transcript_ID)])), sep = "."), 
    unique(gtfDF$Transcript_ID[!grepl("ENST", gtfDF$Transcript_ID)]))), "Other"))), aes(x = Proportion, y = Sample_ID, fill = Transcript_ID)) + 
    geom_bar(stat = "identity", position = "stack", color = NA) + theme_classic() + xlab("Isoform proportion (PIGN)") + 
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_line(color = "black", 
    linewidth = 0.25), axis.text.x = element_text(color = "black", size = 6), axis.title.x = element_text(color = "black", size = 7), legend.text = 
    element_text(color = "black", size = 6), legend.title = element_blank(), legend.key.size = unit(0.3, "cm"),
    plot.margin = margin(r = 20, unit = "pt"), legend.position = "bottom") + scale_fill_manual(values = palette) + guides(fill = guide_legend(ncol = 2))

ggsave(outfile2, plot = p2, width = 2.5, height = 4.5)

# Remove intermediate files
system(paste("rm -rf", file.path(dirname(outfile2), "tmp")))
