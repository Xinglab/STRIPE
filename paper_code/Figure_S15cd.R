#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.01.18
# Supplementary Figure 5c,d

# (c) Full-length structures and (d) isoform-level proportions of NUBPL transcripts detected from TEQUILA-seq data of 
# individual Q1687 (resolved by haplotype) and other cohort individuals. Shaded and unshaded regions within transcript 
# structures represent putative coding sequences and untranslated regions respectively. 

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
target.gene <- read.table(file.path(workdir, "PMD/references/target_genes.bed"), sep = "\t", header = FALSE) %>% filter(V5 == "NUBPL")
gencode.gtf <- "/scr1/users/wangr5/references/gencode.v45.annotation.gtf"

# Create a temporary directory
dir.create(file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp"), showWarnings = FALSE)

# Extract transcript annotations for target.gene from gencode.gtf
system(paste("grep \"", target.gene$V4, "\" ", gencode.gtf, " > ", file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/gene.gtf"), sep = ""))

# Run stringtie on haplotype-specific BAM files of interest for Q1687
dir.create(file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/Haplotype_1"), showWarnings = FALSE)
system(paste(stringtie, file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/Q1687/NUBPL/hap1_reads.bam"), 
    "-G", file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/gene.gtf"), "-o", file.path(workdir, 
    "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/Haplotype_1/output.gtf"), "-L -s 5 -c 5 -u -M 0"))
dir.create(file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/Haplotype_2"), showWarnings = FALSE)
system(paste(stringtie, file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/Q1687/NUBPL/hap2_reads.bam"), 
    "-G", file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/gene.gtf"), "-o", file.path(workdir, 
    "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/Haplotype_2/output.gtf"), "-L -s 5 -c 5 -u -M 0"))

# Merge output.gtf files with gene.gtf via gffcompare
# Disable -X because this unintentionally removes an aberrant NUBPL transcript featuring intron 8 polyadenylation
system(paste(gffcompare, "-M -D -T -C", "-o", file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/merged"), file.path(workdir, 
    "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/Haplotype_1/output.gtf"), file.path(workdir, 
    "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/Haplotype_2/output.gtf"), file.path(workdir, 
    "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/gene.gtf")))

# Used the following Bash code to run SQANTI3 on merged.combined.gtf 
if (FALSE) {
    "/scr1/users/wangr5/tools/SQANTI3-5.2.2/sqanti3_qc.py" --force_id_ignore --skipORF -o sqanti3 \
        -d "/mnt/isilon/lin_lab_share/STRIPE/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp" --report skip \
        "/mnt/isilon/lin_lab_share/STRIPE/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/merged.combined.gtf" \
        "/mnt/isilon/lin_lab_share/STRIPE/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/gene.gtf" \
        "/scr1/users/wangr5/references/GRCh38.primary_assembly.genome.fa"
}

system(paste("python /scr1/users/wangr5/tools/Annotate_ORF.py -i", file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/merged.combined.gtf"), 
    "-a", file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/gene.gtf"), "-f /scr1/users/wangr5/references/GRCh38.primary_assembly.genome.fa", 
    "-o", file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/merged.updated.gtf")))

# Re-run stringtie on TEQUILA-seq BAM files (disable -L because it is providing quantifications that are inconsistent with raw data)
system(paste(stringtie, file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/Q1687/NUBPL/hap1_reads.bam"), 
    "-G", file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/merged.combined.gtf"), "-o", file.path(workdir, 
    "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/Haplotype_1/output.gtf"), "-s 5 -c 5 -u -M 0 -e -B"))
system(paste(stringtie, file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/Q1687/NUBPL/hap2_reads.bam"), 
    "-G", file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/merged.combined.gtf"), "-o", file.path(workdir, 
    "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/Haplotype_2/output.gtf"), "-s 5 -c 5 -u -M 0 -e -B"))

# Run stringtie on TEQUILA-seq BAM files for other cohort individuals
cohort.samples <- read.table(file.path(workdir, "PMD/samples.txt"), sep = "\t", header = TRUE) %>% filter((Provider == "Rebecca Ganetzky" & ID != "Q1687") | 
    (Provider == "Marni Falk" & Status != "Diagnosed" & !grepl("-MF", ID))) %>% pull(ID)
for (sample.id in cohort.samples) {
    dir.create(file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp", sample.id), showWarnings = FALSE)
    system(paste(samtools, "view -hb -F 256 -q 1", file.path(workdir, "PMD", sample.id, "RNA", paste(sample.id, "TEQUILA.bam", sep = "_")), 
        paste(target.gene$V1, ":", target.gene$V2+1, "-", target.gene$V3, sep = ""), " > ", file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp", 
        sample.id, "input.bam")))
    system(paste(samtools, "index", file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp", sample.id, "input.bam")))
    system(paste(stringtie, file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp", sample.id, "input.bam"), "-G", 
        file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/merged.combined.gtf"), "-o", file.path(workdir, 
        "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp", sample.id, "output.gtf"), "-s 5 -c 5 -u -M 0 -e -B"))
}

# Construct FPKM matrix across samples
outDF <- tibble(Transcript_ID = character())
for (sample.id in c("Haplotype_1", "Haplotype_2", cohort.samples)) {
    outDF <- full_join(outDF, read.table(file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp", sample.id, 
        "t_data.ctab"), sep = "\t", header = TRUE) %>% select(t_name, FPKM) %>% setNames(c("Transcript_ID", sample.id)), by = join_by(Transcript_ID))
}

# Rename transcripts in outDF based on SQANTI3 assignments
sqanti3 <- read.table(file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/sqanti3_classification.txt"), 
    sep = "\t", header = TRUE) %>% filter(structural_category == "full-splice_match") %>% select(isoform, associated_transcript)
outDF$Real_ID <- plyr::mapvalues(outDF$Transcript_ID, sqanti3$isoform, sqanti3$associated_transcript, warn_missing = FALSE)
outDF$Transcript_ID <- ifelse(outDF$Real_ID != "novel", outDF$Real_ID, outDF$Transcript_ID)
outDF <- outDF %>% select(-Real_ID)

# Convert TPM matrix into a proportion matrix
propMatrix <- bind_cols(outDF[,1], sweep(outDF[,-1], 2, colSums(outDF[,-1]), `/`))
sampleTx1 <- propMatrix[,c(1,2)] %>% arrange(desc(`Haplotype_1`))
sampleTx2 <- propMatrix[,c(1,3)] %>% arrange(desc(`Haplotype_2`))
cohortTx <- propMatrix[,-c(2,3)] %>% mutate(Total = rowSums(across(where(is.numeric)))) %>% select(Transcript_ID, Total) %>% arrange(desc(Total))
keepTranscripts <- c("ENST00000281081.12", "TCONS_00000004", "TCONS_00000008")
propMatrix$Transcript_ID[!(propMatrix$Transcript_ID %in% keepTranscripts)] <- "Other"
propMatrix <- propMatrix %>% group_by(Transcript_ID) %>% summarise(across(everything(), sum)) %>% ungroup

# Rescale features in gtfDF
gtfDF <- read.table(file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp/merged.updated.gtf"), sep = "\t", header = FALSE, comment = "#") %>%
    mutate(Transcript_ID = unlist(lapply(V9, function(x) gsub(";", "", PullFeature(x, "transcript_id"))))) %>% mutate(Real_ID = plyr::mapvalues(Transcript_ID, 
    sqanti3$isoform, sqanti3$associated_transcript, warn_missing = FALSE)) %>% mutate(Transcript_ID = ifelse(Real_ID == "novel", Transcript_ID, Real_ID)) %>%
    select(-Real_ID) %>% filter(Transcript_ID %in% keepTranscripts) %>% select(Transcript_ID, V3, V4, V5) %>% arrange(Transcript_ID, V4, V5) 
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

p1 <- ggplot() + geom_rect(data = utrDF, fill = "white", xmin = utrDF$V4, xmax = utrDF$V5, ymin = utrDF$Transcript_Number - 0.25,
    ymax = utrDF$Transcript_Number + 0.25, color = "black", linewidth = 0.5) + geom_rect(data = cdsDF %>% mutate(Transcript_Number = factor(Transcript_Number)), 
    aes(fill = Transcript_Number), xmin = cdsDF$V4, xmax = cdsDF$V5, ymin = cdsDF$Transcript_Number - 0.25, ymax = cdsDF$Transcript_Number + 0.25, 
    color = "black", linewidth = 0.5) + geom_segment(data = intronDF, x = intronDF$V4, xend = intronDF$V5, y = intronDF$Transcript_Number, 
    yend = intronDF$Transcript_Number, linewidth = 0.5, color = "black") + theme_classic() + geom_text(data = labelDF, aes(x = -0.05, y = Transcript_Number, 
    label = Transcript_ID), size = 6*5/14, color = "black", hjust = 1) + coord_cartesian(xlim = c(-0.25, 1), ylim = c(0.5, length(keepTranscripts) + 0.5), 
    clip = "off") + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + scale_fill_manual(values = c("#CBB3D0", "#FCAE62", "#3188BD"))

ggsave(file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/Figure_S15c.pdf"), plot = p1, width = 3, height = 1.25)

propDF <- gather(propMatrix, "Sample_ID", "Proportion", -Transcript_ID) %>% mutate(Transcript_ID = factor(Transcript_ID, levels = c(keepTranscripts, "Other")))
propDF <- propDF %>% mutate(Transcript_ID = factor(recode(Transcript_ID, !!!setNames(paste("NovelTx", 1:length(unique(gtfDF$Transcript_ID[!grepl("ENST", 
    gtfDF$Transcript_ID)])), sep = "."), unique(gtfDF$Transcript_ID[!grepl("ENST", gtfDF$Transcript_ID)]))), levels = c(labelDF %>% arrange(desc(Transcript_Number)) %>%
    pull(Transcript_ID),"Other")), Group = factor(ifelse(grepl("Haplotype", Sample_ID), gsub("_", " ", Sample_ID), "Cohort individuals"), 
    levels = c("Cohort individuals", "Haplotype 1", "Haplotype 2")))
summaryDF <- propDF %>% select(-Sample_ID) %>% group_by(Group, Transcript_ID) %>% summarise(Prop_Mean = mean(Proportion, na.rm = T))

p2 <- ggplot() + geom_bar(data = summaryDF, aes(x = Group, y = Prop_Mean * 100, fill = Transcript_ID, group = Transcript_ID), linewidth = 0.5,
    stat = "identity", position = position_dodge(), width = 0.75) + geom_point(data = propDF, aes(x = Group, y = Proportion * 100,
    color = Transcript_ID), stroke = NA, size = 0.75, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2, 
    jitter.height = 0, dodge.width = 0.75)) + theme_classic() + ylab("NUBPL isoform proportion (%)") + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title.y = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    axis.line = element_line(color = "black", linewidth = 0.25), axis.title.x = element_blank(), legend.position = "bottom", 
    legend.key.size = unit(0.3, "cm"), legend.title = element_blank(), legend.text = element_text(color = "black", size = 6)) + 
    scale_fill_manual(values = c("#3188BD", "#FCAE62", "#CBB3D0", "#BEBEBE")) + scale_color_manual(values = c("black", "black", "black", "black")) + 
    ylim(0, 100) + guides(fill = guide_legend(nrow = 2))

# Remove intermediate files
system(paste("rm -rf", file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/tmp")))

ggsave(file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S15/Figure_S15d.pdf"), plot = p2, width = 3, height = 2.25)
