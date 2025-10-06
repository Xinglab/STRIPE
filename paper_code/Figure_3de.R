#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.10.06
# Figure 3d,e

# (d) Full-length structures and (e) isoform-level proportions of ATP5MK transcripts detected from TEQUILA-seq data of 
# individual PMD-P01 and other cohort individuals. Shaded and unshaded regions within transcript structures represent 
# putative coding sequences and untranslated regions respectively. Transcripts featuring exon 3 skipping, exons 2 and 3
# skipping, or usage of a cryptic donor splice site within exon 3, with an isoform-level proportion ≥ 5% in individual 
# PMD-P01 or transcripts with an average isoform-level proportion ≥ 10% across other cohort individuals are displayed in 
# (d, e). 

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
target.gene <- read.table(file.path(workdir, "PMD/references/target_genes.bed"), sep = "\t", header = FALSE) %>% filter(V5 == "ATP5MK")
gencode.gtf <- "/scr1/users/wangr5/references/gencode.v45.annotation.gtf"

# Create a temporary directory
dir.create(file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp"), showWarnings = FALSE)

# Extract transcript annotations for target.gene from gencode.gtf
system(paste("grep \"", target.gene$V4, "\" ", gencode.gtf, " > ", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/gene.gtf"), sep = ""))

# Read in mapping file describing new IDs for cohort samples
sample.map <- read.table(file.path(workdir, "manuscript/Revisions/20250228/sample.map"), sep = "\t", header = FALSE) %>% 
    select(V2, V1) %>% tibble::deframe()

# Run stringtie on TEQUILA-seq BAM file for PMD-P01
# Use "-E 10" during StringTie run because cryptic donor splice site is being "overcorrected" back to canonical donor splice site
dir.create(file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp", sample.map["PMD-P01"]), showWarnings = FALSE)
system(paste(samtools, "view -hb -F 256 -q 1", file.path(workdir, "PMD", sample.map["PMD-P01"], "RNA", paste(sample.map["PMD-P01"], "TEQUILA.bam", sep = "_")), 
    paste(target.gene$V1, ":", target.gene$V2+1, "-", target.gene$V3, sep = ""), " > ", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp", 
    sample.map["PMD-P01"], "input.bam")))
system(paste(samtools, "index", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp", sample.map["PMD-P01"], "input.bam")))
system(paste(stringtie, file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp", sample.map["PMD-P01"], "input.bam"), "-G",
    file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/gene.gtf"), "-E 10", "-o", file.path(workdir, 
    "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp", sample.map["PMD-P01"], "output.gtf"), "-L -s 5 -c 5 -u -M 0"))

# Merge output.gtf with gene.gtf via gffcompare
system(paste(gffcompare, "-M -D -T -C -X", "-o", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/merged"), file.path(workdir, 
    "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp", sample.map["PMD-P01"], "output.gtf"), file.path(workdir, 
    "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/gene.gtf")))

# Used the following Bash code to run SQANTI3 on merged.combined.gtf 
if (FALSE) {
    "/scr1/users/wangr5/tools/SQANTI3-5.2.2/sqanti3_qc.py" --force_id_ignore --skipORF -o sqanti3 \
        -d "/mnt/isilon/lin_lab_share/STRIPE/manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp" --report skip \
        "/mnt/isilon/lin_lab_share/STRIPE/manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/merged.combined.gtf" \
        "/mnt/isilon/lin_lab_share/STRIPE/manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/gene.gtf" \
        "/scr1/users/wangr5/references/GRCh38.primary_assembly.genome.fa"
}

system(paste("python /scr1/users/wangr5/tools/Annotate_ORF.py -i", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/merged.combined.gtf"), 
    "-a", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/gene.gtf"), "-f /scr1/users/wangr5/references/GRCh38.primary_assembly.genome.fa", 
    "-o", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/merged.updated.gtf")))

# Run stringtie on TEQUILA-seq BAM files for cohort individuals
cohort.samples <- read.table(file.path(workdir, "PMD/samples.txt"), sep = "\t", header = TRUE) %>% filter((Provider == "Rebecca Ganetzky") | 
    (Provider == "Marni Falk" & Status != "Diagnosed" & !grepl("-MF", ID))) %>% pull(ID)
for (sample.id in cohort.samples) {
    dir.create(file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp", sample.id), showWarnings = FALSE)
    system(paste(samtools, "view -hb -F 256 -q 1", file.path(workdir, "PMD", sample.id, "RNA", paste(sample.id, "TEQUILA.bam", sep = "_")), 
        paste(target.gene$V1, ":", target.gene$V2+1, "-", target.gene$V3, sep = ""), " > ", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp", 
        sample.id, "input.bam")))
    system(paste(samtools, "index", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp", sample.id, "input.bam")))
    system(paste(stringtie, file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp", sample.id, "input.bam"), "-G", 
        file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/merged.combined.gtf"), "-o", file.path(workdir, 
        "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp", sample.id, "output.gtf"), "-s 5 -c 5 -u -M 0 -e -B -L"))
}

# Construct count matrix across samples
outDF <- tibble(Transcript_ID = character())
for (sample.id in cohort.samples) {
    outDF <- full_join(outDF, read.table(file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp", sample.id, 
        "t_data.ctab"), sep = "\t", header = TRUE) %>% select(t_name, cov) %>% setNames(c("Transcript_ID", sample.id)), by = join_by(Transcript_ID))
}

# Rename transcripts in outDF based on SQANTI3 assignments (also only keep FSM/NIC/NNC transcripts)
sqanti3 <- read.table(file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/sqanti3_classification.txt"), 
    sep = "\t", header = TRUE) %>% filter(structural_category %in% c("full-splice_match", "novel_in_catalog", "novel_not_in_catalog")) %>% 
    select(isoform, associated_transcript)
outDF <- outDF %>% filter(Transcript_ID %in% sqanti3$isoform)
outDF$Real_ID <- plyr::mapvalues(outDF$Transcript_ID, sqanti3$isoform, sqanti3$associated_transcript, warn_missing = FALSE)
outDF$Transcript_ID <- ifelse(outDF$Real_ID != "novel", outDF$Real_ID, outDF$Transcript_ID)
outDF <- outDF %>% select(-Real_ID)

# Convert TPM matrix into a proportion matrix
propMatrix <- bind_cols(outDF[,1], sweep(outDF[,-1], 2, colSums(outDF[,-1]), `/`))
sampleTx <- propMatrix %>% select(Transcript_ID, !!!as.character(sample.map["PMD-P01"])) %>% setNames(c("Transcript_ID", "Patient")) %>% arrange(desc(Patient))
cohortTx <- propMatrix %>% select(-(!!!as.character(sample.map["PMD-P01"]))) %>% mutate(Total = rowMeans(across(where(is.numeric)))) %>% select(Transcript_ID, Total) %>% arrange(desc(Total))

# Only keep transcripts featuring exon 3 skipping, exons 2 and 3 skipping, or usage of a cryptic donor splice site within exon 3, 
# with isoform fraction >= 0.05 in PMD-P01 or average isoform fraction >= 0.1 across other cohort individuals
keepTranscripts <- c("ENST00000369815.6", "ENST00000309579.7", "TCONS_00000007", "TCONS_00000006", "TCONS_00000002")
propMatrix$Transcript_ID[!(propMatrix$Transcript_ID %in% keepTranscripts)] <- "Other"
propMatrix <- propMatrix %>% group_by(Transcript_ID) %>% summarise(across(everything(), sum)) %>% ungroup

# Rescale features in gtfDF
gtfDF <- read.table(file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/merged.updated.gtf"), sep = "\t", header = FALSE, comment = "#") %>%
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

p1 <- ggplot() + geom_rect(data = utrDF, fill = "white", xmin = 1-utrDF$V4, xmax = 1-utrDF$V5, ymin = utrDF$Transcript_Number - 0.25,
    ymax = utrDF$Transcript_Number + 0.25, color = "black", linewidth = 0.5) + geom_rect(data = cdsDF %>% mutate(Transcript_Number = factor(Transcript_Number)), 
    aes(fill = Transcript_Number), xmin = 1-cdsDF$V4, xmax = 1-cdsDF$V5, ymin = cdsDF$Transcript_Number - 0.25, ymax = cdsDF$Transcript_Number + 0.25, 
    color = "black", linewidth = 0.5) + geom_segment(data = intronDF, x = 1-intronDF$V4, xend = 1-intronDF$V5, y = intronDF$Transcript_Number, 
    yend = intronDF$Transcript_Number, linewidth = 0.5, color = "black") + theme_classic() + geom_text(data = labelDF, aes(x = -0.05, y = Transcript_Number, 
    label = Transcript_ID), size = 6*5/14, color = "black", hjust = 1) + coord_cartesian(xlim = c(-0.25, 1), ylim = c(0.5, length(keepTranscripts) + 0.5), 
    clip = "off") + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + scale_fill_manual(values = c("#CBB3D0", "#66C2A5",
    "#3188BD", "#FEAD63", "#F36C45"))

ggsave(file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/Figure_3d.pdf"), plot = p1, width = 3, height = 1.5)

propDF <- gather(propMatrix, "Sample_ID", "Proportion", -Transcript_ID) %>% mutate(Transcript_ID = factor(Transcript_ID, levels = c(keepTranscripts, "Other")))
propDF <- propDF %>% mutate(Transcript_ID = factor(recode(Transcript_ID, !!!setNames(paste("NovelTx", 1:length(unique(gtfDF$Transcript_ID[!grepl("ENST", 
    gtfDF$Transcript_ID)])), sep = "."), unique(gtfDF$Transcript_ID[!grepl("ENST", gtfDF$Transcript_ID)]))), levels = c(labelDF %>% arrange(desc(Transcript_Number)) %>%
    pull(Transcript_ID), "Other")), Group = factor(ifelse(Sample_ID == sample.map["PMD-P01"], "Patient", "Cohort individuals"), 
    levels = c("Cohort individuals", "Patient")))
summaryDF <- propDF %>% select(-Sample_ID) %>% group_by(Group, Transcript_ID) %>% summarise(Prop_Mean = mean(Proportion))

p2 <- ggplot() + geom_bar(data = summaryDF, aes(x = Group, y = Prop_Mean * 100, fill = Transcript_ID, group = Transcript_ID), linewidth = 0.5,
    stat = "identity", position = position_dodge(), width = 0.75) + geom_point(data = propDF, aes(x = Group, y = Proportion * 100,
    color = Transcript_ID), stroke = NA, size = 0.75, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2, 
    jitter.height = 0, dodge.width = 0.75)) + theme_classic() + ylab("ATP5MK isoform proportion (%)") + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title.y = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    axis.line = element_line(color = "black", linewidth = 0.25), axis.title.x = element_blank(), legend.position = "bottom", 
    legend.key.size = unit(0.3, "cm"), legend.title = element_blank(), legend.text = element_text(color = "black", size = 6)) + 
    scale_fill_manual(values = c("#3188BD", "#66C2A5", "#CBB3D0", "#FEAD63", "#F36C45", "#BEBEBE")) + scale_color_manual(values = rep("black", 6)) + 
    ylim(0, 100)

# Remove intermediate files
system(paste("rm -rf", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp")))

ggsave(file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/Figure_3e.pdf"), plot = p2, width = 3, height = 2.5)
