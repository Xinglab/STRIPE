#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.03
# Figure 3c

# (c) Usage frequencies for splice junctions involving ATP5MK exon 3 skipping, exons 2-3 skipping, and activation of a 
# cryptic donor splice site within exon 3 for individual PMD-P01, other cohort individuals, and 504 tissue-matched GTEx 
# controls. 

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

# Used the following Python code to retrieve read support and total coverage for the following splice junctions in cohort
# samples: 
#   * chr10:103,392,284-103,395,745 (exon 3 skipping)
#   * chr10:103,392,284-103,396,408 (exons 2-3 skipping)
#   * chr10:103,392,284-103,392,381 (cryptic donor splice site, exon 3)

if (FALSE) {
    import sys

    # Establish working directories and file paths
    workdir = '/mnt/isilon/lin_lab_share/STRIPE'
    gencode = '/scr1/users/wangr5/references/gencode.v45.annotation.gtf'

    # Import functions used in Run_STRIPE.py
    sys.path.append(workdir + '/scripts/main')
    from Run_STRIPE import *

    _ = os.system('mkdir -p ' + workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp')
    sampleDF = pd.read_csv(workdir + '/PMD/samples.txt', sep = '\t', header = 0)
    sampleDF = sampleDF[(sampleDF['Provider'] == 'Rebecca Ganetzky') | ((sampleDF['Provider'] == 'Marni Falk') & 
        (sampleDF['Status'] != 'Diagnosed') & ~sampleDF['ID'].str.contains('-MF'))]
    targetDF = pd.read_csv(workdir + '/PMD/references/target_genes.bed', sep = '\t', header = None)
    gene_id = targetDF[targetDF[4] == 'ATP5MK'].values[0][3]
    _ = os.system('grep ' + gene_id + ' ' + gencode + ' > ' + workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/input.gtf')
    geneDF = pd.read_csv(workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/input.gtf', sep = '\t', header = None)
    geneDF = geneDF[geneDF[2] == 'exon']
    geneDF['transcript_id'] = geneDF[8].apply(lambda x: next(iter([item.split('"')[1] for item in x.split(';') if 'transcript_id' in item]), 'NA'))
    geneDF = geneDF.groupby(['transcript_id', 0, 6])[[3, 4]].agg(lambda x: sorted(list(x))).reset_index()
    geneDF['splice_sites'] = geneDF.apply(lambda x: [site for idx in range(len(x[3])-1) for site in [x[4][idx]+1, x[3][idx+1]-1]], axis = 1)
    sites = set(geneDF.explode('splice_sites').dropna()['splice_sites'])

    output = dict()
    for sampid in sampleDF['ID']:
        junctions = dict()
        infile = workdir + '/PMD/' + sampid + '/RNA/' + sampid + '_TEQUILA.bam'
        for read in pysam.AlignmentFile(infile, 'rb').fetch(targetDF[targetDF[4] == 'ATP5MK'].values[0][0],
            targetDF[targetDF[4] == 'ATP5MK'].values[0][1], targetDF[targetDF[4] == 'ATP5MK'].values[0][2]):
            for item in ParseCigarJunctions(read.cigarstring, read.reference_name, read.reference_start):
                if any([int(ss) in sites for ss in item.split('_')[1:]]):
                    junctions[item] = junctions.get(item, 0) + 1
        junctions = pd.DataFrame(junctions.items(), columns = ['Junction', 'Read_Counts'])
        junctions[['Chrom', 'Start', 'End']] = junctions['Junction'].str.split('_', n = 2, expand = True)
        junctions = junctions[['Chrom', 'Start', 'End', 'Read_Counts']]
        for events in [('103392284', '103395745'), ('103392284', '103396408'), ('103392284', '103392381')]:
            results = junctions[(junctions['Start'] == events[0]) | (junctions['End'] == events[1])].merge(
                pd.DataFrame([['chr10', events[0], events[1]]], columns = ['Chrom', 'Start', 'End']), how = 'outer').fillna(0)
            results['PSI'] = results['Read_Counts']/(np.nan if results['Read_Counts'].sum() < 20 else results['Read_Counts'].sum())
            output[sampid] = output.get(sampid, []) + [results[(results['Start'] == events[0]) & (results['End'] == events[1])]['PSI'].item()]

    outDF = pd.DataFrame(output.items(), columns = ['Sample_ID', 'PSI'])
    outDF[['Exon_3_Skipping', 'Exon_2_3_Skipping', 'Cryptic_SS']] = pd.DataFrame(outDF.PSI.tolist())
    outDF = outDF.drop('PSI', axis = 1)
    outDF.to_csv(workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/cohort.tsv', sep = '\t', index = False)
}

# Establish working directories and file paths
workdir <- "/mnt/isilon/lin_lab_share/STRIPE"
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/Figure_3c.pdf")
outDF <- tibble(Group = character(), Event = character(), PSI = numeric())

# Read in mapping file describing new IDs for cohort samples
sample.map <- read.table(file.path(workdir, "manuscript/Revisions/20250228/sample.map"), sep = "\t", header = FALSE) %>% 
    select(V2, V1) %>% tibble::deframe()

# Retrieve read support and total coverage for the following splice junctions in GTEx controls:
#   * chr10:103,392,284-103,395,745 (exon 3 skipping)
#   * chr10:103,392,284-103,396,408 (exons 2-3 skipping)
#   * chr10:103,392,284-103,392,381 (cryptic donor splice site, exon 3)

gtex.splice <- full_join(read.table(file.path(workdir, "PMD/references/GTEx_v8/Fibroblast/target_genes.junction_counts.txt"), 
    sep = "\t", header = TRUE, check.names = FALSE) %>% separate(Name, c("Chrom", "Start", "End"), sep = "_") %>% 
    filter(Chrom == "chr10" & (Start == "103392284" | End %in% c("103395745", "103396408", "103392381"))) %>% distinct(), 
    tibble(Chrom = rep("chr10", 3), Start = rep("103392284", 3), End = c("103395745", "103396408", "103392381")), 
    by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)

junction.counts <- as.numeric(gtex.splice %>% filter(Start == "103392284" & End == "103395745") %>% select(-c(Chrom, Start, End)))
junction.coverage <- as.numeric(gtex.splice %>% filter(Start == "103392284" | End == "103395745") %>% select(-c(Chrom, Start, End)) %>% colSums)
outDF <- bind_rows(outDF, tibble(PSI = junction.counts/ifelse(junction.coverage < 20, NA, junction.coverage)) %>%
    mutate(Group = "GTEx v8 controls", Event = "Exon 3 skipping") %>% select(Group, Event, PSI))

junction.counts <- as.numeric(gtex.splice %>% filter(Start == "103392284" & End == "103396408") %>% select(-c(Chrom, Start, End)))
junction.coverage <- as.numeric(gtex.splice %>% filter(Start == "103392284" | End == "103396408") %>% select(-c(Chrom, Start, End)) %>% colSums)
outDF <- bind_rows(outDF, tibble(PSI = junction.counts/ifelse(junction.coverage < 20, NA, junction.coverage)) %>%
    mutate(Group = "GTEx v8 controls", Event = "Exons 2-3 skipping") %>% select(Group, Event, PSI))

junction.counts <- as.numeric(gtex.splice %>% filter(Start == "103392284" & End == "103392381") %>% select(-c(Chrom, Start, End)))
junction.coverage <- as.numeric(gtex.splice %>% filter(Start == "103392284" | End == "103392381") %>% select(-c(Chrom, Start, End)) %>% colSums)
outDF <- bind_rows(outDF, tibble(PSI = junction.counts/ifelse(junction.coverage < 20, NA, junction.coverage)) %>%
    mutate(Group = "GTEx v8 controls", Event = "Cryptic donor splice site (exon 3)") %>% select(Group, Event, PSI))

# Retrieve read support and total coverage for the following splice junctions in cohort samples
#   * chr10:103,392,284-103,395,745 (exon 3 skipping)
#   * chr10:103,392,284-103,396,408 (exons 2-3 skipping)
#   * chr10:103,392,284-103,392,381 (cryptic donor splice site, exon 3)

outDF <- bind_rows(outDF, read.table(file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/cohort.tsv"), sep = "\t", 
    header = TRUE) %>% gather("Event", "PSI", -Sample_ID) %>% mutate(Group = ifelse(Sample_ID == sample.map["PMD-P01"], "Patient", "Cohort individuals"),
    Event = ifelse(Event == "Exon_3_Skipping", "Exon 3 skipping", ifelse(Event == "Exon_2_3_Skipping", "Exons 2-3 skipping", 
    "Cryptic donor splice site (exon 3)"))) %>% select(Group, Event, PSI))
outDF$Group <- factor(outDF$Group, levels = c("Patient", "Cohort individuals", "GTEx v8 controls"))
outDF$Event <- factor(outDF$Event, levels = c("Exon 3 skipping", "Exons 2-3 skipping", "Cryptic donor splice site (exon 3)"))
summaryDF <- outDF %>% group_by(Group, Event) %>% summarise(PSI = mean(PSI, na.rm = TRUE))

p <- ggplot() + geom_bar(data = summaryDF, aes(x = PSI * 100, y = Group, fill = Event, group = Event), stat = "identity", position = position_dodge(), 
    width = 0.75) + geom_point(data = outDF, aes(x = PSI * 100, y = Group, color = Event), stroke = NA, size = 0.75, alpha = 0.5, 
    position = position_jitterdodge(jitter.width = 0, jitter.height = 0.2, dodge.width = 0.75)) + theme_classic() + xlab("ATP5MK splicing event frequency (%)") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title.x = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    axis.line = element_line(color = "black", linewidth = 0.25), axis.title.y = element_blank(), legend.position = "bottom", legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 6), legend.key.size = unit(0.3, "cm")) + guides(fill = guide_legend(ncol = 1)) + xlim(0, 100) +
    scale_fill_manual(values = c("#F36C45", "#FEAD63", "#CBB3D0")) + scale_color_manual(values = c("black", "black", "black"))

ggsave(outfile, plot = p, width = 3, height = 2)

system(paste("rm -rf ", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp")))
