#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.03
# Supplementary Figure 11b

# (b) Usage frequencies for splice junctions involving TRIT1 exons 4-5 skipping and exons 2-5 skipping in individual 
# PMD-P12 (resolved by haplotype), other cohort individuals, and 504 tissue-matched GTEx controls. 

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
#   * chr1:39848098-39853969 (exons 4-5 skipping)
#   * chr1:39848098-39883317 (exons 2-5 skipping)

if (FALSE) {
    import sys

    # Establish working directories and file paths
    workdir = '/mnt/isilon/lin_lab_share/STRIPE'
    gencode = '/scr1/users/wangr5/references/gencode.v45.annotation.gtf'

    # Import functions used in Run_STRIPE.py
    sys.path.append(workdir + '/scripts/main')
    from Run_STRIPE import *

    # Read in mapping file linking patient IDs to their original IDs
    sampleMap = pd.read_csv(workdir + '/manuscript/Revisions/20250228/sample.map', sep = '\t', header = None)
    sampleMap = dict(zip(sampleMap[1], sampleMap[0]))

    _ = os.system('mkdir -p ' + workdir + '/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/tmp')
    sampleDF = pd.read_csv(workdir + '/PMD/samples.txt', sep = '\t', header = 0)
    sampleDF = sampleDF[((sampleDF['Provider'] == 'Rebecca Ganetzky') & (sampleDF['ID'] != sampleMap['PMD-P12'])) | ((sampleDF['Provider'] == 'Marni Falk') & 
        (sampleDF['Status'] != 'Diagnosed') & ~sampleDF['ID'].str.contains('-MF'))]
    targetDF = pd.read_csv(workdir + '/PMD/references/target_genes.bed', sep = '\t', header = None)
    gene_id = targetDF[targetDF[4] == 'TRIT1'].values[0][3]
    _ = os.system('grep ' + gene_id + ' ' + gencode + ' > ' + workdir + '/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/tmp/input.gtf')
    geneDF = pd.read_csv(workdir + '/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/tmp/input.gtf', sep = '\t', header = None)
    geneDF = geneDF[geneDF[2] == 'exon']
    geneDF['transcript_id'] = geneDF[8].apply(lambda x: next(iter([item.split('"')[1] for item in x.split(';') if 'transcript_id' in item]), 'NA'))
    geneDF = geneDF.groupby(['transcript_id', 0, 6])[[3, 4]].agg(lambda x: sorted(list(x))).reset_index()
    geneDF['splice_sites'] = geneDF.apply(lambda x: [site for idx in range(len(x[3])-1) for site in [x[4][idx]+1, x[3][idx+1]-1]], axis = 1)
    sites = set(geneDF.explode('splice_sites').dropna()['splice_sites'])

    output = dict()
    for sampid in sampleDF['ID']:
        junctions = dict()
        infile = workdir + '/PMD/' + sampid + '/RNA/' + sampid + '_TEQUILA.bam'
        for read in pysam.AlignmentFile(infile, 'rb').fetch(targetDF[targetDF[4] == 'TRIT1'].values[0][0],
            targetDF[targetDF[4] == 'TRIT1'].values[0][1], targetDF[targetDF[4] == 'TRIT1'].values[0][2]):
            for item in ParseCigarJunctions(read.cigarstring, read.reference_name, read.reference_start):
                if any([int(ss) in sites for ss in item.split('_')[1:]]):
                    junctions[item] = junctions.get(item, 0) + 1
        junctions = pd.DataFrame(junctions.items(), columns = ['Junction', 'Read_Counts'])
        junctions[['Chrom', 'Start', 'End']] = junctions['Junction'].str.split('_', n = 2, expand = True)
        junctions = junctions[['Chrom', 'Start', 'End', 'Read_Counts']]
        for events in [('39848098', '39853969'), ('39848098', '39883317')]:
            results = junctions[(junctions['Start'] == events[0]) | (junctions['End'] == events[1])].merge(
                pd.DataFrame([['chr1', events[0], events[1]]], columns = ['Chrom', 'Start', 'End']), how = 'outer').fillna(0)
            results['PSI'] = results['Read_Counts']/(np.nan if results['Read_Counts'].sum() < 20 else results['Read_Counts'].sum())
            output[sampid] = output.get(sampid, []) + [results[(results['Start'] == events[0]) & (results['End'] == events[1])]['PSI'].item()]

    for idx in [1, 2]:
        junctions = dict()
        infile = workdir + '/PMD/' + sampleMap['PMD-P12'] + '/RNA/stripe/target_genes/TRIT1/hap' + str(idx) + '_reads.bam'
        for read in pysam.AlignmentFile(infile, 'rb').fetch(targetDF[targetDF[4] == 'TRIT1'].values[0][0],
            targetDF[targetDF[4] == 'TRIT1'].values[0][1], targetDF[targetDF[4] == 'TRIT1'].values[0][2]):
            for item in ParseCigarJunctions(read.cigarstring, read.reference_name, read.reference_start):
                if any([int(ss) in sites for ss in item.split('_')[1:]]):
                    junctions[item] = junctions.get(item, 0) + 1
        junctions = pd.DataFrame(junctions.items(), columns = ['Junction', 'Read_Counts'])
        junctions[['Chrom', 'Start', 'End']] = junctions['Junction'].str.split('_', n = 2, expand = True)
        junctions = junctions[['Chrom', 'Start', 'End', 'Read_Counts']]
        for events in [('39848098', '39853969'), ('39848098', '39883317')]:
            results = junctions[(junctions['Start'] == events[0]) | (junctions['End'] == events[1])].merge(
                pd.DataFrame([['chr1', events[0], events[1]]], columns = ['Chrom', 'Start', 'End']), how = 'outer').fillna(0)
            results['PSI'] = results['Read_Counts']/(np.nan if results['Read_Counts'].sum() < 20 else results['Read_Counts'].sum())
            output['Haplotype_' + str(idx)] = output.get('Haplotype_' + str(idx), []) + [results[(results['Start'] == events[0]) & (results['End'] == events[1])]['PSI'].item()]

    outDF = pd.DataFrame(output.items(), columns = ['Sample_ID', 'PSI'])
    outDF[['Exon_4_5_Skipping', 'Exon_2_5_Skipping']] = pd.DataFrame(outDF.PSI.tolist())
    outDF = outDF.drop('PSI', axis = 1)
    outDF.to_csv(workdir + '/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/tmp/cohort.tsv', sep = '\t', index = False)
}

# Establish working directories and file paths
workdir <- "/mnt/isilon/lin_lab_share/STRIPE"
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/Figure_S11b.pdf")
outDF <- tibble(Group = character(), Event = character(), PSI = numeric())

# Retrieve read support and total coverage for the following splice junctions in GTEx controls:
#   * chr1:39848098-39853969 (exons 4-5 skipping)
#   * chr1:39848098-39883317 (exons 2-5 skipping)

gtex.splice <- full_join(read.table(file.path(workdir, "PMD/references/GTEx_v8/Fibroblast/target_genes.junction_counts.txt"), 
    sep = "\t", header = TRUE, check.names = FALSE) %>% separate(Name, c("Chrom", "Start", "End"), sep = "_") %>% 
    filter(Chrom == "chr1" & (Start == "39848098" | End %in% c("39853969", "39883317"))) %>% distinct(), 
    tibble(Chrom = rep("chr1", 2), Start = rep("39848098", 2), End = c("39853969", "39883317")), 
    by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)

junction.counts <- as.numeric(gtex.splice %>% filter(Start == "39848098" & End == "39853969") %>% select(-c(Chrom, Start, End)))
junction.coverage <- as.numeric(gtex.splice %>% filter(Start == "39848098" | End == "39853969") %>% select(-c(Chrom, Start, End)) %>% colSums)
outDF <- bind_rows(outDF, tibble(PSI = junction.counts/ifelse(junction.coverage < 20, NA, junction.coverage)) %>%
    mutate(Group = "GTEx v8 controls", Event = "Exons 4-5 skipping") %>% select(Group, Event, PSI))

junction.counts <- as.numeric(gtex.splice %>% filter(Start == "39848098" & End == "39883317") %>% select(-c(Chrom, Start, End)))
junction.coverage <- as.numeric(gtex.splice %>% filter(Start == "39848098" | End == "39883317") %>% select(-c(Chrom, Start, End)) %>% colSums)
outDF <- bind_rows(outDF, tibble(PSI = junction.counts/ifelse(junction.coverage < 20, NA, junction.coverage)) %>%
    mutate(Group = "GTEx v8 controls", Event = "Exons 2-5 skipping") %>% select(Group, Event, PSI))

# Retrieve read support and total coverage for the following splice junctions in cohort samples
#   * chr1:39848098-39853969 (exons 4-5 skipping)
#   * chr1:39848098-39883317 (exons 2-5 skipping)

outDF <- bind_rows(outDF, read.table(file.path(workdir, "manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/tmp/cohort.tsv"), sep = "\t", 
    header = TRUE) %>% gather("Event", "PSI", -Sample_ID) %>% mutate(Group = ifelse(grepl("Haplotype", Sample_ID), gsub("_", " ", Sample_ID), 
    "Cohort individuals"), Event = ifelse(Event == "Exon_4_5_Skipping", "Exons 4-5 skipping", "Exons 2-5 skipping"))) %>%
    select(Group, Event, PSI)
outDF$Group <- factor(outDF$Group, levels = c("Haplotype 2", "Haplotype 1", "Cohort individuals", "GTEx v8 controls"))
outDF$Event <- factor(outDF$Event, levels = c("Exons 4-5 skipping", "Exons 2-5 skipping"))
summaryDF <- outDF %>% group_by(Group, Event) %>% summarise(PSI = mean(PSI, na.rm = TRUE))

p <- ggplot() + geom_bar(data = summaryDF, aes(x = PSI * 100, y = Group, fill = Event, group = Event), stat = "identity", position = position_dodge(), 
    width = 0.75) + geom_point(data = outDF, aes(x = PSI * 100, y = Group, color = Event), stroke = NA, size = 0.75, alpha = 0.5, 
    position = position_jitterdodge(jitter.width = 0, jitter.height = 0.2, dodge.width = 0.75)) + theme_classic() + xlab("TRIT1 splicing event frequency (%)") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title.x = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    axis.line = element_line(color = "black", linewidth = 0.25), axis.title.y = element_blank(), legend.position = "bottom", legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 6), legend.key.size = unit(0.3, "cm")) + guides(fill = guide_legend(nrow = 1)) + xlim(0, 100) +
    scale_fill_manual(values = c("#F36C45", "#FEAD63")) + scale_color_manual(values = c("black", "black"))

ggsave(outfile, plot = p, width = 3, height = 2)

system(paste("rm -rf ", file.path(workdir, "manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/tmp")))
