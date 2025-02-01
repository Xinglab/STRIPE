#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.01.18
# Supplementary Figure 16b

# Usage frequencies for a splice junction involving NUBPL exon 10 skipping in individual Q1687 (resolved by haplotype), 
# other cohort individuals, and tissue-matched GTEx controls.

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

# Used the following Python code to retrieve read support and total coverage for splice junction chr14:31846592-31859117
# in cohort samples

if (FALSE) {
    import sys

    # Establish working directories and file paths
    workdir = '/mnt/isilon/lin_lab_share/STRIPE'
    gencode = '/scr1/users/wangr5/references/gencode.v45.annotation.gtf'

    # Import functions used in Run_STRIPE.py
    sys.path.append(workdir + '/scripts/main')
    from Run_STRIPE import *

    _ = os.system('mkdir -p ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/tmp')
    sampleDF = pd.read_csv(workdir + '/PMD/samples.txt', sep = '\t', header = 0)
    sampleDF = sampleDF[((sampleDF['Provider'] == 'Rebecca Ganetzky') & (sampleDF['ID'] != 'Q1687')) | ((sampleDF['Provider'] == 'Marni Falk') & 
        (sampleDF['Status'] != 'Diagnosed') & ~sampleDF['ID'].str.contains('-MF'))]
    targetDF = pd.read_csv(workdir + '/PMD/references/target_genes.bed', sep = '\t', header = None)
    gene_id = targetDF[targetDF[4] == 'NUBPL'].values[0][3]
    _ = os.system('grep ' + gene_id + ' ' + gencode + ' > ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/tmp/input.gtf')
    geneDF = pd.read_csv(workdir + '/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/tmp/input.gtf', sep = '\t', header = None)
    geneDF = geneDF[geneDF[2] == 'exon']
    geneDF['transcript_id'] = geneDF[8].apply(lambda x: next(iter([item.split('"')[1] for item in x.split(';') if 'transcript_id' in item]), 'NA'))
    geneDF = geneDF.groupby(['transcript_id', 0, 6])[[3, 4]].agg(lambda x: sorted(list(x))).reset_index()
    geneDF['splice_sites'] = geneDF.apply(lambda x: [site for idx in range(len(x[3])-1) for site in [x[4][idx]+1, x[3][idx+1]-1]], axis = 1)
    sites = set(geneDF.explode('splice_sites').dropna()['splice_sites'])

    output = dict()
    for sampid in sampleDF['ID']:
        junctions = dict()
        infile = workdir + '/PMD/' + sampid + '/RNA/' + sampid + '_TEQUILA.bam'
        for read in pysam.AlignmentFile(infile, 'rb').fetch(targetDF[targetDF[4] == 'NUBPL'].values[0][0],
            targetDF[targetDF[4] == 'NUBPL'].values[0][1], targetDF[targetDF[4] == 'NUBPL'].values[0][2]):
            for item in ParseCigarJunctions(read.cigarstring, read.reference_name, read.reference_start):
                if any([int(ss) in sites for ss in item.split('_')[1:]]):
                    junctions[item] = junctions.get(item, 0) + 1
        junctions = pd.DataFrame(junctions.items(), columns = ['Junction', 'Read_Counts'])
        junctions[['Chrom', 'Start', 'End']] = junctions['Junction'].str.split('_', n = 2, expand = True)
        junctions = junctions[['Chrom', 'Start', 'End', 'Read_Counts']]
        junctions = junctions[(junctions['Start'] == '31846592') | (junctions['End'] == '31859117')].merge(
            pd.DataFrame([['chr14', '31846592', '31859117']], columns = ['Chrom', 'Start', 'End']), how = 'outer').fillna(0)
        junctions['PSI'] = junctions['Read_Counts']/(np.nan if junctions['Read_Counts'].sum() < 20 else junctions['Read_Counts'].sum())
        output[sampid] = junctions[(junctions['Start'] == '31846592') & (junctions['End'] == '31859117')]['PSI'].item()

    for idx in [1, 2]:
        junctions = dict()
        infile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/Q1687/NUBPL/hap' + str(idx) + '_reads.bam'
        for read in pysam.AlignmentFile(infile, 'rb').fetch(targetDF[targetDF[4] == 'NUBPL'].values[0][0],
            targetDF[targetDF[4] == 'NUBPL'].values[0][1], targetDF[targetDF[4] == 'NUBPL'].values[0][2]):
            for item in ParseCigarJunctions(read.cigarstring, read.reference_name, read.reference_start):
                if any([int(ss) in sites for ss in item.split('_')[1:]]):
                    junctions[item] = junctions.get(item, 0) + 1
        junctions = pd.DataFrame(junctions.items(), columns = ['Junction', 'Read_Counts'])
        junctions[['Chrom', 'Start', 'End']] = junctions['Junction'].str.split('_', n = 2, expand = True)
        junctions = junctions[['Chrom', 'Start', 'End', 'Read_Counts']]
        junctions = junctions[(junctions['Start'] == '31846592') | (junctions['End'] == '31859117')].merge(
            pd.DataFrame([['chr14', '31846592', '31859117']], columns = ['Chrom', 'Start', 'End']), how = 'outer').fillna(0)
        junctions['PSI'] = junctions['Read_Counts']/(np.nan if junctions['Read_Counts'].sum() < 20 else junctions['Read_Counts'].sum())
        output['Haplotype_' + str(idx)] = junctions[(junctions['Start'] == '31846592') & (junctions['End'] == '31859117')]['PSI'].item()

    pd.DataFrame(output.items(), columns = ['Sample_ID', 'PSI']).to_csv(workdir + '/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/tmp/cohort.tsv', 
        sep = '\t', index = False)
}

# Establish working directories and file paths
workdir <- "/mnt/isilon/lin_lab_share/STRIPE"
outfile <- file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/Figure_S16b.pdf")

# Retrieve read support and total coverage for splice junction chr14:31846592-31859117 in GTEx controls
gtex.splice <- full_join(read.table(file.path(workdir, "PMD/references/GTEx_v8/Fibroblast/target_genes.junction_counts.txt"), 
    sep = "\t", header = TRUE, check.names = FALSE) %>% separate(Name, c("Chrom", "Start", "End"), sep = "_") %>% 
    filter(Chrom == "chr14" & (Start == "31846592" | End == "31859117")) %>% distinct(), tibble(Chrom = "chr14", 
    Start = "31846592", End = "31859117"), by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)
junction.counts <- as.numeric(gtex.splice %>% filter(Chrom == "chr14" & Start == "31846592" & End == "31859117") %>% 
    select(-c(Chrom, Start, End)))
junction.coverage <- as.numeric(gtex.splice %>% select(-c(Chrom, Start, End)) %>% colSums)
outDF <- tibble(PSI = junction.counts/ifelse(junction.coverage < 20, NA, junction.coverage)) %>% mutate(Group = "GTEx v8 controls")

# Retrieve read support and total coverage for splice junction chr14:31846592-31859117 in cohort samples
outDF <- bind_rows(outDF, read.table(file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/tmp/cohort.tsv"), sep = "\t", 
    header = TRUE) %>% mutate(Group = ifelse(grepl("Haplotype", Sample_ID), gsub("_", " ", Sample_ID), "Cohort individuals")) %>%
    select(PSI, Group))
outDF$Group <- factor(outDF$Group, levels = c("Haplotype 2", "Haplotype 1", "Cohort individuals", "GTEx v8 controls"))
summaryDF <- outDF %>% group_by(Group) %>% summarise(PSI = mean(PSI, na.rm = TRUE))

p <- ggplot() + geom_bar(data = summaryDF, aes(x = PSI * 100, y = Group, fill = Group), stat = "identity", width = 0.75) + geom_point(data = outDF,
    aes(x = PSI * 100, y = Group), color = "black", stroke = NA, size = 0.75, alpha = 0.5, position = position_jitter(width = 0, 
    height = 0.2)) + theme_classic() + xlab("NUBPL exon 10 skipping (%)") + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title.x = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    axis.line = element_line(color = "black", linewidth = 0.25), axis.title.y = element_blank(), legend.position = "none") + 
    scale_fill_manual(values = c("#90AFB1", "#EC7B02", "#94C1B6", "#5A8BA8")) + xlim(0, 100)

ggsave(outfile, plot = p, width = 3, height = 1.25)

system(paste("rm -rf ", file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/tmp")))
