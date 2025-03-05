#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.05
# Figure 6b

# (b) Usage frequencies for splice junctions involving MTFMT exon 4 skipping and inclusion of a 53 bp pseudoexon within MTFMT 
# intron 6 in individual PMD-P23 (resolved by haplotype), other cohort individuals, and 504 tissue-matched GTEx controls. 

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
#   * chr15:65020273-65023671 (exon 4 skipping)
#   * chr15:65008954-65016435 (psuedoexon inclusion, intron 6, upstream junction)
#   * chr15:65006192-65008900 (pseudoexon inclusion, intron 6, downstream junction)

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

    _ = os.system('mkdir -p ' + workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/tmp')
    sampleDF = pd.read_csv(workdir + '/PMD/samples.txt', sep = '\t', header = 0)
    sampleDF = sampleDF[(sampleDF['Provider'] == 'Rebecca Ganetzky') | ((sampleDF['Provider'] == 'Marni Falk') & 
        (sampleDF['Status'] != 'Diagnosed') & ~sampleDF['ID'].str.contains('-MF') & (sampleDF['ID'] != sampleMap['PMD-P23']))]
    targetDF = pd.read_csv(workdir + '/PMD/references/target_genes.bed', sep = '\t', header = None)
    gene_id = targetDF[targetDF[4] == 'MTFMT'].values[0][3]
    _ = os.system('grep ' + gene_id + ' ' + gencode + ' > ' + workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/tmp/input.gtf')
    geneDF = pd.read_csv(workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/tmp/input.gtf', sep = '\t', header = None)
    geneDF = geneDF[geneDF[2] == 'exon']
    geneDF['transcript_id'] = geneDF[8].apply(lambda x: next(iter([item.split('"')[1] for item in x.split(';') if 'transcript_id' in item]), 'NA'))
    geneDF = geneDF.groupby(['transcript_id', 0, 6])[[3, 4]].agg(lambda x: sorted(list(x))).reset_index()
    geneDF['splice_sites'] = geneDF.apply(lambda x: [site for idx in range(len(x[3])-1) for site in [x[4][idx]+1, x[3][idx+1]-1]], axis = 1)
    sites = set(geneDF.explode('splice_sites').dropna()['splice_sites'])

    output = dict()
    for sampid in sampleDF['ID']:
        junctions = dict()
        infile = workdir + '/PMD/' + sampid + '/RNA/' + sampid + '_TEQUILA.bam'
        for read in pysam.AlignmentFile(infile, 'rb').fetch(targetDF[targetDF[4] == 'MTFMT'].values[0][0],
            targetDF[targetDF[4] == 'MTFMT'].values[0][1], targetDF[targetDF[4] == 'MTFMT'].values[0][2]):
            for item in ParseCigarJunctions(read.cigarstring, read.reference_name, read.reference_start):
                if any([int(ss) in sites for ss in item.split('_')[1:]]):
                    junctions[item] = junctions.get(item, 0) + 1
        junctions = pd.DataFrame(junctions.items(), columns = ['Junction', 'Read_Counts'])
        junctions[['Chrom', 'Start', 'End']] = junctions['Junction'].str.split('_', n = 2, expand = True)
        junctions = junctions[['Chrom', 'Start', 'End', 'Read_Counts']]
        for events in [('65020273', '65023671'), ('65008954', '65016435'), ('65006192', '65008900')]:
            results = junctions[(junctions['Start'] == events[0]) | (junctions['End'] == events[1])].merge(
                pd.DataFrame([['chr15', events[0], events[1]]], columns = ['Chrom', 'Start', 'End']), how = 'outer').fillna(0)
            results['PSI'] = results['Read_Counts']/(np.nan if results['Read_Counts'].sum() < 20 else results['Read_Counts'].sum())
            output[sampid] = output.get(sampid, []) + [results[(results['Start'] == events[0]) & (results['End'] == events[1])]['PSI'].item()]

    for idx in [1, 2]:
        junctions = dict()
        infile = workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/hap' + str(idx) + '_reads.bam'
        for read in pysam.AlignmentFile(infile, 'rb').fetch(targetDF[targetDF[4] == 'MTFMT'].values[0][0],
            targetDF[targetDF[4] == 'MTFMT'].values[0][1], targetDF[targetDF[4] == 'MTFMT'].values[0][2]):
            for item in ParseCigarJunctions(read.cigarstring, read.reference_name, read.reference_start):
                if any([int(ss) in sites for ss in item.split('_')[1:]]):
                    junctions[item] = junctions.get(item, 0) + 1
        junctions = pd.DataFrame(junctions.items(), columns = ['Junction', 'Read_Counts'])
        junctions[['Chrom', 'Start', 'End']] = junctions['Junction'].str.split('_', n = 2, expand = True)
        junctions = junctions[['Chrom', 'Start', 'End', 'Read_Counts']]
        for events in [('65020273', '65023671'), ('65008954', '65016435'), ('65006192', '65008900')]:
            results = junctions[(junctions['Start'] == events[0]) | (junctions['End'] == events[1])].merge(
                pd.DataFrame([['chr15', events[0], events[1]]], columns = ['Chrom', 'Start', 'End']), how = 'outer').fillna(0)
            results['PSI'] = results['Read_Counts']/(np.nan if results['Read_Counts'].sum() < 20 else results['Read_Counts'].sum())
            output['Haplotype_' + str(idx)] = output.get('Haplotype_' + str(idx), []) + [results[(results['Start'] == events[0]) & (results['End'] == events[1])]['PSI'].item()]

    outDF = pd.DataFrame(output.items(), columns = ['Sample_ID', 'PSI'])
    outDF[['Exon_4_Skipping', 'Pseudoexon_Upstream', 'Pseudoexon_Downstream']] = pd.DataFrame(outDF.PSI.tolist())
    outDF = outDF.drop('PSI', axis = 1)
    outDF.to_csv(workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/tmp/cohort.tsv', sep = '\t', index = False)
}

# Establish working directories and file paths
workdir <- "/mnt/isilon/lin_lab_share/STRIPE"
outfile <- file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_6/Figure_6b.pdf")
outDF <- tibble(Group = character(), Event = character(), PSI = numeric())

# Retrieve read support and total coverage for the following splice junctions in GTEx controls:
#   * chr15:65020273-65023671 (exon 4 skipping)
#   * chr15:65008954-65016435 (psuedoexon inclusion, intron 6, upstream junction)
#   * chr15:65006192-65008900 (pseudoexon inclusion, intron 6, downstream junction)

gtex.splice <- full_join(read.table(file.path(workdir, "PMD/references/GTEx_v8/Fibroblast/target_genes.junction_counts.txt"), 
    sep = "\t", header = TRUE, check.names = FALSE) %>% separate(Name, c("Chrom", "Start", "End"), sep = "_") %>% 
    filter(Chrom == "chr15" & (Start %in% c("65020273", "65008954", "65006192") | End %in% c("65023671", "65016435", "65008900"))) %>% distinct(), 
    tibble(Chrom = rep("chr15", 3), Start = c("65020273", "65008954", "65006192"), End =c("65023671", "65016435", "65008900")), 
    by = join_by(Chrom, Start, End)) %>% replace(is.na(.), 0)

junction.counts <- as.numeric(gtex.splice %>% filter(Start == "65020273" & End == "65023671") %>% select(-c(Chrom, Start, End)))
junction.coverage <- as.numeric(gtex.splice %>% filter(Start == "65020273" | End == "65023671") %>% select(-c(Chrom, Start, End)) %>% colSums)
outDF <- bind_rows(outDF, tibble(PSI = junction.counts/ifelse(junction.coverage < 20, NA, junction.coverage)) %>%
    mutate(Group = "GTEx v8 controls", Event = "Exon 4 skipping") %>% select(Group, Event, PSI))

junction.counts.upstream <- as.numeric(gtex.splice %>% filter(Start == "65008954" & End == "65016435") %>% select(-c(Chrom, Start, End)))
junction.coverage.upstream <- as.numeric(gtex.splice %>% filter(Start == "65008954" | End == "65016435") %>% select(-c(Chrom, Start, End)) %>% colSums)
junction.counts.downstream <- as.numeric(gtex.splice %>% filter(Start == "65006192" & End == "65008900") %>% select(-c(Chrom, Start, End)))
junction.coverage.downstream <- as.numeric(gtex.splice %>% filter(Start == "65006192" | End == "65008900") %>% select(-c(Chrom, Start, End)) %>% colSums)

outDF <- bind_rows(outDF, tibble(PSI = 0.5*((junction.counts.upstream/ifelse(junction.coverage.upstream < 20, NA, junction.coverage.upstream)) +
    (junction.counts.downstream/ifelse(junction.coverage.downstream < 20, NA, junction.coverage.downstream)))) %>%
    mutate(Group = "GTEx v8 controls", Event = "Pseudoexon inclusion (intron 6)") %>% select(Group, Event, PSI))

# Retrieve read support and total coverage for the following splice junctions in cohort samples
#   * chr15:65020273-65023671 (exon 4 skipping)
#   * chr15:65008954-65016435 (psuedoexon inclusion, intron 6, upstream junction)
#   * chr15:65006192-65008900 (pseudoexon inclusion, intron 6, downstream junction)

outDF <- bind_rows(outDF, read.table(file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_6/tmp/cohort.tsv"), sep = "\t", 
    header = TRUE) %>% gather("Event", "PSI", -Sample_ID) %>% mutate(Group = ifelse(grepl("Haplotype", Sample_ID), gsub("_", " ", Sample_ID), 
    "Cohort individuals"), Event = ifelse(Event == "Exon_4_Skipping", "Exon 4 skipping", "Pseudoexon inclusion (intron 6)")) %>%
    group_by(Sample_ID, Event, Group) %>% summarise(PSI = mean(PSI))) %>% select(Group, Event, PSI)

outDF$Group <- factor(outDF$Group, levels = c("Haplotype 2", "Haplotype 1", "Cohort individuals", "GTEx v8 controls"))
outDF$Event <- factor(outDF$Event, levels = c("Exon 4 skipping", "Pseudoexon inclusion (intron 6)"))
summaryDF <- outDF %>% group_by(Group, Event) %>% summarise(PSI = mean(PSI, na.rm = TRUE))

p <- ggplot() + geom_bar(data = summaryDF, aes(x = PSI * 100, y = Group, fill = Event, group = Event), stat = "identity", position = position_dodge(), 
    width = 0.75) + geom_point(data = outDF, aes(x = PSI * 100, y = Group, color = Event), stroke = NA, size = 0.75, alpha = 0.5, 
    position = position_jitterdodge(jitter.width = 0, jitter.height = 0.2, dodge.width = 0.75)) + theme_classic() + xlab("MTFMT splicing event frequency (%)") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", size = 6), 
    axis.title.x = element_text(color = "black", size = 7), axis.ticks = element_line(color = "black", linewidth = 0.25), 
    axis.line = element_line(color = "black", linewidth = 0.25), axis.title.y = element_blank(), legend.position = "bottom", legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 6), legend.key.size = unit(0.3, "cm")) + guides(fill = guide_legend(nrow = 1)) + xlim(0, 100) +
    scale_fill_manual(values = c("#FEAD63", "#CBB3D0")) + scale_color_manual(values = c("black", "black"))

ggsave(outfile, plot = p, width = 3, height = 2)

system(paste("rm -rf ", file.path(workdir, "manuscript/Revisions/20250228/Main_Figures/Figure_6/tmp")))
