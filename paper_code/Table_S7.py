#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2025.01.03
Supplementary Table 7

Usage frequencies for splice junctions linked to haplotypes carrying pathogenic variants based on patient TEQUILA-seq reads.
P values were computed using an outlier splice junction test, which assesses whether the usage frequency of a splice junction 
in a sample is consistent with usage frequencies observed in tissue-matched GTEx controls (Methods).
'''

# =====================================================================================================================
#                                                   HELPER FUNCTIONS
# =====================================================================================================================

def CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome):
    _ = os.system('mkdir -p ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/SpliceAI')
    _ = os.system('echo "##fileformat=VCFv4.2" > ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/SpliceAI/input.vcf')
    _ = os.system(bcftools + ' view -h ' + workdir + '/' + gene_panel.split('-')[0] + '/' + patient_id + '/RNA/stripe/target_genes/longcallR/output.vcf.gz | grep "##contig=" >> ' + 
        workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/SpliceAI/input.vcf')
    item = variant_id.replace(' ', ':').replace('>', ':').split(':')
    pd.DataFrame([[item[0], item[1], '.', item[2], item[3], '.', '.', '.']], columns = ['#CHROM', 'POS', 'ID',
        'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']).to_csv(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/SpliceAI/input.vcf', 
        mode = 'a', sep = '\t', index = False)
    _ = os.system('(spliceai -I ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/SpliceAI/input.vcf -O ' +
        workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/SpliceAI/output.vcf -R ' + genome + ' -A grch38) > /dev/null 2>&1')
    if os.path.isfile(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/SpliceAI/output.vcf'):
        spliceai_results = pd.read_csv(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/SpliceAI/output.vcf', header = None, sep = '\t',
            names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'], comment = '#')
        spliceai_score = float(list(spliceai_results['INFO'].apply(lambda x: np.max([float(val) if val != '.' else 0 for val in x.split('|')[2:6]]) if 'SpliceAI' in x else np.nan))[0])
        _ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/SpliceAI/input.vcf && rm ' + workdir + 
            '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/SpliceAI/output.vcf && rmdir ' + workdir + 
            '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/SpliceAI')
    else:
        spliceai_score = 'NA'
        _ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/SpliceAI/input.vcf && rmdir ' + workdir + 
            '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/SpliceAI')
    return spliceai_score

def GetTotalCoverage(junctions):
    coverageDF = junctions[['Chrom', 'Start', 'End']].copy()
    lssDF = junctions.drop('End', axis = 1).groupby(['Chrom', 'Start']).sum().reset_index()
    rssDF = junctions.drop('Start', axis = 1).groupby(['Chrom', 'End']).sum().reset_index()
    samples = [sampid for sampid in junctions.columns if sampid not in {'Chrom', 'Start', 'End'}]
    for clm in samples:
        appendDF = pd.DataFrame(((coverageDF['Chrom'] + '_' + coverageDF['Start'].astype(str)).map(dict(zip(
            lssDF['Chrom'] + '_' + lssDF['Start'].astype(str), lssDF[clm]))).fillna(0) +
            (coverageDF['Chrom'] + '_' + coverageDF['End'].astype(str)).map(dict(zip(
            rssDF['Chrom'] + '_' + rssDF['End'].astype(str), rssDF[clm]))).fillna(0) - junctions[clm]), columns = [clm])
        coverageDF = pd.concat([coverageDF, appendDF], axis = 1)
    return coverageDF

def CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel):
    targetDF = pd.read_csv(workdir + '/' + gene_panel.split('-')[0] + '/references/target_genes.bed', sep = '\t', header = None)
    gene_id = targetDF[targetDF[4] == gene_name].values[0][3]
    _ = os.system('grep ' + gene_id + ' ' + gencode + ' > ' + workdir + '/tmp/input.gtf')
    geneDF = pd.read_csv(workdir + '/tmp/input.gtf', sep = '\t', header = None)
    geneDF = geneDF[geneDF[2] == 'exon']
    geneDF['transcript_id'] = geneDF[8].apply(lambda x: next(iter([item.split('"')[1] for item in x.split(';') if 'transcript_id' in item]), 'NA'))
    geneDF = geneDF.groupby(['transcript_id', 0, 6])[[3, 4]].agg(lambda x: sorted(list(x))).reset_index()
    geneDF['splice_sites'] = geneDF.apply(lambda x: [site for idx in range(len(x[3])-1) for site in [x[4][idx]+1, x[3][idx+1]-1]], axis = 1)
    sites = set(geneDF.explode('splice_sites').dropna()['splice_sites'])
    junctions = dict()
    for read in pysam.AlignmentFile(bamfile, 'rb').fetch():
        for item in ParseCigarJunctions(read.cigarstring, read.reference_name, read.reference_start):
            if any([int(ss) in sites for ss in item.split('_')[1:]]):
                junctions[item] = junctions.get(item, 0) + 1
    junctions = pd.DataFrame(junctions.items(), columns = ['Junction', 'Read_Counts'])
    junctions[['Chrom', 'Start', 'End']] = junctions['Junction'].str.split('_', n = 2, expand = True)
    junctions = junctions[['Chrom', 'Start', 'End', 'Read_Counts']]
    junctions['Coverage'] = GetTotalCoverage(junctions)['Read_Counts']
    junctions = junctions[(junctions['Read_Counts'] >= 5) & (junctions['Coverage'] >= 20) & 
        (junctions['Read_Counts']/junctions['Coverage'] >= 0.1)]
    junctions['Usage'] = junctions['Read_Counts']/junctions['Coverage']
    gtexDF = pd.read_csv(workdir + '/' + gene_panel.split('-')[0] + '/references/GTEx_v8/Fibroblast/target_genes.junction_counts.txt', 
        sep = '\t', header = 0)
    gtexDF[['Chrom', 'Start', 'End']] = gtexDF['Name'].str.split('_', n = 2, expand = True)
    gtexDF = gtexDF[['Chrom', 'Start', 'End'] + list(gtexDF.columns[1:-3])]
    ctrlCounts = gtexDF[(gtexDF['Chrom'].isin(set(junctions['Chrom']))) & (gtexDF['Start'].isin(set(junctions['Start'])) | 
        gtexDF['End'].isin(set(junctions['End'])))].drop_duplicates().reset_index(drop = True)
    ctrlCounts = junctions[['Chrom', 'Start', 'End']].merge(ctrlCounts, on = ['Chrom', 'Start', 'End'], 
        how = 'outer').sort_values(by = ['Chrom', 'Start', 'End']).fillna(0)
    ctrlCoverage = GetTotalCoverage(ctrlCounts)
    ctrlCounts = junctions[['Chrom', 'Start', 'End']].merge(ctrlCounts, on = ['Chrom', 'Start', 'End'], how = 'left')
    ctrlCoverage = junctions[['Chrom', 'Start', 'End']].merge(ctrlCoverage, on = ['Chrom', 'Start', 'End'], how = 'left')
    ctrlPSI = pd.concat([junctions[['Chrom', 'Start', 'End']].reset_index(drop = True), (ctrlCounts.iloc[:, 3:]/
        ctrlCoverage.iloc[:, 3:].mask(ctrlCoverage.iloc[:, 3:] < 20)).reset_index(drop = True)], axis = 1)
    ctrlPSI.iloc[:, 3:] = pd.DataFrame(ctrlPSI.iloc[:, 3:].apply(lambda x: MaskOutliers(np.array(x, dtype = float)), axis = 1).tolist())
    na_filter = list((~ctrlPSI.iloc[:, 3:].isna()).sum(axis = 1) >= 100)
    ctrlPSI, junctions = ctrlPSI[na_filter], junctions[na_filter]
    ctrlPSI.iloc[:, 3:] = np.round(RescaleValues(ctrlPSI.iloc[:, 3:], 0, 1, 1e-3, 1-1e-3), 3)
    parameters = list(ctrlPSI.iloc[:, 3:].apply(lambda x: FitBetaDist(np.array(x, dtype = float), 1e-3), axis = 1))
    junctions['GTEx_Usage'] = [param[0]/sum(param) for param in parameters]
    junctions['P_Value'] = ['{:,.3e}'.format(BetaBinomialTest(list(junctions['Read_Counts'])[idx], list(junctions['Coverage'])[idx],
        parameters[idx])) for idx in range(len(parameters))]
    _ = os.system('rm ' + workdir + '/tmp/input.gtf')
    return junctions

# =====================================================================================================================
#                                                         MAIN
# =====================================================================================================================

import sys

# Establish working directories and file paths
workdir = '/mnt/isilon/lin_lab_share/STRIPE'
bcftools = '/scr1/users/wangr5/tools/bcftools-1.21/bcftools'
genome = '/scr1/users/wangr5/references/GRCh38.primary_assembly.genome.fa'
gencode = '/scr1/users/wangr5/references/gencode.v45.annotation.gtf'
outfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S7/Table_S7.txt'

# Import functions used in Run_STRIPE.py
sys.path.append(workdir + '/scripts/main')
from Run_STRIPE import *
outDF = pd.DataFrame(columns = ['Individual', 'Gene', 'Variant', 'Variant_Class', 'SpliceAI', 'Junction_Coord',
    'Junction_Distance', 'Usage', 'GTEx_Usage', 'P_Value'])

# =====================================================================================================================
#                                                       BS2-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'BS2-1', 'PMD-359', 'EARS2'
variant_id, variant_class = 'chr16:Δ23525826-23551747', 'Structural deletion'

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 23525826) - min(int(x['End']), 23551747)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], 'NA'
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                       Q1695
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1695', 'PMD-359', 'ELAC2'
variant_id, variant_class = 'chr17:13016933 CT>A', 'Splice acceptor'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap2_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 13016933) - min(int(x['End']), 13016933)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                       Q1513
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1513', 'PMD-359', 'FBXL4'
variant_id, variant_class = 'chr6:98874306 A>T', 'Missense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 98874306) - min(int(x['End']), 98874306)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

variant_id, variant_class = 'chr6:98880546 CACTT>C', 'Splice donor region'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap2_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 98880546) - min(int(x['End']), 98880546)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                      FH-2209
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'FH-2209', 'PMD-359', 'FH'
variant_id, variant_class = 'chr1:241519682 A>AG', 'Frameshift duplication'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap2_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 241519682) - min(int(x['End']), 241519682)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                        JaWe
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'JaWe', 'CDG-466', 'MPDU1'
variant_id, variant_class = 'chr17:7583881 G>T', 'Nonsense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 7583881) - min(int(x['End']), 7583881)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                    FCDGC-02003
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'FCDGC-02003', 'CDG-466', 'NGLY1'
variant_id, variant_class = 'chr3:25729207 TTTGA>T', 'Frameshift deletion'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 25729207) - min(int(x['End']), 25729207)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                        Q1687
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1687', 'PMD-359', 'NUBPL'
variant_id, variant_class = 'chr14:31826715 G>A', 'Splice donor'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap2_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 31826715) - min(int(x['End']), 31826715)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                        AnJa
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'AnJa', 'CDG-466', 'PIGN'
variant_id, variant_class = 'chr18:62106985 C>G', 'Splice donor'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 62106985) - min(int(x['End']), 62106985)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

variant_id, variant_class = 'chr18:62113129 C>T', 'Splice donor region'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap2_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 62113129) - min(int(x['End']), 62113129)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                     CDG-110-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-110-1', 'CDG-466', 'PIGN'
variant_id, variant_class = 'chr18:62095902 C>T', 'Missense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap2_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 62095902) - min(int(x['End']), 62095902)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                     CDG-137-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-137-1', 'CDG-466', 'PIGN'
variant_id, variant_class = 'chr18:62106985 C>G', 'Splice donor'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 62106985) - min(int(x['End']), 62106985)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

variant_id, variant_class = 'chr18:62143337 A>C', 'Missense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap2_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 62143337) - min(int(x['End']), 62143337)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                     CDG-147-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-147-1', 'CDG-466', 'PIGN'
variant_id, variant_class = 'chr18:62113277 C>T', 'Missense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 62113277) - min(int(x['End']), 62113277)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

variant_id, variant_class = 'chr18:62143306 C>T', 'Splice donor region'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap2_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 62143306) - min(int(x['End']), 62143306)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                     CDG-161-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-161-1', 'CDG-466', 'PIGN'
variant_id, variant_class = 'chr18:62113129 C>T', 'Splice donor region'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 62113129) - min(int(x['End']), 62113129)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                     CDG-183-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-183-1', 'CDG-466', 'PIGN'
variant_id, variant_class = 'chr18:62106797 G>A', 'Nonsense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap2_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 62106797) - min(int(x['End']), 62106797)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                     CDG-184-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-184-1', 'CDG-466', 'PIGN'
variant_id, variant_class = 'chr18:62045973 G>C', 'Missense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 62045973) - min(int(x['End']), 62045973)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

variant_id, variant_class = 'chr18:62143337 A>C', 'Missense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap2_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 62143337) - min(int(x['End']), 62143337)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                     CDG-132-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-132-1', 'CDG-466', 'PIGQ'
variant_id, variant_class = 'chr16:576255 G>A', 'Splice donor'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 576255) - min(int(x['End']), 576255)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

variant_id, variant_class = 'chr16:578911 TCTA>T', 'Inframe deletion'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap2_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 578911) - min(int(x['End']), 578911)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                       11547
# =====================================================================================================================

patient_id, gene_panel, gene_name = '11547', 'CDG-466', 'PMM2'
variant_id, variant_class = 'chr16:8811146 G>A', 'Missense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 8811146) - min(int(x['End']), 8811146)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

variant_id, variant_class = 'chr16:8811153 G>A', 'Missense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap2_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 8811153) - min(int(x['End']), 8811153)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                     CDG-125-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-125-1', 'CDG-466', 'PMM2'
variant_id, variant_class = 'chr16:8811088 C>A', 'Missense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 8811088) - min(int(x['End']), 8811088)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

variant_id, variant_class = 'chr16:8811153 G>A', 'Missense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap2_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 8811153) - min(int(x['End']), 8811153)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                     CDG-152-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-152-1', 'CDG-466', 'PMM2'
variant_id, variant_class = 'chr16:Δ8805603-8820989', 'Structural deletion'

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap2_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 8805603) - min(int(x['End']), 8820989)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], 'NA'
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

variant_id, variant_class = 'chr16:8806398 C>T', 'Missense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 8806398) - min(int(x['End']), 8806398)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                     CDG-115-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-115-1', 'CDG-466', 'SLC35A2'
variant_id, variant_class = 'chrX:48909877 C>T', 'Missense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 48909877) - min(int(x['End']), 48909877)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                       08013
# =====================================================================================================================

patient_id, gene_panel, gene_name = '08013', 'CDG-466', 'SLC35C1'
variant_id, variant_class = 'chr11:45806301 CCTT>C', 'Inframe deletion'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 45806301) - min(int(x['End']), 45806301)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                       Q1819
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1819', 'PMD-359', 'TRIT1'
variant_id, variant_class = 'chr1:39883470 G>A', 'Nonsense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap2_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 39883470) - min(int(x['End']), 39883470)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# =====================================================================================================================
#                                                       Q1642
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1642', 'PMD-359', 'TRMU'
variant_id, variant_class = 'chr22:46356055 G>A', 'Missense'

# Compute SpliceAI score for variant
spliceai_score = CheckSpliceAI(workdir, patient_id, gene_name, gene_panel, bcftools, variant_id, genome)

# Retrieve junctions linked to haplotype of interest and characterize these junctions
bamfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/hap1_reads.bam'
junctions = CheckJunctions(workdir, bamfile, gencode, gene_name, gene_panel)
junctions['Junction_Distance'] = junctions.apply(lambda x: max(0, max(int(x['Start']), 46356055) - min(int(x['End']), 46356055)), axis = 1)
junctions['Individual'], junctions['Gene'], junctions['Variant'], junctions['Variant_Class'] = patient_id, gene_name, variant_id, variant_class
junctions['Junction_Coord'], junctions['SpliceAI'] = junctions['Chrom'] + ':' + junctions['Start'] + '-' + junctions['End'], spliceai_score
outDF = pd.concat([outDF, junctions[outDF.columns]]).reset_index(drop = True)

print('Done with ' + patient_id, flush = True)

# Save outDF to outfile
outDF.to_csv(outfile, sep = '\t', index = False)
