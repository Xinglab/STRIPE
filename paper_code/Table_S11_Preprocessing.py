#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2025.01.13
Supplementary Table 11 (Data Pre-Processing)
'''

# =====================================================================================================================
#                                                       PACKAGES 
# =====================================================================================================================

# Load required packages
import sys

# =====================================================================================================================
#                                                         MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir = '/mnt/isilon/lin_lab_share/STRIPE'
samtools = '/scr1/users/wangr5/tools/samtools-1.21/samtools'
gencode = '/scr1/users/wangr5/references/gencode.v45.annotation.gtf'
outdir = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S11'

# Parse command line arguments
panel, gene = sys.argv[1], sys.argv[2]

# Import functions used in Run_STRIPE.py
sys.path.append(workdir + '/scripts/main')
from Run_STRIPE import *

# Create folder for input panel/gene
_ = os.system('mkdir -p ' + outdir + '/' + panel)
_ = os.system('mkdir -p ' + outdir + '/' + panel + '/' + gene)
_ = os.system('mkdir -p ' + outdir + '/' + panel + '/' + gene + '/tmp')
currdir = outdir + '/' + panel + '/' + gene
output = []

def ExtractSpliceJunctions(infile, samtools, tmpfile, sites, targets = None, region = None):
    '''
    Function designed to extract splice junctions from primary alignments (MAPQ >= 1) in infile
    '''
    if targets is None and region is None:
        _ = os.system(samtools + ' view -hb -F 256 -q 1 ' + infile + ' > ' + tmpfile + 
            ' && ' + samtools + ' index ' + tmpfile)
    elif targets is None:
        _ = os.system(samtools + ' view -hb -F 256 -q 1 ' + infile + ' ' + region + ' > ' + tmpfile +
            ' && ' + samtools + ' index ' + tmpfile)
    elif region is None:
        _ = os.system(samtools + ' view -hb -L ' + targets + ' -F 256 -q 1 ' + infile + 
            ' > ' + tmpfile + ' && ' + samtools + ' index ' + tmpfile)
    junctions = dict()
    for read in pysam.AlignmentFile(tmpfile, 'rb').fetch():
        for item in ParseCigarJunctions(read.cigarstring, read.reference_name, read.reference_start):
            if any([int(ss) in sites for ss in item.split('_')[1:]]):
                junctions[item] = junctions.get(item, 0) + 1
    _ = os.system('rm ' + tmpfile + ' && rm ' + tmpfile + '.bai')
    junctions = pd.DataFrame(junctions.items(), columns = ['Junction', 'Read_Counts'])
    if junctions.shape[0] > 0:
        junctions[['Chrom', 'Start', 'End']] = junctions['Junction'].str.split('_', n = 2, expand = True)
        return junctions.sort_values(by = ['Chrom', 'Start', 'End'], ascending = True).reset_index(drop = True)[['Chrom', 'Start', 'End', 'Read_Counts']]
    else:
        return pd.DataFrame(columns = ['Chrom', 'Start', 'End', 'Read_Counts'])

# Establish information on cohort individuals
metaDF = pd.DataFrame([['18-10792', 'Healthy', '', 'Control'], ['CFB-CHOP-1', 'Healthy', '', 'Control'], ['E1877', 'Healthy', '', 'Control'],
    ['E1878', 'Healthy', '', 'Control'], ['E4722', 'Healthy', '', 'Control'], ['E4723', 'Healthy', '', 'Control'],
    ['GM05399', 'Healthy', '', 'Control'], ['GM05565', 'Healthy', '', 'Control'], ['GM08398', 'Healthy', '', 'Control'],
    ['PCS-201-010', 'Healthy', '', 'Control'], ['Q1007p1', 'Healthy', '', 'Control'], ['Q1007p2', 'Healthy', '', 'Control'],
    ['Q1269p1', 'Healthy', '', 'Control'], ['Q1508p1', 'Healthy', '', 'Control'], ['Q1508p2', 'Healthy', '', 'Control'],
    ['Q1508g2', 'Healthy', '', 'Control'], ['Q1877p1', 'Healthy', '', 'Control'], ['Q1877p2', 'Healthy', '', 'Control'],
    ['Q1979p1', 'Healthy', '', 'Control'], ['Q1979p2', 'Healthy', '', 'Control'],
    ['CDG-114-1', 'Affected', 'Congenital disorders of glycosylation', 'Known genetic cause'],
    ['AnJa', 'Affected', 'Congenital disorders of glycosylation', 'Known genetic cause'],
    ['CDG-137-1', 'Affected', 'Congenital disorders of glycosylation', 'Known genetic cause'],
    ['CDG-147-1', 'Affected', 'Congenital disorders of glycosylation', 'Known genetic cause'],
    ['CDG-184-1', 'Affected', 'Congenital disorders of glycosylation', 'Known genetic cause'],
    ['CDG-132-1', 'Affected', 'Congenital disorders of glycosylation', 'Known genetic cause'],
    ['11547', 'Affected', 'Congenital disorders of glycosylation', 'Known genetic cause'],
    ['CDG-125-1', 'Affected', 'Congenital disorders of glycosylation', 'Known genetic cause'],
    ['CDG-152-1', 'Affected', 'Congenital disorders of glycosylation', 'Known genetic cause'],
    ['CDG-115-1', 'Affected', 'Congenital disorders of glycosylation', 'Known genetic cause'],
    ['ATP6AP2', 'Affected', 'Congenital disorders of glycosylation', 'Variant(s) of uncertain significance'],
    ['CDG-110-1', 'Affected', 'Congenital disorders of glycosylation', 'Variant(s) of uncertain significance'],
    ['CDG-161-1', 'Affected', 'Congenital disorders of glycosylation', 'Variant(s) of uncertain significance'],
    ['CDG-183-1', 'Affected', 'Congenital disorders of glycosylation', 'Variant(s) of uncertain significance'],
    ['JaWe', 'Affected', 'Congenital disorders of glycosylation', 'Single variant, recessive condition'],
    ['FCDGC-02003', 'Affected', 'Congenital disorders of glycosylation', 'Single variant, recessive condition'],
    ['08013', 'Affected', 'Congenital disorders of glycosylation', 'Single variant, recessive condition'],
    ['CDG-150-1', 'Affected', 'Congenital disorders of glycosylation', 'Single variant, recessive condition'],
    ['CDG-IIx-1', 'Affected', 'Congenital disorders of glycosylation', 'No candidate variants'],
    ['CDG-IIx-2', 'Affected', 'Congenital disorders of glycosylation', 'No candidate variants'],
    ['CDG-0346', 'Affected', 'Congenital disorders of glycosylation', 'No candidate variants'],
    ['CDG-0367', 'Affected', 'Congenital disorders of glycosylation', 'No candidate variants'],
    ['CDG-0414', 'Affected', 'Congenital disorders of glycosylation', 'No candidate variants'],
    ['CDG-167-4', 'Affected', 'Congenital disorders of glycosylation', 'No candidate variants'],
    ['CDG-170-1', 'Affected', 'Congenital disorders of glycosylation', 'No candidate variants'],
    ['CDG-176-1', 'Affected', 'Congenital disorders of glycosylation', 'No candidate variants'],
    ['CDG-182-1', 'Affected', 'Congenital disorders of glycosylation', 'No candidate variants'],
    ['CDG-186-1', 'Affected', 'Congenital disorders of glycosylation', 'No candidate variants'],
    ['SuHu', 'Affected', 'Congenital disorders of glycosylation', 'No candidate variants'],
    ['Q1663', 'Affected', 'Primary mitochondrial diseases', 'Known genetic cause'],
    ['Q1964', 'Affected', 'Primary mitochondrial diseases', 'Known genetic cause'],
    ['Q1513', 'Affected', 'Primary mitochondrial diseases', 'Known genetic cause'],
    ['Q1363', 'Affected', 'Primary mitochondrial diseases', 'Known genetic cause'],
    ['Q1642', 'Affected', 'Primary mitochondrial diseases', 'Known genetic cause'],
    ['Q1695', 'Affected', 'Primary mitochondrial diseases', 'Variant(s) of uncertain significance'],
    ['FH-2209', 'Affected', 'Primary mitochondrial diseases', 'Variant(s) of uncertain significance'],
    ['Q2032s1', 'Affected', 'Primary mitochondrial diseases', 'Variant(s) of uncertain significance'],
    ['Q2319', 'Affected', 'Primary mitochondrial diseases', 'Variant(s) of uncertain significance'],
    ['Q1819', 'Affected', 'Primary mitochondrial diseases', 'Variant(s) of uncertain significance'],
    ['BS2-1', 'Affected', 'Primary mitochondrial diseases', 'Single variant, recessive condition'],
    ['Q1687', 'Affected', 'Primary mitochondrial diseases', 'Single variant, recessive condition'],
	['Q1043a1', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1161', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1279', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1339', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
    ['Q1343', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1347', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1349', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1357', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1423', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1449', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1507', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1544', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1564', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1574p1', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1577s1', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1616', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1617', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1648', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1724', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1807', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1816', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1832', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1842', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q1849', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q2002', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q2012', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants'],
	['Q2071', 'Affected', 'Primary mitochondrial diseases', 'No candidate variants']], 
    columns = ['ID', 'Status', 'Disease', 'Category'])

targets = pd.read_csv(workdir + '/' + panel.split('-')[0] + '/references/target_genes.bed', sep = '\t', header = None)
target = targets[targets[4] == gene].iloc[0]
sampleDF = pd.read_csv(workdir + '/' + panel.split('-')[0] + '/samples.txt', sep = '\t', header = 0)
if panel.split('-')[0] == 'CDG':
    sampleDF = sampleDF[(sampleDF['Provider'] != 'Lan Lin') & (sampleDF['ID'] != 'E1877')]
else:
    sampleDF = sampleDF[(sampleDF['Provider'] == 'Rebecca Ganetzky') | ((sampleDF['Provider'] == 'Marni Falk') & 
        (sampleDF['Status'] != 'Diagnosed') & ~sampleDF['ID'].str.contains('-MF'))]
cohort = metaDF[metaDF['ID'].isin(set(sampleDF['ID']))]
gtexDF = pd.read_csv(workdir + '/' + panel.split('-')[0] + '/references/GTEx_v8/Fibroblast/target_genes.junction_counts.txt', sep = '\t', header = 0)
gtexDF[['Chrom', 'Start', 'End']] = gtexDF['Name'].str.split('_', n = 2, expand = True)
gtexDF = gtexDF[['Chrom', 'Start', 'End'] + list(gtexDF.columns[1:-3])]

_ = os.system('grep ' + target[3] + ' ' + gencode + ' > ' + currdir + '/tmp/input.gtf')
geneDF = pd.read_csv(currdir + '/tmp/input.gtf', sep = '\t', header = None)
geneDF = geneDF[geneDF[2] == 'exon']
geneDF['transcript_id'] = geneDF[8].apply(lambda x: next(iter([item.split('"')[1] for item in x.split(';') if 'transcript_id' in item]), 'NA'))
geneDF = geneDF.groupby(['transcript_id', 0, 6])[[3, 4]].agg(lambda x: sorted(list(x))).reset_index()
geneDF['splice_sites'] = geneDF.apply(lambda x: [site for idx in range(len(x[3])-1) for site in [x[4][idx]+1, x[3][idx+1]-1]], axis = 1)
sites = set(geneDF.explode('splice_sites').dropna()['splice_sites'])
for patient in cohort.values:
    infile = workdir + '/' + panel.split('-')[0] + '/' + patient[0] + '/RNA/' + patient[0] + '_TEQUILA.bam'
    allJunctions = ExtractSpliceJunctions(infile, samtools, currdir + '/tmp/all_reads.bam', sites, region = target[0] + ':' + str(target[1]+1) + '-' + 
        str(target[2]))
    allJunctions.columns = ['Chrom', 'Start', 'End', 'All']
    if allJunctions.shape[0] > 0:
        if os.path.isfile(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient[0] + '/' + target[4] + '/phasing_stats.txt'):
            phasedir = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient[0] + '/' + target[4]
            phasingQuality = float(pd.read_csv(phasedir + '/phasing_stats.txt', sep = '\t', header = None).iloc[0, 1])
        elif os.path.isfile(workdir + '/' + panel.split('-')[0] + '/' + patient[0] + '/RNA/stripe/target_genes/' + target[4] + '/phasing_stats.txt'):
            phasedir = workdir + '/' + panel.split('-')[0] + '/' + patient[0] + '/RNA/stripe/target_genes/' + target[4]
            phasingQuality = float(pd.read_csv(phasedir + '/phasing_stats.txt', sep = '\t', header = None).iloc[0, 1])
        else:
            phasedir = None
            phasingQuality = ''
        if phasedir is not None:
            hapJunctions1 = ExtractSpliceJunctions(phasedir + '/hap1_reads.bam', samtools, currdir + '/tmp/hap1_reads.bam', sites, region = target[0] + 
                ':' + str(target[1]+1) + '-' + str(target[2]))
            hapJunctions1.columns = ['Chrom', 'Start', 'End', 'Hap1']
            hapJunctions2 = ExtractSpliceJunctions(phasedir + '/hap2_reads.bam', samtools, currdir + '/tmp/hap2_reads.bam', sites, region = target[0] + 
                ':' + str(target[1]+1) + '-' + str(target[2]))
            hapJunctions2.columns = ['Chrom', 'Start', 'End', 'Hap2']
            allJunctions = allJunctions.merge(hapJunctions1, on = ['Chrom', 'Start', 'End'], how = 'outer').merge(hapJunctions2, 
                on = ['Chrom', 'Start', 'End'], how = 'outer').sort_values(by = ['Chrom', 'Start', 'End']).fillna(0)
        allCoverage = GetTotalCoverage(allJunctions)
        countFilter = list((allCoverage[allCoverage.columns[3:]].max(axis = 1) >= 20) & (allJunctions[allJunctions.columns[3:]].max(axis = 1) >= 20))
        allJunctions, allCoverage = allJunctions[countFilter], allCoverage[countFilter]
        if sum(countFilter) > 0:
            ctrlCounts = gtexDF[(gtexDF['Chrom'] == target[0]) & (gtexDF['Start'].isin(set(allJunctions['Start'])) | 
                gtexDF['End'].isin(set(allJunctions['End'])))].drop_duplicates().reset_index(drop = True)
            ctrlCounts = allJunctions[['Chrom', 'Start', 'End']].merge(ctrlCounts, on = ['Chrom', 'Start', 'End'], 
                how = 'outer').sort_values(by = ['Chrom', 'Start', 'End']).fillna(0)
            ctrlCoverage = GetTotalCoverage(ctrlCounts)
            ctrlCounts = allJunctions[['Chrom', 'Start', 'End']].merge(ctrlCounts, on = ['Chrom', 'Start', 'End'], how = 'left')
            ctrlCoverage = allJunctions[['Chrom', 'Start', 'End']].merge(ctrlCoverage, on = ['Chrom', 'Start', 'End'], how = 'left')
            ctrlPSI = pd.concat([allJunctions[['Chrom', 'Start', 'End']].reset_index(drop = True), 
                (ctrlCounts.iloc[:, 3:]/ctrlCoverage.iloc[:, 3:].mask(ctrlCoverage.iloc[:, 3:] < 20)).reset_index(drop = True)], axis = 1)
            ctrlPSI.iloc[:, 3:] = pd.DataFrame(ctrlPSI.iloc[:, 3:].apply(lambda x: MaskOutliers(np.array(x, dtype = float)), axis = 1).tolist())
            na_filter = list((~ctrlPSI.iloc[:, 3:].isna()).sum(axis = 1) >= 100)
            ctrlPSI, allJunctions, allCoverage = ctrlPSI[na_filter], allJunctions[na_filter], allCoverage[na_filter]
            if sum(na_filter) > 0:
                ctrlPSI.iloc[:, 3:] = np.round(RescaleValues(ctrlPSI.iloc[:, 3:], 0, 1, 1e-3, 1-1e-3), 3)
                parameters = list(ctrlPSI.iloc[:, 3:].apply(lambda x: FitBetaDist(np.array(x, dtype = float), 1e-3), axis = 1))
                allDF = pd.DataFrame({
                    'Junction': list(allJunctions['Chrom'] + ':' + allJunctions['Start'].astype(str) + '-' + allJunctions['End'].astype(str)),
                    'Count': list(allJunctions['All'].astype(int).astype(object)), 'Coverage': list(allCoverage['All'].astype(int).astype(object)),
                    'Population': ['{:,.3f}'.format(item[0]/(item[0] + item[1])) for item in parameters],
                    'PVal': ['{:,.3e}'.format(BetaBinomialTest(allJunctions.iloc[idx]['All'], allCoverage.iloc[idx]['All'],
                        parameters[idx])) for idx in range(len(parameters))]})
                output += [[patient[0], patient[1], patient[2], patient[3], panel, target[4], phasingQuality, 'All'] + 
                    list(item) for item in allDF.values]
                if 'Hap1' in allJunctions.columns and 'Hap2' in allJunctions.columns:
                    hap1DF = pd.DataFrame({
                        'Junction': list(allJunctions['Chrom'] + ':' + allJunctions['Start'].astype(str) + '-' + allJunctions['End'].astype(str)),
                        'Count': list(allJunctions['Hap1'].astype(int).astype(object)), 'Coverage': list(allCoverage['Hap1'].astype(int).astype(object)),
                        'Population': ['{:,.3f}'.format(item[0]/(item[0] + item[1])) for item in parameters],
                        'PVal': ['{:,.3e}'.format(BetaBinomialTest(allJunctions.iloc[idx]['Hap1'], allCoverage.iloc[idx]['Hap1'],
                            parameters[idx])) for idx in range(len(parameters))]})
                    output += [[patient[0], patient[1], patient[2], patient[3], panel, target[4], phasingQuality, 'Haplotype 1'] + 
                        list(item) for item in hap1DF.values]
                    hap2DF = pd.DataFrame({
                        'Junction': list(allJunctions['Chrom'] + ':' + allJunctions['Start'].astype(str) + '-' + allJunctions['End'].astype(str)),
                        'Count': list(allJunctions['Hap2'].astype(int).astype(object)), 'Coverage': list(allCoverage['Hap2'].astype(int).astype(object)),
                        'Population': ['{:,.3f}'.format(item[0]/(item[0] + item[1])) for item in parameters],
                        'PVal': ['{:,.3e}'.format(BetaBinomialTest(allJunctions.iloc[idx]['Hap2'], allCoverage.iloc[idx]['Hap2'],
                            parameters[idx])) for idx in range(len(parameters))]})
                    output += [[patient[0], patient[1], patient[2], patient[3], panel, target[4], phasingQuality, 'Haplotype 2'] + 
                        list(item) for item in hap2DF.values]

# Convert output into a dataframe and save to output folder
outDF = pd.DataFrame(output, columns = ['id', 'status', 'disease_group', 'category', 'gene_panel', 'gene_name', 'phasing_quality', 'read_group', 
    'junction_coord', 'junction_read_count', 'splice_site_coverage', 'gtex_junction_usage', 'p_value'])
outDF.to_csv(currdir + '/output.tsv', sep = '\t', index = False)
_ = os.system('rm ' + currdir + '/tmp/input.gtf' + ' && rmdir ' + currdir + '/tmp')
