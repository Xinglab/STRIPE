#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2025.01.13
Supplementary Table 10

Target genes with unusually strong haplotype dosage imbalance that were identified and prioritized from TEQUILA-seq data for cohort individuals.
We filtered target genes based on the following sequential criteria: 
    (i) have a minimum of 20 reads assigned to a haplotype in a given cohort individual and in at least 100 GTEx controls
    (ii) called a haplotype dosage outlier (FDR < 1%, > 10% shift from balanced expression of both gene haplotypes) in no more than 3 cohort individuals
We also report a "phasing quality score" to represent the proportion of reads that could be confidently assigned to a haplotype.
P values were computed using a gene expression dosage outlier test, which assesses whether haplotype expression ratios in a sample are consistent 
with those observed in tissue-matched GTEx controls (Methods).
'''

# =====================================================================================================================
#                                                       PACKAGES 
# =====================================================================================================================

# Load required packages
import sys
from statsmodels.stats.multitest import fdrcorrection

# =====================================================================================================================
#                                                         MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir = '/mnt/isilon/lin_lab_share/STRIPE'
outfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S10/Table_S10.txt'

# Import functions used in Run_STRIPE.py
sys.path.append(workdir + '/scripts/main')
from Run_STRIPE import *
output = []

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

# Iterate over all CDG patients
targets = pd.read_csv(workdir + '/CDG/references/target_genes.bed', sep = '\t', header = None)
sampleDF = pd.read_csv(workdir + '/CDG/samples.txt', sep = '\t', header = 0)
sampleDF = sampleDF[(sampleDF['Provider'] != 'Lan Lin') & (sampleDF['ID'] != 'E1877')]
cdg_cohort = metaDF[metaDF['ID'].isin(set(sampleDF['ID']))]
gtexDF = pd.read_csv(workdir + '/CDG/references/GTEx_v8/Fibroblast/target_genes.haplotype_expression.txt', sep = '\t', header = 0)
gtexDF.loc[:, 'name'] = gtexDF['name'].str.split('.').str[0]

for patient in cdg_cohort.values:
    for target in targets.values:
        # Check if patient already has phasing results for target gene in folder for Supplementary Table 6
        if os.path.isfile(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient[0] + '/' + target[4] + '/phasing_stats.txt'):
            infile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient[0] + '/' + target[4] + '/phasing_stats.txt'
        elif os.path.isfile(workdir + '/CDG/' + patient[0] + '/RNA/stripe/target_genes/' + target[4] + '/phasing_stats.txt'):
            infile = workdir + '/CDG/' + patient[0] + '/RNA/stripe/target_genes/' + target[4] + '/phasing_stats.txt'
        else:
            infile = None
        if infile is not None:
            phasingInfo = list(pd.read_csv(infile, sep = '\t', header = None).iloc[:3, 1])
            if phasingInfo[1] + phasingInfo[2] >= 20:
                ctrlDF = gtexDF[gtexDF['name'] == target[3]].iloc[:, 1:]
                if ctrlDF.shape[0] > 0:
                    ctrlCounts = [tuple(map(int, item.split('|'))) for item in ctrlDF.values[0]]
                    ctrlRatios = [min(item)/sum(item) if sum(item) >= 20 else np.nan for item in ctrlCounts]
                    ctrlRatios = MaskOutliers(np.array(ctrlRatios)) # mask any haplotype expression ratios that are "boxplot" outliers
                    if np.sum(~np.isnan(ctrlRatios)) >= 100:
                        ctrlRatios = np.round(RescaleValues(ctrlRatios, 0, 1, 1e-3, 1-1e-3), 3)
                        myFunc = lambda x: 2*(digamma(x)-digamma(2*x)) - np.nanmean(np.log(ctrlRatios*(1-ctrlRatios)))
                        if np.sign(myFunc(1)) != np.sign(myFunc(999)):
                            fit = brentq(myFunc, 1, 999, full_output = True, disp = False)
                            if fit[1].converged:
                                pval = BetaBinomialTest(min(phasingInfo[1], phasingInfo[2]), phasingInfo[1] + phasingInfo[2], [fit[0]] * 2)
                                output += [[patient[0], patient[1], patient[2], patient[3], 'CDG-466', target[4], phasingInfo[0],
                                    phasingInfo[1], phasingInfo[2], 1/(8*fit[0]+4), pval]]

# Iterate over all PMD patients
targets = pd.read_csv(workdir + '/PMD/references/target_genes.bed', sep = '\t', header = None)
sampleDF = pd.read_csv(workdir + '/PMD/samples.txt', sep = '\t', header = 0)
sampleDF = sampleDF[(sampleDF['Provider'] == 'Rebecca Ganetzky') | ((sampleDF['Provider'] == 'Marni Falk') & 
    (sampleDF['Status'] != 'Diagnosed') & ~sampleDF['ID'].str.contains('-MF'))]
pmd_cohort = metaDF[metaDF['ID'].isin(set(sampleDF['ID']))]
gtexDF = pd.read_csv(workdir + '/PMD/references/GTEx_v8/Fibroblast/target_genes.haplotype_expression.txt', sep = '\t', header = 0)
gtexDF.loc[:, 'name'] = gtexDF['name'].str.split('.').str[0]

for patient in pmd_cohort.values:
    for target in targets.values:
        # Check if patient already has phasing results for target gene in folder for Supplementary Table 6
        if os.path.isfile(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient[0] + '/' + target[4] + '/phasing_stats.txt'):
            infile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient[0] + '/' + target[4] + '/phasing_stats.txt'
        elif os.path.isfile(workdir + '/PMD/' + patient[0] + '/RNA/stripe/target_genes/' + target[4] + '/phasing_stats.txt'):
            infile = workdir + '/PMD/' + patient[0] + '/RNA/stripe/target_genes/' + target[4] + '/phasing_stats.txt'
        else:
            infile = None
        if infile is not None:
            phasingInfo = list(pd.read_csv(infile, sep = '\t', header = None).iloc[:3, 1])
            if phasingInfo[1] + phasingInfo[2] >= 20:
                ctrlDF = gtexDF[gtexDF['name'] == target[3]].iloc[:, 1:]
                if ctrlDF.shape[0] > 0:
                    ctrlCounts = [tuple(map(int, item.split('|'))) for item in ctrlDF.values[0]]
                    ctrlRatios = [min(item)/sum(item) if sum(item) >= 20 else np.nan for item in ctrlCounts]
                    ctrlRatios = MaskOutliers(np.array(ctrlRatios)) # mask any haplotype expression ratios that are "boxplot" outliers
                    if np.sum(~np.isnan(ctrlRatios)) >= 100:
                        ctrlRatios = np.round(RescaleValues(ctrlRatios, 0, 1, 1e-3, 1-1e-3), 3)
                        myFunc = lambda x: 2*(digamma(x)-digamma(2*x)) - np.nanmean(np.log(ctrlRatios*(1-ctrlRatios)))
                        if np.sign(myFunc(1)) != np.sign(myFunc(999)):
                            fit = brentq(myFunc, 1, 999, full_output = True, disp = False)
                            if fit[1].converged:
                                pval = BetaBinomialTest(min(phasingInfo[1], phasingInfo[2]), phasingInfo[1] + phasingInfo[2], [fit[0]] * 2)
                                output += [[patient[0], patient[1], patient[2], patient[3], 'PMD-359', target[4], phasingInfo[0],
                                    phasingInfo[1], phasingInfo[2], 1/(8*fit[0]+4), pval]]

# Convert output into a dataframe
outDF = pd.DataFrame(output, columns = ['id', 'status', 'disease_group', 'category', 'gene_panel', 'gene_name', 'phasing_quality', 'hap1_count', 
    'hap2_count', 'gtex_dosage_variation', 'p_value'])

# Filter out target genes that were called haplotype dosage outliers (FDR < 1%, > 10% shift from balanced expression of both gene haplotypes) in more than 3 cohort individuals
outDF['fdr'] = fdrcorrection(outDF['p_value'], alpha = 0.01, is_sorted = False)[1]
filterDF = outDF[outDF['fdr'] < 0.01]
filterDF = filterDF[abs(0.5 - (filterDF['hap1_count']/(filterDF['hap1_count'] + filterDF['hap2_count']))) > 0.1]
gene_counts = filterDF['gene_name'].value_counts()
filterDF = filterDF[filterDF['gene_name'].map(dict(zip(gene_counts.index, gene_counts.values))) <= 3]

# Save filterDF to outfile
filterDF.drop('fdr', axis = 1).to_csv(outfile, sep = '\t', index = False)

# Calculate median number of target genes prioritized among 44 cohort individuals for whom prior genetic testing revealed only a single 
# pathogenic variant in a gene associated with a recessive condition (n = 6) or no candidate disease-causing variants (n = 38)
counts = filterDF['id'].value_counts()
np.median([dict(zip(counts.index, counts.values)).get(sample_id, 0) for sample_id in metaDF[metaDF['Category'].isin({'Single variant, recessive condition',
    'No candidate variants'})]['ID']]) # 1 target gene per patient
