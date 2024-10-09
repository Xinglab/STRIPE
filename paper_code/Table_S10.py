#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2024.10.04
Supplementary Table 10

Aberrant splice junctions detected and prioritized from TEQUILA-seq data for 15 undiagnosed patients.
'''

# =====================================================================================================================
#                                                       PACKAGES 
# =====================================================================================================================

# Load required packages
import os
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

# =====================================================================================================================
#                                                         MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir = '/mnt/isilon/lin_lab_share/STRIPE'
outfile = workdir + '/manuscript/Supplementary_Tables/Table_S10/Table_S10.txt'

# Create an empty data frame with columns of interest
outDF = pd.DataFrame(columns = ['patient_id', 'gene_panel', 'gene_name', 'read_group', 'splice_junction', 'junction_count', 
    'junction_coverage', 'junction_usage', 'gtex_usage', 'usage_shift', 'p_value'])

# Iterate over aberrant splice junction from STRIPE for each undiagnosed patient
targets = list(pd.read_csv(workdir + '/CDG/references/target_genes.bed', sep = '\t', header = None)[4])
sampleDF = pd.read_csv(workdir + '/CDG/samples.txt', sep = '\t', header = 0)
cohort = list(sampleDF[sampleDF['Provider'] != 'Lan Lin']['ID'])
undiagnosed = list(sampleDF[(sampleDF['Status'] == 'Undiagnosed') & (sampleDF['Provider'] != 'Lan Lin')]['ID'])

for patient in cohort:
    for target in targets:
        currdir, mode = workdir + '/CDG/' + patient + '/RNA/stripe/target_genes/' + target, 'all'
        if os.path.isfile(currdir + '/phasing_stats.txt'):
            # Only use phased reads if phasing quality score is greater than 0.75
            if pd.read_csv(currdir + '/phasing_stats.txt', sep = '\t', header = None).iloc[0,1] > 0.75:
                mode = 'phased'
        if mode == 'phased':
            for i in range(1,3):
                infile = currdir + '/aberrant_splicing/output_hap' + str(i) + '.txt'
                if os.path.exists(infile) and os.path.getsize(infile) > 0:
                    outDF = pd.concat([outDF, pd.DataFrame([[patient, 'CDG-466', target, 'Haplotype ' + str(i)] + list(item) for item in 
                        pd.read_csv(infile, sep = '\t', header = 0).values], columns = outDF.columns)])
        else:
            infile = currdir + '/aberrant_splicing/output_all.txt'
            if os.path.exists(infile) and os.path.getsize(infile) > 0:
                outDF = pd.concat([outDF, pd.DataFrame([[patient, 'CDG-466', target, 'All'] + list(item) for item in 
                    pd.read_csv(infile, sep = '\t', header = 0).values], columns = outDF.columns)])

# Keep junctions with FDR-adjusted p-value < 1%
_, fdr_values, _, _ = multipletests(outDF['p_value'], method = 'fdr_bh')
filterDF = outDF[fdr_values < 0.01].copy()

# Compute usage frequencies for junctions in filterDF among cohort samples
sjDF = filterDF[['splice_junction']].drop_duplicates().reset_index(drop = True)
sjDF[['Chrom', 'Start', 'End']] = sjDF['splice_junction'].str.split(r':|-', expand = True)

for patient in cohort:
    infile = workdir + '/CDG/' + patient + '/RNA/stripe/quality_control/' + patient + '_TEQUILA.junctions.tsv'
    currDF = pd.read_csv(infile, sep = '\t', header = 0)
    lssDF = currDF.drop(['End'], axis = 1).groupby(['Chrom', 'Start']).sum().reset_index()
    rssDF = currDF.drop(['Start'], axis = 1).groupby(['Chrom', 'End']).sum().reset_index()
    juncDict = dict(zip(currDF['Chrom'] + '_' + currDF['Start'].astype(str) + '_' + currDF['End'].astype(str), currDF['Read_Counts']))
    lssDict = dict(zip(lssDF['Chrom'] + '_' + lssDF['Start'].astype(str), lssDF['Read_Counts']))
    rssDict = dict(zip(rssDF['Chrom'] + '_' + rssDF['End'].astype(str), rssDF['Read_Counts']))
    coverage = ((sjDF['Chrom'] + '_' + sjDF['Start']).map(lssDict).fillna(0) + (sjDF['Chrom'] + '_' + sjDF['End']).map(rssDict).fillna(0) - 
        (sjDF['Chrom'] + '_' + sjDF['Start'] + '_' + sjDF['End']).map(juncDict).fillna(0))
    sjDF[patient] = (sjDF['Chrom'] + '_' + sjDF['Start'] + '_' + sjDF['End']).map(juncDict).fillna(0)/coverage.mask(coverage <= 20, np.nan)

# Iterate over individuals in cohort and filter for junctions meeting the following criteria:
#   * Usage frequency for splice junction is the highest in a given individual
#   * Usage frequency for splice junction is twice that of the next highest individual

cohort_filter = []
for patient in cohort:
    currDF = filterDF[filterDF['patient_id'] == patient]
    currPSI = dict(zip(sjDF['splice_junction'], sjDF.drop(['splice_junction', 'Chrom', 'Start', 'End', patient], axis = 1).max(axis = 1)))
    currFilter = ((filterDF[filterDF['patient_id'] == patient]['junction_usage'] > filterDF[filterDF['patient_id'] == patient]['splice_junction'].map(currPSI)) &
        (filterDF[filterDF['patient_id'] == patient]['junction_usage'] > 2*filterDF[filterDF['patient_id'] == patient]['splice_junction'].map(currPSI)))
    cohort_filter += list(currFilter)

filterDF = filterDF[cohort_filter]

# Keep junctions called as outliers in no more than 1 individual and have a usage frequency shift > 20% relative to GTEx controls
junction_counts = filterDF[['patient_id', 'splice_junction']].drop_duplicates()['splice_junction'].value_counts()
filterDF = filterDF[(filterDF['splice_junction'].isin(set(junction_counts[junction_counts == 1].index))) & (filterDF['usage_shift'] > 0.2)]

# Filter for outlier junctions in undiagnosed patients
filterDF = filterDF[filterDF['patient_id'].isin(set(undiagnosed))]

# Save filterDF to outfile
filterDF.to_csv(outfile, sep = '\t', index = False)
