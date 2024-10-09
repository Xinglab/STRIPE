#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2024.10.04
Supplementary Table 11

Genes with unusually strong haplotype dosage imbalance detected and prioritized from TEQUILA-seq data for 15 undiagnosed patients.
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
outfile = workdir + '/manuscript/Supplementary_Tables/Table_S11/Table_S11.txt'

# Create an empty data frame with columns of interest
outDF = pd.DataFrame(columns = ['patient_id', 'gene_panel', 'gene_name', 'hap1_count', 'hap2_count', 'p_value'])

# Iterate over haplotype dosage imbalance results from STRIPE for each undiagnosed patient
targets = list(pd.read_csv(workdir + '/CDG/references/target_genes.bed', sep = '\t', header = None)[4])
sampleDF = pd.read_csv(workdir + '/CDG/samples.txt', sep = '\t', header = 0)
cohort = list(sampleDF[sampleDF['Provider'] != 'Lan Lin']['ID'])
undiagnosed = list(sampleDF[(sampleDF['Status'] == 'Undiagnosed') & (sampleDF['Provider'] != 'Lan Lin')]['ID'])

for patient in cohort:
    for target in targets:
        currdir = workdir + '/CDG/' + patient + '/RNA/stripe/target_genes/' + target
        if os.path.isfile(currdir + '/phasing_stats.txt'):
            # Only use phased reads if phasing quality score is greater than 0.75
            if pd.read_csv(currdir + '/phasing_stats.txt', sep = '\t', header = None).iloc[0,1] > 0.75:
                infile = currdir + '/haplotype_dosage/output.txt'
                if os.path.exists(infile) and os.path.getsize(infile) > 0:
                    inDF = pd.read_csv(infile, sep = '\t', header = None)
                    outDF = pd.concat([outDF, pd.DataFrame([[patient, 'CDG-466', target, inDF.iloc[0,1], inDF.iloc[1,1],
                        inDF.iloc[3,1]]], columns = outDF.columns)])

# Keep genes with FDR-adjusted p-value < 5%
_, fdr_values, _, _ = multipletests(outDF['p_value'], method = 'fdr_bh')
filterDF = outDF[fdr_values < 0.05].copy()

# Compute minor haplotype expression ratios for genes in filterDF among cohort samples
geneDF = filterDF[['gene_name']].drop_duplicates().reset_index(drop = True)

for patient in cohort:
    minor_ratios = []
    for target in list(geneDF['gene_name'].values):
        minor_ratio = np.nan
        currdir = workdir + '/CDG/' + patient + '/RNA/stripe/target_genes/' + target
        if os.path.isfile(currdir + '/phasing_stats.txt'):
            if pd.read_csv(currdir + '/phasing_stats.txt', sep = '\t', header = None).iloc[0,1] > 0.75:
                infile = currdir + '/haplotype_dosage/output.txt'
                if os.path.exists(infile) and os.path.getsize(infile) > 0:
                    inDF = pd.read_csv(infile, sep = '\t', header = None)
                    minor_ratio = min(inDF.iloc[0,1], inDF.iloc[1,1])/sum([inDF.iloc[0,1], inDF.iloc[1,1]])
        minor_ratios.append(minor_ratio)
    geneDF[patient] = minor_ratios
    
# Iterate over individuals in cohort and filter for genes meeting the following criteria:
#   * Minor haplotype expression ratio is lowest in a given individual
#   * Departure from allelic balance is more than 1.5 times that of the next lowest individual

cohort_filter = []
filterDF['minor_ratio'] = filterDF[['hap1_count', 'hap2_count']].min(axis = 1)/filterDF[['hap1_count', 'hap2_count']].sum(axis = 1)
for patient in cohort:
    currDF = filterDF[filterDF['patient_id'] == patient]
    if currDF.shape[0] > 0:
        currGene = dict(zip(geneDF['gene_name'], geneDF.drop(['gene_name', patient], axis = 1).min(axis = 1)))
        currFilter = (((filterDF[filterDF['patient_id'] == patient]['minor_ratio'] < filterDF[filterDF['patient_id'] == patient]['gene_name'].map(currGene))) &
            ((0.5 - filterDF[filterDF['patient_id'] == patient]['minor_ratio'])/(0.5 - filterDF[filterDF['patient_id'] == patient]['gene_name'].map(currGene)) > 1.5))
        cohort_filter += list(currFilter)

filterDF = filterDF[cohort_filter]

# Filter for results in undiagnosed patients
filterDF = filterDF[filterDF['patient_id'].isin(set(undiagnosed))]

# Save filterDF to outfile
filterDF = filterDF.drop(['minor_ratio'], axis = 1)
filterDF.to_csv(outfile, sep = '\t', index = False)
