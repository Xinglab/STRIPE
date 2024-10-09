#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2024.10.04
Supplementary Table 9

Candidate disease-causing variants detected and prioritized from TEQUILA-seq data for 15 undiagnosed patients.
'''

# =====================================================================================================================
#                                                       PACKAGES 
# =====================================================================================================================

# Load required packages
import os
import pandas as pd

# =====================================================================================================================
#                                                         MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir = '/mnt/isilon/lin_lab_share/STRIPE'
outfile = workdir + '/manuscript/Supplementary_Tables/Table_S9/Table_S9.txt'

# Create an empty data frame with columns of interest
outDF = pd.DataFrame(columns = ['patient_id', 'gene_panel', 'gene_name', 'variant_grch38'] + 
    [tool + '_' + field for tool in ['longcallR', 'clair3', 'deepvariant', 'phasing'] for field in ['dp', 'af', 'gt']] +
    ['cadd', 'spliceai', 'gnomad_af_grpmax_joint', 'gnomad_nhomalt_grpmax_joint', 'clinvar', 'homopolymer',
    'end_distance_pct', 'skipping_frequency'])

# Iterate over variant calling results from STRIPE for each undiagnosed patient
targets = list(pd.read_csv(workdir + '/CDG/references/target_genes.bed', sep = '\t', header = None)[4])
sampleDF = pd.read_csv(workdir + '/CDG/samples.txt', sep = '\t', header = 0)
cohort = list(sampleDF[sampleDF['Provider'] != 'Lan Lin']['ID'])
undiagnosed = list(sampleDF[(sampleDF['Status'] == 'Undiagnosed') & (sampleDF['Provider'] != 'Lan Lin')]['ID'])

for patient in cohort:
    for target in targets:
        infile = workdir + '/CDG/' + patient + '/RNA/stripe/target_genes/' + target + '/variant_calling/merged_variants.tsv'
        if os.path.isfile(infile):
            inDF = pd.read_csv(infile, sep = '\t', header = 0)
            # Retain variants that meet the following criteria:
            #   * No homozygous individuals in gnomAD
            #   * Not annotated as Benign or Likely benign in ClinVar
            #   * Does not fall within a homopolymer region
            inDF = inDF[(inDF['GNOMAD_NHOMALT_GRPMAX_JOINT'] == 0) & ~inDF['CLINVAR'].str.contains('enign') & ~inDF['HOMOPOLYMER']]
            if inDF.shape[0] > 0:
                outDF = pd.concat([outDF, pd.DataFrame([[patient, 'CDG-466', target, item[0] + ':' + str(item[1]) + ' ' + item[2] + '>' + item[3]] + list(item[4:]) 
                    for item in inDF.values], columns = outDF.columns)])

# Only keep variants that appear in no more than one individual and meet the following criteria:
#   * Proportion of reads with variant found at end < 0.5
#   * Proportion of reads skipping variant locus < 0.5
variant_counts = outDF['variant_grch38'].value_counts()
outDF = outDF[outDF['variant_grch38'].isin(set(variant_counts[variant_counts == 1].index))]
outDF = outDF[(outDF['end_distance_pct'] < 0.5) & (outDF['skipping_frequency'] < 0.5)]

# Filter outDF for variant calls in undiagnosed cases
outDF = outDF[outDF['patient_id'].isin(set(undiagnosed))]

# Save outDF to outfile
outDF = outDF.drop(['gnomad_nhomalt_grpmax_joint', 'homopolymer'], axis = 1)
outDF.to_csv(outfile, sep = '\t', index = False)
