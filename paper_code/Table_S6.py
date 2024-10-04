#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2024.10.03
Supplementary Table 6

Haplotype-resolved read counts for genes with heterozygous pathogenic variants in previously diagnosed patients. P values 
are computed from a haplotype dosage outlier test, which assesses whether haplotype expression ratios in a sample are 
consistent with those observed in tissue-matched GTEx controls (Methods).
'''

# =====================================================================================================================
#                                                       PACKAGES 
# =====================================================================================================================

# Load required packages
import os
import numpy as np
import pandas as pd

# =====================================================================================================================
#                                                   HELPER FUNCTIONS
# =====================================================================================================================

def QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes):
    results = [patient_id, gene_panel, gene_name]
    currdir = workdir + '/' + gene_panel.split('-')[0] + '/' + patient_id + '/RNA/stripe/target_genes/' + gene_name
    inDF = pd.read_csv(currdir + '/haplotype_dosage/output.txt', sep = '\t', header = None)
    vcfDF = pd.read_csv(currdir + '/gene.vcf.gz', sep = '\t', header = None, comment = '#', compression = 'gzip')
    hap1Variants = set(vcfDF[vcfDF[9] == '1|0'].apply(lambda x: x[0] + ':' + str(x[1]) + ' ' + x[3] + '>' + x[4], axis = 1))
    hap2Variants = set(vcfDF[vcfDF[9] == '0|1'].apply(lambda x: x[0] + ':' + str(x[1]) + ' ' + x[3] + '>' + x[4], axis = 1))
    if variants[0] in hap1Variants or variants[1] in hap2Variants:
        results += [variants[0], classes[0], int(inDF.iloc[0, 1]), variants[1], classes[1], int(inDF.iloc[1, 1]),
            format(inDF.iloc[2, 1], '.4f'), format(inDF.iloc[3, 1], '.3e')]
    elif variants[0] in hap2Variants or variants[1] in hap1Variants:
        results += [variants[1], classes[1], int(inDF.iloc[0, 1]), variants[0], classes[0], int(inDF.iloc[1, 1]),
            format(inDF.iloc[2, 1], '.4f'), format(inDF.iloc[3, 1], '.3e')]
    else:
        vcfDF = pd.read_csv(currdir + '/variant_calling/bcftools.vcf.gz', sep = '\t', header = None, comment = '#', compression = 'gzip')
        vcfDF['H1'], vcfDF['H2'] = vcfDF[9].str.split(':').str[0], vcfDF[10].str.split(':').str[0]
        hap1Variants = set(vcfDF[vcfDF['H1'] == '1'].apply(lambda x: x[0] + ':' + str(x[1]) + ' ' + x[3] + '>' + x[4], axis = 1))
        hap2Variants = set(vcfDF[vcfDF['H2'] == '1'].apply(lambda x: x[0] + ':' + str(x[1]) + ' ' + x[3] + '>' + x[4], axis = 1))
        if variants[0] in hap1Variants or variants[1] in hap2Variants:
            results += [variants[0], classes[0], int(inDF.iloc[0, 1]), variants[1], classes[1], int(inDF.iloc[1, 1]),
                format(inDF.iloc[2, 1], '.4f'), format(inDF.iloc[3, 1], '.3e')]
        elif variants[0] in hap2Variants or variants[1] in hap1Variants:
            results += [variants[1], classes[1], int(inDF.iloc[0, 1]), variants[0], classes[0], int(inDF.iloc[1, 1]),
                format(inDF.iloc[2, 1], '.4f'), format(inDF.iloc[3, 1], '.3e')]
    if len(results) == 3:
        return results + ['.'] * 8
    else:
        return results

# =====================================================================================================================
#                                                         MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir = '/mnt/isilon/lin_lab_share/STRIPE'
outfile = workdir + '/manuscript/Supplementary_Tables/Table_S6/Table_S6.txt'

# Create an empty data frame with columns of interest
outDF = pd.DataFrame(columns = ['patient_id', 'gene_panel', 'gene_name', 'hap1_variant_grch38', 'hap1_variant_class',
    'hap1_read_count', 'hap2_variant_grch38', 'hap2_variant_class', 'hap2_read_count', 'phasing_quality', 'p_value'])

# =====================================================================================================================
#                                                     11547 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = '11547', 'CDG-466', 'PMM2'
variants, classes = ['chr16:8811146 G>A', 'chr16:8811153 G>A'], ['Missense', 'Missense']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                     AnJa (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'AnJa', 'CDG-466', 'PIGN'
variants, classes = ['chr18:62106985 C>G', 'chr18:62113129 C>T'], ['Splice donor', 'Splice donor']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-132-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-132-1', 'CDG-466', 'PIGQ'
variants, classes = ['chr16:576255 G>A', 'chr16:578911 TCTA>T'], ['Splice donor', 'Inframe deletion']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-137-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-137-1', 'CDG-466', 'PIGN'
variants, classes = ['chr18:62106985 C>G', 'chr18:62143337 A>C'], ['Splice donor', 'Missense']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-147-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-147-1', 'CDG-466', 'PIGN'
variants, classes = ['chr18:62113277 C>T', 'chr18:62143306 C>T'], ['Missense', 'Splice donor']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-152-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-152-1', 'CDG-466', 'PMM2'
variants, classes = ['chr16:Δ8805603-8820989', 'chr16:8806398 C>T'], ['Structural deletion', 'Missense']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-161-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-161-1', 'CDG-466', 'PIGN'
variants, classes = ['chr18:62113129 C>T', 'chr18:62161218 ATCT>A'], ['Splice donor', 'Inframe deletion']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-183-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-183-1', 'CDG-466', 'PIGN'
variants, classes = ['chr18:62084583 C>T', 'chr18:62106797 G>A'], ['Missense', 'Nonsense']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-110-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-110-1', 'CDG-466', 'PIGN'
variants, classes = ['chr18:62095902 C>T', 'chr18:62157743 G>C'], ['Missense', 'Missense']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-125-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-125-1', 'CDG-466', 'PMM2'
variants, classes = ['chr16:8811088 C>A', 'chr16:8811153 G>A'], ['Missense', 'Missense']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-184-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-184-1', 'CDG-466', 'PIGN'
variants, classes = ['chr18:62045973 G>C', 'chr18:62143337 A>C'], ['Missense', 'Missense']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                    BS2-1 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'BS2-1', 'PMD-359', 'EARS2'
variants, classes = ['chr16:Δ23525826-23551747', 'chr16:23535176 C>T'], ['Structural deletion', 'Missense']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                   FH-2209 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'FH-2209', 'PMD-359', 'FH'
variants, classes = ['chr1:241508631 T>C', 'chr1:241519682 A>AG'], ['Missense', 'Frameshift duplication']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                   Q1513 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1513', 'PMD-359', 'FBXL4'
variants, classes = ['chr6:98874306 A>T', 'chr6:98880546 CACTT>C'], ['Missense', 'Splice donor']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                   Q1642 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1642', 'PMD-359', 'TRMU'
variants, classes = ['chr22:Δ46334556-46336380', 'chr22:46356055 G>A'], ['Structural deletion', 'Missense']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                   Q1687 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1687', 'PMD-359', 'NUBPL'
variants, classes = ['chr14:31562125 G>A', 'chr14:31826715 G>A'], ['Missense', 'Splice donor']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                   Q1695 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1695', 'PMD-359', 'ELAC2'
variants, classes = ['chr17:13005072 TGGA>T', 'chr17:13016933 CT>A'], ['Inframe deletion', 'Splice acceptor']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                   Q1819 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1819', 'PMD-359', 'TRIT1'
variants, classes = ['chr1:39850120 T>C', 'chr1:39883470 G>A'], ['Splice donor', 'Nonsense']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# =====================================================================================================================
#                                                  Q2032s1 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q2032s1', 'PMD-359', 'OPA1'
variants, classes = ['chr3:193643448 AG>A', '.'], ['Splice donor', '.']

currdir = workdir + '/' + gene_panel.split('-')[0] + '/' + patient_id + '/RNA/stripe/target_genes/' + gene_name
inDF = pd.read_csv(currdir + '/haplotype_dosage/output.txt', sep = '\t', header = None)
results = [patient_id, gene_panel, gene_name, variants[0], classes[0], int(inDF.iloc[1, 1]), variants[1], classes[1], 
    int(inDF.iloc[0, 1]), format(inDF.iloc[2, 1], '.4f'), format(inDF.iloc[3, 1], '.3e')]
outDF = pd.concat([outDF, pd.DataFrame([results], columns = outDF.columns)])

# =====================================================================================================================
#                                                  Q2319 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q2319', 'PMD-359', 'SUCLG1'
variants, classes = ['chr2:84425503 C>G', 'chr2:84449761 T>C'], ['Missense', 'Splice acceptor']

outDF = pd.concat([outDF, pd.DataFrame([QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes)],
    columns = outDF.columns)])

# Save outDF to outfile
outDF.to_csv(outfile, sep = '\t', index = False)
