#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2024.10.04
Supplementary Table 7

Haplotype-resolved splice junctions (supported by at least 20 reads) in genes with heterozygous or hemizgyous pathogenic 
variants in previously diagnosed patients. P values are computed from a outlier splice junction test, which assesses whether 
the usage frequency of a splice junction in a sample is consistent with usage frequencies observed in tissue-matched 
GTEx controls (Methods).
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

def RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome):
    _ = os.system('mkdir -p ' + workdir + '/manuscript/Supplementary_Tables/Table_S7/tmp')
    _ = os.system('echo "##fileformat=VCFv4.2" > ' + workdir + '/manuscript/Supplementary_Tables/Table_S7/tmp/input.vcf')
    _ = os.system(bcftools + ' view -h ' + workdir + '/' + gene_panel.split('-')[0] + '/' + patient_id + 
        '/RNA/stripe/target_genes/longcallR/output.vcf.gz | grep "##contig=" >> ' + workdir +
        '/manuscript/Supplementary_Tables/Table_S7/tmp/input.vcf')
    varDF = pd.DataFrame([[item for item in var.replace(':', '_').replace(' ', '_').replace('>', '_').split('_')] for var in variants],
        columns = ['#CHROM', 'POS', 'REF', 'ALT'])
    varDF['ID'], varDF['QUAL'], varDF['FILTER'], varDF['INFO'] = '.', '.', '.', '.'
    varDF[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']].to_csv(workdir +
        '/manuscript/Supplementary_Tables/Table_S7/tmp/input.vcf', mode = 'a', sep = '\t', index = False)
    _ = os.system('(spliceai -I ' + workdir + '/manuscript/Supplementary_Tables/Table_S7/tmp/input.vcf -O ' +
        workdir + '/manuscript/Supplementary_Tables/Table_S7/tmp/output.vcf -R ' + genome + ' -A grch38) > /dev/null 2>&1')
    if os.path.isfile(workdir + '/manuscript/Supplementary_Tables/Table_S7/tmp/output.vcf'):
        spliceai_results = pd.read_csv(workdir + '/manuscript/Supplementary_Tables/Table_S7/tmp/output.vcf', header = None, sep = '\t',
            names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'], comment = '#')
        spliceai_results['SPLICEAI'] = spliceai_results['INFO'].apply(lambda x: np.max([float(val) if val != '.' else 0 for val in x.split('|')[2:6]]) 
            if 'SpliceAI' in x else np.nan)
        spliceai_dict = dict(zip(spliceai_results['CHROM'] + ':' + spliceai_results['POS'].astype(str) + ' ' + spliceai_results['REF'] + '>' + 
            spliceai_results['ALT'], spliceai_results['SPLICEAI']))
        _ = os.system('rm ' + workdir + '/manuscript/Supplementary_Tables/Table_S7/tmp/input.vcf && rm ' + workdir + 
            '/manuscript/Supplementary_Tables/Table_S7/tmp/output.vcf && rmdir ' + workdir + '/manuscript/Supplementary_Tables/Table_S7/tmp')
        return [spliceai_dict.get(var, np.nan) for var in variants]
    else:
        _ = os.system('rm ' + workdir + '/manuscript/Supplementary_Tables/Table_S7/tmp/input.vcf && rmdir ' + workdir + 
            '/manuscript/Supplementary_Tables/Table_S7/tmp')
        return [np.nan] * length(variants)

def QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, spliceai_scores):
    results, order = [patient_id, gene_panel, gene_name], []
    currdir = workdir + '/' + gene_panel.split('-')[0] + '/' + patient_id + '/RNA/stripe/target_genes/' + gene_name
    vcfDF = pd.read_csv(currdir + '/gene.vcf.gz', sep = '\t', header = None, comment = '#', compression = 'gzip')
    hap1Variants = set(vcfDF[vcfDF[9] == '1|0'].apply(lambda x: x[0] + ':' + str(x[1]) + ' ' + x[3] + '>' + x[4], axis = 1))
    hap2Variants = set(vcfDF[vcfDF[9] == '0|1'].apply(lambda x: x[0] + ':' + str(x[1]) + ' ' + x[3] + '>' + x[4], axis = 1))
    if variants[0] in hap1Variants or variants[1] in hap2Variants:
        order = ['hap1', 'hap2']
    elif variants[0] in hap2Variants or variants[1] in hap1Variants:
        order = ['hap2', 'hap1']
    else:
        vcfDF = pd.read_csv(currdir + '/variant_calling/bcftools.vcf.gz', sep = '\t', header = None, comment = '#', compression = 'gzip')
        vcfDF['H1'], vcfDF['H2'] = vcfDF[9].str.split(':').str[0], vcfDF[10].str.split(':').str[0]
        hap1Variants = set(vcfDF[vcfDF['H1'] == '1'].apply(lambda x: x[0] + ':' + str(x[1]) + ' ' + x[3] + '>' + x[4], axis = 1))
        hap2Variants = set(vcfDF[vcfDF['H2'] == '1'].apply(lambda x: x[0] + ':' + str(x[1]) + ' ' + x[3] + '>' + x[4], axis = 1))
        if variants[0] in hap1Variants or variants[1] in hap2Variants:
            order = ['hap1', 'hap2']
        elif variants[0] in hap2Variants or variants[1] in hap1Variants:
            order = ['hap2', 'hap1']
    if len(order) > 0:
        output = []
        return ([results + [variants[order.index('hap1')], classes[order.index('hap1')], spliceai_scores[order.index('hap1')]] + 
            list(item) for item in pd.read_csv(currdir + '/aberrant_splicing/output_hap1.txt', sep = '\t', header = 0).values] +
            [results + [variants[order.index('hap2')], classes[order.index('hap2')], spliceai_scores[order.index('hap2')]] + 
            list(item) for item in pd.read_csv(currdir + '/aberrant_splicing/output_hap2.txt', sep = '\t', header = 0).values])
    else:
        return [results + ['.'] * 10]
    
# =====================================================================================================================
#                                                         MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir = '/mnt/isilon/lin_lab_share/STRIPE'
bcftools = '/scr1/users/wangr5/tools/bcftools-1.21/bcftools'
genome = '/scr1/users/wangr5/references/GRCh38.primary_assembly.genome.fa'
outfile = workdir + '/manuscript/Supplementary_Tables/Table_S7/Table_S7.txt'

# Create an empty data frame with columns of interest
outDF = pd.DataFrame(columns = ['patient_id', 'gene_panel', 'gene_name', 'variant_grch38', 'variant_class', 'spliceai',
    'hap_junction_coord', 'hap_junction_count', 'hap_junction_coverage', 'hap_junction_usage_sample', 'hap_junction_usage_gtex', 
    'hap_junction_usage_shift', 'p_value']) 

# =====================================================================================================================
#                                                     11547 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = '11547', 'CDG-466', 'PMM2'
variants, classes = ['chr16:8811146 G>A', 'chr16:8811153 G>A'], ['Missense', 'Missense']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                     AnJa (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'AnJa', 'CDG-466', 'PIGN'
variants, classes = ['chr18:62106985 C>G', 'chr18:62113129 C>T'], ['Splice donor', 'Splice donor']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-132-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-132-1', 'CDG-466', 'PIGQ'
variants, classes = ['chr16:576255 G>A', 'chr16:578911 TCTA>T'], ['Splice donor', 'Inframe deletion']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-137-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-137-1', 'CDG-466', 'PIGN'
variants, classes = ['chr18:62106985 C>G', 'chr18:62143337 A>C'], ['Splice donor', 'Missense']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-147-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-147-1', 'CDG-466', 'PIGN'
variants, classes = ['chr18:62113277 C>T', 'chr18:62143306 C>T'], ['Missense', 'Splice donor']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-152-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-152-1', 'CDG-466', 'PMM2'
variants, classes = ['chr16:Δ8805603-8820989', 'chr16:8806398 C>T'], ['Structural deletion', 'Missense']

spliceai_scores = [np.nan] + RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants[1:], genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-161-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-161-1', 'CDG-466', 'PIGN'
variants, classes = ['chr18:62113129 C>T', 'chr18:62161218 ATCT>A'], ['Splice donor', 'Inframe deletion']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-183-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-183-1', 'CDG-466', 'PIGN'
variants, classes = ['chr18:62084583 C>T', 'chr18:62106797 G>A'], ['Missense', 'Nonsense']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                  ATP6AP2 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'ATP6AP2', 'CDG-466', 'ATP6AP2'
variants, classes = ['chrX:40589011 A>G'], ['Missense']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
currdir = workdir + '/' + gene_panel.split('-')[0] + '/' + patient_id + '/RNA/stripe/target_genes/' + gene_name
results = [[patient_id, gene_panel, gene_name] + [variants[0], classes[0], spliceai_scores[0]] + 
    list(item) for item in pd.read_csv(currdir + '/aberrant_splicing/output_all.txt', sep = '\t', header = 0).values]
outDF = pd.concat([outDF, pd.DataFrame(results, columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-110-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-110-1', 'CDG-466', 'PIGN'
variants, classes = ['chr18:62095902 C>T', 'chr18:62157743 G>C'], ['Missense', 'Missense']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-115-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-115-1', 'CDG-466', 'SLC35A2'
variants, classes = ['chrX:48909877 C>T'], ['Missense']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
currdir = workdir + '/' + gene_panel.split('-')[0] + '/' + patient_id + '/RNA/stripe/target_genes/' + gene_name
results = [[patient_id, gene_panel, gene_name] + [variants[0], classes[0], spliceai_scores[0]] + 
    list(item) for item in pd.read_csv(currdir + '/aberrant_splicing/output_hap1.txt', sep = '\t', header = 0).values]
outDF = pd.concat([outDF, pd.DataFrame(results, columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-125-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-125-1', 'CDG-466', 'PMM2'
variants, classes = ['chr16:8811088 C>A', 'chr16:8811153 G>A'], ['Missense', 'Missense']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                  CDG-184-1 (CDG)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-184-1', 'CDG-466', 'PIGN'
variants, classes = ['chr18:62045973 G>C', 'chr18:62143337 A>C'], ['Missense', 'Missense']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                    BS2-1 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'BS2-1', 'PMD-359', 'EARS2'
variants, classes = ['chr16:Δ23525826-23551747', 'chr16:23535176 C>T'], ['Structural deletion', 'Missense']

spliceai_scores = [np.nan] + RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants[1:], genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                   FH-2209 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'FH-2209', 'PMD-359', 'FH'
variants, classes = ['chr1:241508631 T>C', 'chr1:241519682 A>AG'], ['Missense', 'Frameshift duplication']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                   Q1513 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1513', 'PMD-359', 'FBXL4'
variants, classes = ['chr6:98874306 A>T', 'chr6:98880546 CACTT>C'], ['Missense', 'Splice donor']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                   Q1642 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1642', 'PMD-359', 'TRMU'
variants, classes = ['chr22:Δ46334556-46336380', 'chr22:46356055 G>A'], ['Structural deletion', 'Missense']

spliceai_scores = [np.nan] + RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants[1:], genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                   Q1687 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1687', 'PMD-359', 'NUBPL'
variants, classes = ['chr14:31562125 G>A', 'chr14:31826715 G>A'], ['Missense', 'Splice donor']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                   Q1695 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1695', 'PMD-359', 'ELAC2'
variants, classes = ['chr17:13005072 TGGA>T', 'chr17:13016933 CT>A'], ['Inframe deletion', 'Splice acceptor']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                   Q1819 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1819', 'PMD-359', 'TRIT1'
variants, classes = ['chr1:39850120 T>C', 'chr1:39883470 G>A'], ['Splice donor', 'Nonsense']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# =====================================================================================================================
#                                                  Q2032s1 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q2032s1', 'PMD-359', 'OPA1'
variants, classes = ['chr3:193643448 AG>A'], ['Splice donor']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
currdir = workdir + '/' + gene_panel.split('-')[0] + '/' + patient_id + '/RNA/stripe/target_genes/' + gene_name
results = [[patient_id, gene_panel, gene_name] + [variants[0], classes[0], spliceai_scores[0]] + 
    list(item) for item in pd.read_csv(currdir + '/aberrant_splicing/output_hap2.txt', sep = '\t', header = 0).values]
outDF = pd.concat([outDF, pd.DataFrame(results, columns = outDF.columns)])

# =====================================================================================================================
#                                                  Q2319 (PMD)
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q2319', 'PMD-359', 'SUCLG1'
variants, classes = ['chr2:84425503 C>G', 'chr2:84449761 T>C'], ['Missense', 'Splice acceptor']

spliceai_scores = RunSpliceAI(workdir, bcftools, gene_panel, patient_id, variants, genome)
outDF = pd.concat([outDF, pd.DataFrame(QueryVariant(workdir, gene_panel, patient_id, gene_name, variants, classes, 
    spliceai_scores), columns = outDF.columns)])

# Save outDF to outfile
outDF.to_csv(outfile, sep = '\t', index = False)
