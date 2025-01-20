#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2025.01.11
Supplementary Table 9

Rare deleterious variants in target genes identified and prioritized from TEQUILA-seq data for cohort individuals. 
For each variant, we report the following quality control metrics:
* Phasing quality: Value from 0 to 1 that represents the proportion of reads that could be confidently assigned to a haplotype
* Homopolymer: TRUE if the variant lies within a homopolymer region (i.e. stretch of identical nucleotides ≥ 4 bp in length within a 9 bp window centered around the variant position)
    * We tend to observe false positive variant calls within homopolymer regions
* End distance percentage: Percentage of reads supporting the variant in which the variant is found at the end of the read (i.e. last 30 nucleotides)
    * Variants corresponding to read alignment artifacts often have a high "end distance percentage"
* Skipping frequency: Percentage of reads covering variant locus that are intronic
    * Variants corresponding to read alignment artifacts (or variants present on extremely minor isoforms) tend to have a high skipping frequency
We also only retained variants that were identified and prioritized in no more than two cohort individuals.
'''

# =====================================================================================================================
#                                                       PACKAGES 
# =====================================================================================================================

# Load required packages
import os, re
import pandas as pd
import numpy as np
import hgvs.parser, hgvs.assemblymapper, hgvs.dataproviders.uta

# =====================================================================================================================
#                                                         MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir = '/mnt/isilon/lin_lab_share/STRIPE'
outfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S9/Table_S9.txt'
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

# Read in data frame mapping Ensembl gene IDs to RefSeq MANE select transcript IDs
mappingDF = pd.read_csv('/scr1/users/wangr5/references/gene2refseq.txt.gz', sep = '\t', header = 0, compression = 'gzip')
mapDict = dict(zip(mappingDF.iloc[:,0], mappingDF.iloc[:,2]))

# Read in data frame mapping RefSeq accessions to chromosome numbers
assemblyDF = pd.read_csv('/scr1/users/wangr5/references/sequence_report.tsv', sep = '\t', header = 0)
assemblyDF = assemblyDF[assemblyDF['Role'] == 'assembled-molecule']
assemblyDict = dict(zip('chr' + assemblyDF['Chromosome name'], assemblyDF['RefSeq seq accession']))

# Initialize an HGVS parser and variant mapper
hp, hdp = hgvs.parser.Parser(), hgvs.dataproviders.uta.connect()
am38 = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name = 'GRCh38')

# Iterate over variant calling results from STRIPE for all CDG patients
targets = pd.read_csv(workdir + '/CDG/references/target_genes.bed', sep = '\t', header = None)
targets[6], targets[7] = targets[3].map(mapDict), targets[0].map(assemblyDict)
sampleDF = pd.read_csv(workdir + '/CDG/samples.txt', sep = '\t', header = 0)
sampleDF = sampleDF[(sampleDF['Provider'] != 'Lan Lin') & (sampleDF['ID'] != 'E1877')]
cdg_cohort = metaDF[metaDF['ID'].isin(set(sampleDF['ID']))]

for patient in cdg_cohort.values:
    for target in targets.values:
        infile = workdir + '/CDG/' + patient[0] + '/RNA/stripe/target_genes/' + target[4] + '/variant_calling/merged_variants.tsv'
        if os.path.isfile(infile):
            inDF = pd.read_csv(infile, sep = '\t', header = 0)
            inDF = inDF[(inDF['GNOMAD_AF_GRPMAX_JOINT'] < 0.01) & ~inDF['CLINVAR'].str.contains('enign')]
            phasing_file = workdir + '/CDG/' + patient[0] + '/RNA/stripe/target_genes/' + target[4] + '/phasing_stats.txt'
            if os.path.isfile(phasing_file):
                phaseDF = pd.read_csv(phasing_file, sep = '\t', header = None)
                phasing_score = phaseDF.iloc[0, 1]
            else:
                phasing_score = ''
            output += [[patient[0], patient[1], patient[2], patient[3], 'CDG-466', target[4], item[0] + ':' + str(item[1]) + ' ' + 
                item[2] + '>' + item[3], target[7] + ':g.' + str(item[1]) + item[2] + '>' + item[3], target[6]] + list(item[4:-3]) + 
                [phasing_score] + list(item[-3:]) for item in inDF.values]

# Iterate over variant calling results from STRIPE for all PMD patients
targets = pd.read_csv(workdir + '/PMD/references/target_genes.bed', sep = '\t', header = None)
targets[6], targets[7] = targets[3].map(mapDict), targets[0].map(assemblyDict)
sampleDF = pd.read_csv(workdir + '/PMD/samples.txt', sep = '\t', header = 0)
sampleDF = sampleDF[(sampleDF['Provider'] == 'Rebecca Ganetzky') | ((sampleDF['Provider'] == 'Marni Falk') & 
    (sampleDF['Status'] != 'Diagnosed') & ~sampleDF['ID'].str.contains('-MF'))]
pmd_cohort = metaDF[metaDF['ID'].isin(set(sampleDF['ID']))]

for patient in pmd_cohort.values:
    for target in targets.values:
        infile = workdir + '/PMD/' + patient[0] + '/RNA/stripe/target_genes/' + target[4] + '/variant_calling/merged_variants.tsv'
        if os.path.isfile(infile):
            inDF = pd.read_csv(infile, sep = '\t', header = 0)
            inDF = inDF[(inDF['GNOMAD_AF_GRPMAX_JOINT'] < 0.01) & ~inDF['CLINVAR'].str.contains('enign')]
            phasing_file = workdir + '/PMD/' + patient[0] + '/RNA/stripe/target_genes/' + target[4] + '/phasing_stats.txt'
            if os.path.isfile(phasing_file):
                phaseDF = pd.read_csv(phasing_file, sep = '\t', header = None)
                phasing_score = phaseDF.iloc[0, 1]
            else:
                phasing_score = ''
            output += [[patient[0], patient[1], patient[2], patient[3], 'PMD-359', target[4], item[0] + ':' + str(item[1]) + ' ' + 
                item[2] + '>' + item[3], target[7] + ':g.' + str(item[1]) + item[2] + '>' + item[3], target[6]] + list(item[4:-3]) + 
                [phasing_score] + list(item[-3:]) for item in inDF.values]

# Convert output into a dataframe
outDF = pd.DataFrame(output, columns = ['id', 'status', 'disease_group', 'category', 'gene_panel', 'gene_name', 'variant_grch38',
    'variant_hgvs', 'mane_select'] + [tool + '_' + field for tool in ['longcallR', 'clair3', 'deepvariant', 'phasing'] for field in ['dp', 'af', 'gt']] + 
    ['cadd', 'spliceai', 'gnomad_af_grpmax_joint', 'gnomad_nhomalt_grpmax_joint', 'clinvar', 'phasing_quality', 'homopolymer', 
    'end_distance_pct', 'skipping_frequency'])

# Only keep variants called in no more than 2 individuals
# Variants appearing in multiple individuals likely represent variant calling artifacts since the remaining variants are supposed to be rare in the population
variant_counts = outDF['variant_grch38'].value_counts()
filterDF = outDF[outDF['variant_grch38'].map(dict(zip(variant_counts.index, variant_counts.values))) <= 2]

# Retrieve all unique variants represented in filterDF and convert them to HGVS nomenclature
keepVariants = filterDF[['variant_grch38', 'variant_hgvs', 'mane_select']].drop_duplicates().values
hgvsVariants, ctr = [], 1
for kv in keepVariants:
    ref, alt = re.split('\d+', kv[1])[-1].split('>')
    if len(ref) == 1 and len(alt) == 1:
        queryVariant = kv[1]
    elif len(ref) > 1 and len(alt) == 1 and ref[0] == alt:
        position = int(re.findall(r'\d+', kv[1])[-1])
        if len(ref) - len(alt) == 1:
            queryVariant = kv[1].split(':')[0] + ':g.' + str(position+1) + 'del'
        else:
            queryVariant = kv[1].split(':')[0] + ':g.' + str(position+1) + '_' + str(position+len(ref)-len(alt)) + 'del'
    elif len(ref) == 1 and len(alt) > 1 and ref == alt[0]:
        position = int(re.findall(r'\d+', kv[1])[-1])
        queryVariant = kv[1].split(':')[0] + ':g.' + str(position) + '_' + str(position+1) + 'ins' + alt[1:]
    else:
        queryVariant = None
    if queryVariant is not None and isinstance(kv[2], str):
        v = hp.parse_hgvs_variant(queryVariant)
        relevant_tx = am38.relevant_transcripts(v)
        mane_select = [tx for tx in relevant_tx if tx.split('.')[0] == kv[2].split('.')[0]]
        if len(mane_select) == 1:
            cVar = am38.g_to_c(v, mane_select[0])
            pVar = am38.c_to_p(cVar)
            if '?' not in str(pVar):
                hgvsVariants.append([kv[0], str(cVar) + ', ' + str(pVar).split(':')[1].replace('(', '').replace(')', '')])
            else:
                hgvsVariants.append([kv[0], str(cVar)])
        else:
            hgvsVariants.append([kv[0], '.'])
    else:
        hgvsVariants.append([kv[0], '.'])
    print('Done with ' + str(ctr) + ' / ' + str(len(keepVariants)) + ' variants...')
    ctr += 1

# Save filterDF to outfile
filterDF['variant_hgvs'] = filterDF['variant_grch38'].map({item[0]: item[1] for item in hgvsVariants}).fillna('.')
filterDF.drop('mane_select', axis = 1).to_csv(outfile, sep = '\t', index = False)

# Calculate median number of variants prioritized among 44 cohort individuals for whom prior genetic testing revealed only a single 
# pathogenic variant in a gene associated with a recessive condition (n = 6) or no candidate disease-causing variants (n = 38)
counts = filterDF['id'].value_counts()
np.median([dict(zip(counts.index, counts.values)).get(sample_id, 0) for sample_id in metaDF[metaDF['Category'].isin({'Single variant, recessive condition',
    'No candidate variants'})]['ID']]) # 7 variants per patient
