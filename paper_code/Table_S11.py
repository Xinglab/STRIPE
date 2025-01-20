#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2025.01.14
Supplementary Table 11

Aberrant splice junctions in target genes that were identified and prioritized from TEQUILA-seq data for cohort individuals.
We filtered splice junctions in target genes based on the following sequential criteria:
    (i) feature at least one annotated splice site for the target gene
    (ii) have a minimum of 20 reads covering the corresponding splice sites in a given cohort individual and in at least 100 GTEx controls
    (iii) supported by a minimum of 20 reads in a given cohort individual
    (iv) called an aberrant splice junction (FDR < 1%, > 10% shift from the average usage frequency observed in GTEx controls) in no more than 3 cohort individuals
We also report a "phasing quality score" to represent the proportion of reads that could be confidently assigned to a haplotype.
P values were computed using an outlier splice junction test, which assesses whether the usage frequency of a splice junction in 
a sample is consistent with usage frequencies observed in tissue-matched GTEx controls (Methods).
'''

'''
We used the following Bash code to submit individual jobs for running our aberrant splice junction detection test per target gene:

WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"

# Iterate over genes on CDG-466 gene panel
cut -f5 "$WORKDIR/CDG/references/target_genes.bed" | while read GENE; do
    mkdir -p "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S11/CDG-466"
    mkdir -p "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S11/CDG-466/$GENE"
    cd "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S11/CDG-466/$GENE"
    SCRIPT="$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S11/CDG-466/$GENE/CDG-466_${GENE}.sh"
    echo -e "#!""/bin/bash" > "$SCRIPT"
    echo -e "source /scr1/users/wangr5/miniconda3/bin/activate sandbox" >> "$SCRIPT"
    echo -e "python $WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S11/Table_S11_Preprocessing.py CDG-466 $GENE" >> "$SCRIPT"

    sbatch --time=24:00:00 --mem=20G "$SCRIPT"
    sleep 1
done

# Iterate over genes on PMD-359 gene panel
cut -f5 "$WORKDIR/PMD/references/target_genes.bed" | while read GENE; do
    mkdir -p "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S11/PMD-359"
    mkdir -p "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S11/PMD-359/$GENE"
    cd "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S11/PMD-359/$GENE"
    SCRIPT="$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S11/PMD-359/$GENE/PMD-359_${GENE}.sh"
    echo -e "#!""/bin/bash" > "$SCRIPT"
    echo -e "source /scr1/users/wangr5/miniconda3/bin/activate sandbox" >> "$SCRIPT"
    echo -e "python $WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S11/Table_S11_Preprocessing.py PMD-359 $GENE" >> "$SCRIPT"

    sbatch --time=24:00:00 --mem=20G "$SCRIPT"
    sleep 1
done
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
outfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S11/Table_S11.txt'

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

# Iterate over all genes in CDG-466 gene panel
targets = pd.read_csv(workdir + '/CDG/references/target_genes.bed', sep = '\t', header = None)

for target in targets.values:
    infile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S11/CDG-466/' + target[4] + '/output.tsv'
    inDF = pd.read_csv(infile, sep = '\t', header = 0).fillna('')
    inDF = inDF[(inDF['splice_site_coverage'] >= 20) & (inDF['junction_read_count'] >= 20)]
    output += [list(item) for item in inDF.values]

# Iterate over all genes in PMD-359 gene panel
targets = pd.read_csv(workdir + '/PMD/references/target_genes.bed', sep = '\t', header = None)

for target in targets.values:
    infile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S11/PMD-359/' + target[4] + '/output.tsv'
    inDF = pd.read_csv(infile, sep = '\t', header = 0).fillna('')
    inDF = inDF[(inDF['splice_site_coverage'] >= 20) & (inDF['junction_read_count'] >= 20)]
    output += [list(item) for item in inDF.values]

# Convert output into a dataframe
outDF = pd.DataFrame(output, columns = ['id', 'status', 'disease_group', 'category', 'gene_panel', 'gene_name', 'phasing_quality', 'read_group', 
    'junction_coord', 'junction_read_count', 'splice_site_coverage', 'gtex_junction_usage', 'p_value'])

# Filter out aberrant splice junctions in target genes (FDR < 1%, > 10% shift from the average usage frequency observed in GTEx controls) that were called in more than 3 cohort individuals
outDF['fdr'] = fdrcorrection(outDF['p_value'], alpha = 0.01, is_sorted = False)[1]
filterDF = outDF[outDF['fdr'] < 0.01]
filterDF['shift'] = (filterDF['junction_read_count']/filterDF['splice_site_coverage']) - filterDF['gtex_junction_usage']
filterDF = filterDF[filterDF['shift'] > 0.1]

countDF = filterDF[['id', 'gene_panel', 'gene_name', 'junction_coord']].drop_duplicates()
junction_counts = (countDF['gene_panel'] + '_' + countDF['gene_name'] + '_' + countDF['junction_coord']).value_counts()
filterDF = filterDF[(filterDF['gene_panel'] + '_' + filterDF['gene_name'] + '_' + filterDF['junction_coord']).map(dict(zip(junction_counts.index, 
    junction_counts.values))) <= 3]

# Sort rows of filterDF by id (using the order in metaDF) and gene name
filterDF['id_num'] = filterDF['id'].map(dict(zip(metaDF.ID, metaDF.index)))
filterDF = filterDF.sort_values(by = ['id_num', 'gene_name'])

# Save filterDF to outfile
filterDF.drop(['fdr', 'shift', 'id_num'], axis = 1).to_csv(outfile, sep = '\t', index = False)

# Calculate median number of aberrant splice junctions in target genes prioritized among 44 cohort individuals for whom prior genetic testing revealed only a single 
# pathogenic variant in a gene associated with a recessive condition (n = 6) or no candidate disease-causing variants (n = 38)
filterDF['junction_id'] = filterDF['gene_panel'] + '_' + filterDF['gene_name'] + '_' + filterDF['junction_coord']
counts = filterDF[['id', 'junction_id']].drop_duplicates()['id'].value_counts()
np.median([dict(zip(counts.index, counts.values)).get(sample_id, 0) for sample_id in metaDF[metaDF['Category'].isin({'Single variant, recessive condition',
    'No candidate variants'})]['ID']]) # 7 aberrant junctions (target genes) per patient
