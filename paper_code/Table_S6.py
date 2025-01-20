#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2024.12.29
Supplementary Table 6

Expression ratios for haplotypes carrying pathogenic variants based on patient TEQUILA-seq reads. 
P values were computed using a gene expression dosage outlier test, which assesses whether haplotype expression ratios 
in a sample are consistent with those observed in tissue-matched GTEx controls (Methods).
'''

# =====================================================================================================================
#                                                   HELPER FUNCTIONS
# =====================================================================================================================

def PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, dna_variants):
    _ = os.system('mkdir -p ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id)
    if os.path.isdir(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name):
        _ = os.system('rm -rf ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name)
    _ = os.system('mkdir -p ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name)
    targetDF = pd.read_csv(workdir + '/' + gene_panel.split('-')[0] + '/references/target_genes.bed', sep = '\t', header = None)
    target = targetDF[targetDF[4] == gene_name].values[0]
    command = bcftools + ' view -H -r ' + target[0] + ':' + str(target[1]+1) + '-' + str(target[2]) + ' -m2 -M2 '
    rna_variants = pd.DataFrame([item.split() for item in subprocess.run(command + workdir + '/' + gene_panel.split('-')[0] +
        '/' + patient_id + '/RNA/stripe/target_genes/longcallR/output.vcf.gz', shell = True, capture_output = True, 
        text = True).stdout.splitlines()], columns = list(range(10)))
    rna_variants['RNA'] = [PullFeature(item, 'GT') for item in rna_variants[[8, 9]].values.tolist()]
    rna_variants['DP'] = [PullFeature(item, 'DP') for item in rna_variants[[8, 9]].values.tolist()]
    rna_variants['AF'] = [PullFeature(item, 'AF').split(',')[0] for item in rna_variants[[8, 9]].values.tolist()]
    rna_variants['PS'] = [PullFeature(item, 'PS') for item in rna_variants[[8, 9]].values.tolist()]
    rna_variants['Dense'] = CheckDensity(sorted(rna_variants[1].astype(int).tolist()))
    rna_variants = pd.concat([rna_variants[(rna_variants[6] == 'PASS') & ~rna_variants['Dense']], 
        rna_variants[rna_variants[6].isin({'LowQual', 'RnaEdit'}) & (rna_variants['DP'].astype(int) >= 100) & 
        (rna_variants['AF'].astype(float) >= 0.3) & ~rna_variants['Dense']]])
    rna_variants = rna_variants[rna_variants['RNA'].isin({'0/1', '1|0', '0|1'})][[0, 1, 3, 4, 'RNA', 'PS']]
    rna_variants.loc[rna_variants['RNA'] == '0/1', 'PS'] = np.nan
    if dna_variants is not None:
        merge_variants = MergeVariants(rna_variants, dna_variants)
    else:
        merge_variants = MergeVariants(rna_variants, pd.DataFrame(columns = list(range(10)) + ['DNA']))
    _ = os.system(samtools + ' view -hb -F 256 -q 1 ' + workdir + '/' + gene_panel.split('-')[0] +
        '/' + patient_id + '/RNA/' + patient_id + '_TEQUILA.bam' + ' ' + target[0] + ':' + str(target[1]+1) + 
        '-' + str(target[2]) + ' > ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + 
        '/' + gene_name + '/all_reads.bam && ' + samtools + ' index ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
        patient_id + '/' + gene_name + '/all_reads.bam')
    if merge_variants.shape[0] > 0:
        samfile = pysam.AlignmentFile(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + 
            '/' + gene_name + '/all_reads.bam', mode = 'rb')
        merge_variants['Coverage'] = merge_variants.apply(lambda x: np.sum(samfile.count_coverage(x[0], int(x[1])-1, int(x[1]))), axis = 1)
        merge_variants = merge_variants[merge_variants['Coverage'] >= 20]
        if merge_variants.shape[0] > 0:
            bestPS = merge_variants[['PS', 'Coverage']].groupby('PS').max().reset_index().sort_values(by = 'Coverage', ascending = False).head(1)['PS'].item()
            merge_variants = merge_variants[merge_variants['PS'] == bestPS][[0, 1, 3, 4, 'RNA']]
            MakeHaplotypeFasta(samtools, bcftools, merge_variants, genome, target[0] + ':' + str(target[1]+1) + '-' + str(target[2]), 
                workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name)
            _ = os.system('(' + samtools + ' fastq -n ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + 
                '/' + gene_name + '/all_reads.bam | gzip > ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + 
                '/' + gene_name + '/all_reads.fastq.gz) > /dev/null 2>&1')
            _ = os.system('(' + minimap2 + ' -ax splice --secondary=no ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
                patient_id + '/' + gene_name + '/gene.fa ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + 
                '/' + gene_name + '/all_reads.fastq.gz > ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + 
                '/' + gene_name + '/all_reads.sam) > /dev/null 2>&1')
            _ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + 
                '/all_reads.fastq.gz && rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
                gene_name + '/gene.fa && rm ' +  workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
                gene_name + '/gene.fa.fai')
            assignments = HaplotypeAssignment(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.sam')
            _ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.sam')
            if len({key for key in assignments if assignments[key] != 'NA'}) >= 100:
                qscore = PhasingQuality(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + 
                    '/all_reads.bam', assignments, target[1], target[2])
                with open(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/phasing_stats.txt', 'w') as quality_file:
                    _ = quality_file.write('Phasing quality score\t%s\n' % format(qscore, '.4f'))
                    _ = quality_file.write('Read count (haplotype 1)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap1'}))
                    _ = quality_file.write('Read count (haplotype 2)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap2'}))
                    _ = quality_file.write('Read count (unassigned)\t%s\n' % len({key for key in assignments if assignments[key] == 'NA'}))
                WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap1'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
                    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
                    gene_name + '/hap1_reads.bam')
                WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap2'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
                    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
                    gene_name + '/hap2_reads.bam')
                WriteBam(samtools, {key for key in assignments if assignments[key] == 'NA'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
                    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
                    gene_name + '/unassigned_reads.bam')
                _ = os.system('(' + bcftools + ' mpileup -Ou -f ' + genome + ' -r ' + target[0] + ':' + str(target[1]+1) + '-' + 
                    str(target[2]) + ' -B -Q 5 ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + 
                    '/hap1_reads.bam ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + 
                    '/hap2_reads.bam | ' + bcftools + ' call -mv -Oz ' + '--ploidy 1 -o ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
                    patient_id + '/' + gene_name + '/bcftools.vcf.gz) > /dev/null 2>&1')
            else:
                _ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/gene.vcf.gz && rm ' + 
                    workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/gene.vcf.gz.tbi')
    _ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.bam && rm ' + 
        workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.bam.bai')

def DosageOutlierTest(workdir, patient_id, gene_panel, gene_name):
    gtexDF = pd.read_csv(workdir + '/' + gene_panel.split('-')[0] + '/references/GTEx_v8/Fibroblast/target_genes.haplotype_expression.txt', sep = '\t', header = 0)
    gtexDF.loc[:, 'name'] = gtexDF['name'].str.split('.').str[0]
    phasingInfo = list(pd.read_csv(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/phasing_stats.txt', 
        sep = '\t', header = None).iloc[:3, 1])
    targetDF = pd.read_csv(workdir + '/' + gene_panel.split('-')[0] + '/references/target_genes.bed', sep = '\t', header = None)
    ctrlDF = gtexDF[gtexDF['name'] == targetDF[targetDF[4] == gene_name][3].item()].iloc[:, 1:]
    ctrlCounts = [tuple(map(int, item.split('|'))) for item in ctrlDF.values[0]]
    ctrlRatios = [min(item)/sum(item) if sum(item) >= 20 else np.nan for item in ctrlCounts]
    ctrlRatios = MaskOutliers(np.array(ctrlRatios))
    ctrlRatios = np.round(RescaleValues(ctrlRatios, 0, 1, 1e-3, 1-1e-3), 3)
    myFunc = lambda x: 2*(digamma(x)-digamma(2*x)) - np.nanmean(np.log(ctrlRatios*(1-ctrlRatios)))
    fit = brentq(myFunc, 1, 999, full_output = True, disp = False)
    pval = BetaBinomialTest(min(phasingInfo[1], phasingInfo[2]), phasingInfo[1] + phasingInfo[2], [fit[0]] * 2)
    return pval

# =====================================================================================================================
#                                                         MAIN
# =====================================================================================================================

import sys

# Establish working directories and file paths
workdir = '/mnt/isilon/lin_lab_share/STRIPE'
bcftools = '/scr1/users/wangr5/tools/bcftools-1.21/bcftools'
samtools = '/scr1/users/wangr5/tools/samtools-1.21/samtools'
minimap2 = '/scr1/users/wangr5/tools/minimap2-2.26_x64-linux/minimap2'
genome = '/scr1/users/wangr5/references/GRCh38.primary_assembly.genome.fa'

# Import functions used in Run_STRIPE.py
sys.path.append(workdir + '/scripts/main')
from Run_STRIPE import *

# =====================================================================================================================
#                                                       BS2-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'BS2-1', 'PMD-359', 'EARS2'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name)    

# =====================================================================================================================
#                                                       Q1695
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1695', 'PMD-359', 'ELAC2'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Try phasing reads based the following heterozygous variants identified from previous testing:
#   * chr17:13005072 TGGA>T
#   * chr17:13016933 CT>A
_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/*')
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2,
    pd.DataFrame([['chr17', '13005072', '.', 'TGGA', 'T', '.', 'PASS', '.', 'GT', '1|0', '1|0'],
    ['chr17', '13016933', '.', 'CT', 'A', '.', 'PASS', '.', 'GT', '0|1', '0|1']], columns = list(range(10)) + ['DNA']))

# From our initial round of read phasing, reads assigned to haplotype 1 featured usage of the canonical acceptor splice site 
# chr17:13016933 while reads assigned to haplotype 2 featured usage of a cryptic acceptor splice site chr17:13016926
# Many unphased reads appeared to either use the canonical acceptor splice site chr17:13016933 (haplotype 1) or a cryptic
# acceptor splice site chr17:13016926 (haplotype 2). Based on this information, we will adopt the following strategy for phasing reads:
#   Haplotype 1: Reads that harbor the 3 bp in frame deletion chr17:13005072 TGGA>T OR feature the canonical acceptor splice site chr17:13016933
#   Haplotype 2: Reads that feature the cryptic acceptor splice site chr17:13016926
_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/*')
targetDF = pd.read_csv(workdir + '/' + gene_panel.split('-')[0] + '/references/target_genes.bed', sep = '\t', header = None)
target, assignments = targetDF[targetDF[4] == gene_name].values[0], dict()
_ = os.system(samtools + ' view -hb -F 256 -q 1 ' + workdir + '/' + gene_panel.split('-')[0] +
    '/' + patient_id + '/RNA/' + patient_id + '_TEQUILA.bam' + ' ' + target[0] + ':' + str(target[1]+1) + 
    '-' + str(target[2]) + ' > ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + 
    '/' + gene_name + '/all_reads.bam && ' + samtools + ' index ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam')
samfile = pysam.AlignmentFile(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.bam', 'rb')

for read in samfile.fetch():
    blocks = re.findall('(\d+)(\D+)', read.cigarstring)
    cumlen = np.cumsum([0] + [int(item[0]) if item[1] in {'M', 'D', 'N'} else 0 for item in blocks][:-1])
    blocks = [(blocks[i][0], blocks[i][1], cumlen[i] + int(read.reference_start) + 1) for i in range(len(blocks))]
    junctions = [(cumlen[i] + int(read.reference_start) + 1, cumlen[i+1] + int(read.reference_start)) for i in range(len(blocks)) if blocks[i][1] == 'N']
    if ('3', 'D', 13005073) in blocks or any([item[0] == 13016933 for item in junctions]):
        assignments[read.query_name] = 'hap1'
    elif any([item[0] == 13016926 for item in junctions]):
        assignments[read.query_name] = 'hap2'
    else:
        assignments[read.query_name] = 'NA'

qscore = PhasingQuality(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + 
    '/all_reads.bam', assignments, target[1], target[2])
with open(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/phasing_stats.txt', 'w') as quality_file:
    _ = quality_file.write('Phasing quality score\t%s\n' % format(qscore, '.4f'))
    _ = quality_file.write('Read count (haplotype 1)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap1'}))
    _ = quality_file.write('Read count (haplotype 2)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap2'}))
    _ = quality_file.write('Read count (unassigned)\t%s\n' % len({key for key in assignments if assignments[key] == 'NA'}))

WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap1'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
    gene_name + '/hap1_reads.bam')

WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap2'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
    gene_name + '/hap2_reads.bam')

WriteBam(samtools, {key for key in assignments if assignments[key] == 'NA'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
    gene_name + '/unassigned_reads.bam')

_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.bam && rm ' + 
    workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.bam.bai')

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name)  

# =====================================================================================================================
#                                                       Q1513
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1513', 'PMD-359', 'FBXL4'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name) 

# =====================================================================================================================
#                                                      FH-2209
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'FH-2209', 'PMD-359', 'FH'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Try phasing reads based the following heterozygous variants identified from previous testing:
#   * chr1:241508631 T>C
#   * chr1:241519682 A>AG
_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/*')
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2,
    pd.DataFrame([['chr1', '241508631', '.', 'T', 'C', '.', 'PASS', '.', 'GT', '1|0', '1|0'],
    ['chr1', '241519682', '.', 'A', 'AG', '.', 'PASS', '.', 'GT', '0|1', '0|1']], columns = list(range(10)) + ['DNA']))

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name) 

# =====================================================================================================================
#                                                        JaWe
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'JaWe', 'CDG-466', 'MPDU1'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name) 

# =====================================================================================================================
#                                                    FCDGC-02003
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'FCDGC-02003', 'CDG-466', 'NGLY1'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Try phasing reads based the following heterozygous variants identified from previous testing:
#   * chr3:25729207 TTTGA>T
_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/*')
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2,
    pd.DataFrame([['chr3', '25729207', '.', 'TTTGA', 'T', '.', 'PASS', '.', 'GT', '1|0', '1|0']], columns = list(range(10)) + ['DNA']))

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name) 

# =====================================================================================================================
#                                                        Q1687
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1687', 'PMD-359', 'NUBPL'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name) 

# =====================================================================================================================
#                                                        AnJa
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'AnJa', 'CDG-466', 'PIGN'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Try phasing reads based the following heterozygous variants identified from previous testing:
#   * chr18:62106985 C>G
#   * chr18:62113129 C>T
_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/*')
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2,
    pd.DataFrame([['chr18', '62106985', '.', 'C', 'G', '.', 'PASS', '.', 'GT', '1|0', '1|0'],
    ['chr18', '62113129', '.', 'C', 'T', '.', 'PASS', '.', 'GT', '0|1', '0|1']], columns = list(range(10)) + ['DNA']))

# From our initial round of read phasing, reads assigned to haplotype 1 featured complete skipping of exon 18 whereas
# reads assigned to haplotype 2 featured inclusion of exon 18. Many unphased reads featured exon 18 inclusion as well as
# exon 16 skipping. Based on this information, we will adopt the following strategy for phasing reads:
#   Haplotype 1: Reads that feature the 'A' allele at chr18:62113160 
#   Haplotype 2: Reads that feature the 'T' allele at chr18:62113160 or features a junction with donor splice site chr18:62106985
_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/*')
targetDF = pd.read_csv(workdir + '/' + gene_panel.split('-')[0] + '/references/target_genes.bed', sep = '\t', header = None)
target, assignments = targetDF[targetDF[4] == gene_name].values[0], dict()
_ = os.system(samtools + ' view -hb -F 256 -q 1 ' + workdir + '/' + gene_panel.split('-')[0] +
    '/' + patient_id + '/RNA/' + patient_id + '_TEQUILA.bam' + ' ' + target[0] + ':' + str(target[1]+1) + 
    '-' + str(target[2]) + ' > ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + 
    '/' + gene_name + '/all_reads.bam && ' + samtools + ' index ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam')
samfile = pysam.AlignmentFile(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.bam', 'rb')

for read in samfile.fetch():
    blocks = re.findall('(\d+)(\D+)', read.cigarstring)
    cumlen = np.cumsum([0] + [int(item[0]) if item[1] in {'M', 'D', 'N'} else 0 for item in blocks][:-1])
    blocks = [(blocks[i][0], blocks[i][1], cumlen[i] + int(read.reference_start) + 1) for i in range(len(blocks))]
    junctions = [(cumlen[i] + int(read.reference_start) + 1, cumlen[i+1] + int(read.reference_start)) for i in range(len(blocks)) if blocks[i][1] == 'N']
    if any([item[1] == 62106985 for item in junctions]):
        assignments[read.query_name] = 'hap2'
    else:
        assignments[read.query_name] = 'NA'

for pileupcolumn in pileup_truncated(samfile, 'chr18', 62113159, 62113160):
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            if pileupread.alignment.query_sequence[pileupread.query_position] == 'A':
                if assignments.get(pileupread.alignment.query_name, 'NA') == 'NA':
                    assignments[pileupread.alignment.query_name] = 'hap1'
            elif pileupread.alignment.query_sequence[pileupread.query_position] == 'T':
                if assignments.get(pileupread.alignment.query_name, 'NA') == 'NA':
                    assignments[pileupread.alignment.query_name] = 'hap2'

qscore = PhasingQuality(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + 
    '/all_reads.bam', assignments, target[1], target[2])
with open(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/phasing_stats.txt', 'w') as quality_file:
    _ = quality_file.write('Phasing quality score\t%s\n' % format(qscore, '.4f'))
    _ = quality_file.write('Read count (haplotype 1)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap1'}))
    _ = quality_file.write('Read count (haplotype 2)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap2'}))
    _ = quality_file.write('Read count (unassigned)\t%s\n' % len({key for key in assignments if assignments[key] == 'NA'}))

WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap1'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
    gene_name + '/hap1_reads.bam')

WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap2'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
    gene_name + '/hap2_reads.bam')

WriteBam(samtools, {key for key in assignments if assignments[key] == 'NA'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
    gene_name + '/unassigned_reads.bam')

_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.bam && rm ' + 
    workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.bam.bai')

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name) 

# =====================================================================================================================
#                                                     CDG-110-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-110-1', 'CDG-466', 'PIGN'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name) 

# =====================================================================================================================
#                                                     CDG-137-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-137-1', 'CDG-466', 'PIGN'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name) 

# =====================================================================================================================
#                                                     CDG-147-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-147-1', 'CDG-466', 'PIGN'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name)

# =====================================================================================================================
#                                                     CDG-161-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-161-1', 'CDG-466', 'PIGN'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Try phasing reads based the following heterozygous variants identified from previous testing:
#   * chr18:62113129 C>T
#   * chr18:62161218 ATCT>A
_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/*')
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2,
    pd.DataFrame([['chr18', '62113129', '.', 'C', 'T', '.', 'PASS', '.', 'GT', '1|0', '1|0'],
    ['chr18', '62161218', '.', 'ATCT', 'A', '.', 'PASS', '.', 'GT', '0|1', '0|1']], columns = list(range(10)) + ['DNA']))

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name)

# =====================================================================================================================
#                                                     CDG-183-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-183-1', 'CDG-466', 'PIGN'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Try phasing reads based the following heterozygous variants identified from previous testing:
#   * chr18:62084583 C>T
#   * chr18:62106797 G>A
_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/*')
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2,
    pd.DataFrame([['chr18', '62084583', '.', 'C', 'T', '.', 'PASS', '.', 'GT', '1|0', '1|0'],
    ['chr18', '62106797', '.', 'G', 'A', '.', 'PASS', '.', 'GT', '0|1', '0|1']], columns = list(range(10)) + ['DNA']))

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name) 

# =====================================================================================================================
#                                                     CDG-184-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-184-1', 'CDG-466', 'PIGN'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name) 

# =====================================================================================================================
#                                                     CDG-132-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-132-1', 'CDG-466', 'PIGQ'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Try phasing reads based the following heterozygous variants identified from previous testing:
#   * chr16:576255 G>A
#   * chr16:578911 TCTA>T
_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/*')
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2,
    pd.DataFrame([['chr16', '576255', '.', 'G', 'A', '.', 'PASS', '.', 'GT', '1|0', '1|0'],
    ['chr16', '578911', '.', 'TCTA', 'T', '.', 'PASS', '.', 'GT', '0|1', '0|1']], columns = list(range(10)) + ['DNA']))

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name) 

# =====================================================================================================================
#                                                       11547
# =====================================================================================================================

patient_id, gene_panel, gene_name = '11547', 'CDG-466', 'PMM2'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name) 

# =====================================================================================================================
#                                                     CDG-125-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-125-1', 'CDG-466', 'PMM2'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name) 

# =====================================================================================================================
#                                                     CDG-152-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-152-1', 'CDG-466', 'PMM2'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Try phasing reads based the following heterozygous variants identified from previous testing:
#   * chr16:8806398 C>T
_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/*')
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2,
    pd.DataFrame([['chr16', '8806398', '.', 'C', 'T', '.', 'PASS', '.', 'GT', '1|0', '1|0']], columns = list(range(10)) + ['DNA']))

# From our initial round of read phasing, reads assigned to haplotype 1 featured the missense variant chr16:8806398 C>T
# but there were little to no reads assigned to haplotype 2, which carries a structural deletion that removes exons 3-7
# Nearly all unphased reads featured skipping of exons 3-7, which is consistent with haplotype 2. Based on this information, 
# we will adopt the following strategy for phasing reads:
#   Haplotype 1: Reads that feature any exon in the region chr16:8805603-8820989
#   Haplotype 2: Reads that feature the junction skipping exons 3 to 7
_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/*')
targetDF = pd.read_csv(workdir + '/' + gene_panel.split('-')[0] + '/references/target_genes.bed', sep = '\t', header = None)
target, assignments = targetDF[targetDF[4] == gene_name].values[0], dict()
_ = os.system(samtools + ' view -hb -F 256 -q 1 ' + workdir + '/' + gene_panel.split('-')[0] +
    '/' + patient_id + '/RNA/' + patient_id + '_TEQUILA.bam' + ' ' + target[0] + ':' + str(target[1]+1) + 
    '-' + str(target[2]) + ' > ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + 
    '/' + gene_name + '/all_reads.bam && ' + samtools + ' index ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam')
samfile = pysam.AlignmentFile(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.bam', 'rb')

for read in samfile.fetch():
    blocks = re.findall('(\d+)(\D+)', read.cigarstring)
    cumlen = np.cumsum([0] + [int(item[0]) if item[1] in {'M', 'D', 'N'} else 0 for item in blocks][:-1])
    blocks = [(blocks[i][0], blocks[i][1], cumlen[i] + int(read.reference_start) + 1) for i in range(len(blocks))]
    junctions = [(cumlen[i] + int(read.reference_start) + 1, cumlen[i+1] + int(read.reference_start)) for i in range(len(blocks)) if blocks[i][1] == 'N']
    if any([item[1] in {'M', 'D', 'I'} and item[2] > 8805603 and item[2] < 8820989 for item in blocks]):
        assignments[read.query_name] = 'hap1'
    elif (8801911, 8847723) in junctions:
        assignments[read.query_name] = 'hap2'
    else:
        assignments[read.query_name] = 'NA'

qscore = PhasingQuality(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + 
    '/all_reads.bam', assignments, target[1], target[2])
with open(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/phasing_stats.txt', 'w') as quality_file:
    _ = quality_file.write('Phasing quality score\t%s\n' % format(qscore, '.4f'))
    _ = quality_file.write('Read count (haplotype 1)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap1'}))
    _ = quality_file.write('Read count (haplotype 2)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap2'}))
    _ = quality_file.write('Read count (unassigned)\t%s\n' % len({key for key in assignments if assignments[key] == 'NA'}))

WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap1'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
    gene_name + '/hap1_reads.bam')

WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap2'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
    gene_name + '/hap2_reads.bam')

WriteBam(samtools, {key for key in assignments if assignments[key] == 'NA'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
    gene_name + '/unassigned_reads.bam')

_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.bam && rm ' + 
    workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.bam.bai')

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name)

# =====================================================================================================================
#                                                     CDG-115-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-115-1', 'CDG-466', 'SLC35A2'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name)

# =====================================================================================================================
#                                                       08013
# =====================================================================================================================

patient_id, gene_panel, gene_name = '08013', 'CDG-466', 'SLC35C1'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Try phasing reads based the following heterozygous variants identified from previous testing:
#   * chr11:45806301 CCTT>C
_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/*')
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2,
    pd.DataFrame([['chr11', '45806301', '.', 'CCTT', 'C', '.', 'PASS', '.', 'GT', '1|0', '1|0']], columns = list(range(10)) + ['DNA']))

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name)

# =====================================================================================================================
#                                                       Q1819
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1819', 'PMD-359', 'TRIT1'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Try phasing reads based the following heterozygous variants identified from previous testing:
#   * chr1:39850120 T>C
#   * chr1:39883470 G>A
_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/*')
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2,
    pd.DataFrame([['chr1', '39850120', '.', 'T', 'C', '.', 'PASS', '.', 'GT', '1|0', '1|0'],
    ['chr1', '39883470', '.', 'G', 'A', '.', 'PASS', '.', 'GT', '0|1', '0|1']], columns = list(range(10)) + ['DNA']))

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name)

# =====================================================================================================================
#                                                       Q1642
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1642', 'PMD-359', 'TRMU'

# Try phasing reads without using any patient genotyping information
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2, None)

# Try phasing reads based the following heterozygous variants identified from previous testing:
#   * chr22:46356055 G>A
_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/*')
PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2,
    pd.DataFrame([['chr22', '46356055', '.', 'G', 'A', '.', 'PASS', '.', 'GT', '1|0', '1|0']], columns = list(range(10)) + ['DNA']))

# From our initial round of read phasing, reads assigned to haplotype 1 featured the missense variant chr22:46356055 G>A
# but there were little to no reads assigned to haplotype 2, which carries a structural deletion that removes exon 1
# A handful of unphased reads featured inclusion of exon 1, which is consistent with haplotype 1. Based on this information, 
# we will adopt the following strategy for phasing reads:
#   Haplotype 1: Reads that feature any exon in the region chr22:46334556-46336380 OR feature the 'A' allele at chr22:46356055
#   Haplotype 2: Reads that feature the 'G' allele at chr22:46356055
_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/*')
targetDF = pd.read_csv(workdir + '/' + gene_panel.split('-')[0] + '/references/target_genes.bed', sep = '\t', header = None)
target, assignments = targetDF[targetDF[4] == gene_name].values[0], dict()
_ = os.system(samtools + ' view -hb -F 256 -q 1 ' + workdir + '/' + gene_panel.split('-')[0] +
    '/' + patient_id + '/RNA/' + patient_id + '_TEQUILA.bam' + ' ' + target[0] + ':' + str(target[1]+1) + 
    '-' + str(target[2]) + ' > ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + 
    '/' + gene_name + '/all_reads.bam && ' + samtools + ' index ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam')
samfile = pysam.AlignmentFile(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.bam', 'rb')

for read in samfile.fetch():
    blocks = re.findall('(\d+)(\D+)', read.cigarstring)
    cumlen = np.cumsum([0] + [int(item[0]) if item[1] in {'M', 'D', 'N'} else 0 for item in blocks][:-1])
    blocks = [(blocks[i][0], blocks[i][1], cumlen[i] + int(read.reference_start) + 1) for i in range(len(blocks))]
    junctions = [(cumlen[i] + int(read.reference_start) + 1, cumlen[i+1] + int(read.reference_start)) for i in range(len(blocks)) if blocks[i][1] == 'N']
    if any([item[1] in {'M', 'D', 'I'} and item[2] > 46334556 and item[2] < 46336380 for item in blocks]):
        assignments[read.query_name] = 'hap1'
    else:
        assignments[read.query_name] = 'NA'

for pileupcolumn in pileup_truncated(samfile, 'chr18', 46356054, 46356055):
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            if pileupread.alignment.query_sequence[pileupread.query_position] == 'A':
                if assignments.get(pileupread.alignment.query_name, 'NA') == 'NA':
                    assignments[pileupread.alignment.query_name] = 'hap1'
            elif pileupread.alignment.query_sequence[pileupread.query_position] == 'G':
                if assignments.get(pileupread.alignment.query_name, 'NA') == 'NA':
                    assignments[pileupread.alignment.query_name] = 'hap2'

qscore = PhasingQuality(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + 
    '/all_reads.bam', assignments, target[1], target[2])
with open(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/phasing_stats.txt', 'w') as quality_file:
    _ = quality_file.write('Phasing quality score\t%s\n' % format(qscore, '.4f'))
    _ = quality_file.write('Read count (haplotype 1)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap1'}))
    _ = quality_file.write('Read count (haplotype 2)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap2'}))
    _ = quality_file.write('Read count (unassigned)\t%s\n' % len({key for key in assignments if assignments[key] == 'NA'}))

WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap1'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
    gene_name + '/hap1_reads.bam')

WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap2'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
    gene_name + '/hap2_reads.bam')

WriteBam(samtools, {key for key in assignments if assignments[key] == 'NA'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + 
    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + 
    gene_name + '/unassigned_reads.bam')

_ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.bam && rm ' + 
    workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/' + patient_id + '/' + gene_name + '/all_reads.bam.bai')

# Compute p-value for gene expression dosage outlier test
DosageOutlierTest(workdir, patient_id, gene_panel, gene_name)
