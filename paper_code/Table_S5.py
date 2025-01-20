#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2024.12.29
Supplementary Table 5

Detection and phasing of previously identified pathogenic variants within exonic regions from patient TEQUILA-seq reads. 
Genotyping information from patient DNA sequencing data was not used to guide read phasing or variant detection.
'''

# =====================================================================================================================
#                                                   HELPER FUNCTIONS
# =====================================================================================================================

def PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2):
    _ = os.system('mkdir -p ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id)
    if os.path.isdir(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name):
        _ = os.system('rm -rf ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name)
    _ = os.system('mkdir -p ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name)
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
    merge_variants = MergeVariants(rna_variants, pd.DataFrame(columns = list(range(10)) + ['DNA']))
    _ = os.system(samtools + ' view -hb -F 256 -q 1 ' + workdir + '/' + gene_panel.split('-')[0] +
        '/' + patient_id + '/RNA/' + patient_id + '_TEQUILA.bam' + ' ' + target[0] + ':' + str(target[1]+1) + 
        '-' + str(target[2]) + ' > ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + 
        '/' + gene_name + '/all_reads.bam && ' + samtools + ' index ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + 
        patient_id + '/' + gene_name + '/all_reads.bam')
    if merge_variants.shape[0] > 0:
        samfile = pysam.AlignmentFile(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + 
            '/' + gene_name + '/all_reads.bam', mode = 'rb')
        merge_variants['Coverage'] = merge_variants.apply(lambda x: np.sum(samfile.count_coverage(x[0], int(x[1])-1, int(x[1]))), axis = 1)
        merge_variants = merge_variants[merge_variants['Coverage'] >= 20]
        if merge_variants.shape[0] > 0:
            bestPS = merge_variants[['PS', 'Coverage']].groupby('PS').max().reset_index().sort_values(by = 'Coverage', ascending = False).head(1)['PS'].item()
            merge_variants = merge_variants[merge_variants['PS'] == bestPS][[0, 1, 3, 4, 'RNA']]
            MakeHaplotypeFasta(samtools, bcftools, merge_variants, genome, target[0] + ':' + str(target[1]+1) + '-' + str(target[2]), 
                workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name)
            _ = os.system('(' + samtools + ' fastq -n ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + 
                '/' + gene_name + '/all_reads.bam | gzip > ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + 
                '/' + gene_name + '/all_reads.fastq.gz) > /dev/null 2>&1')
            _ = os.system('(' + minimap2 + ' -ax splice --secondary=no ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + 
                patient_id + '/' + gene_name + '/gene.fa ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + 
                '/' + gene_name + '/all_reads.fastq.gz > ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + 
                '/' + gene_name + '/all_reads.sam) > /dev/null 2>&1')
            _ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + 
                '/all_reads.fastq.gz && rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + 
                gene_name + '/gene.fa && rm ' +  workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + 
                gene_name + '/gene.fa.fai')
            assignments = HaplotypeAssignment(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/all_reads.sam')
            _ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/all_reads.sam')
            if len({key for key in assignments if assignments[key] != 'NA'}) >= 100:
                qscore = PhasingQuality(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + 
                    '/all_reads.bam', assignments, target[1], target[2])
                with open(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/phasing_stats.txt', 'w') as quality_file:
                    _ = quality_file.write('Phasing quality score\t%s\n' % format(qscore, '.4f'))
                    _ = quality_file.write('Read count (haplotype 1)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap1'}))
                    _ = quality_file.write('Read count (haplotype 2)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap2'}))
                    _ = quality_file.write('Read count (unassigned)\t%s\n' % len({key for key in assignments if assignments[key] == 'NA'}))
                WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap1'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + 
                    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + 
                    gene_name + '/hap1_reads.bam')
                WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap2'}, workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + 
                    patient_id + '/' + gene_name + '/all_reads.bam', workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + 
                    gene_name + '/hap2_reads.bam')
                _ = os.system('(' + bcftools + ' mpileup -Ou -f ' + genome + ' -r ' + target[0] + ':' + str(target[1]+1) + '-' + 
                    str(target[2]) + ' -B -Q 5 ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + 
                    '/hap1_reads.bam ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + 
                    '/hap2_reads.bam | ' + bcftools + ' call -mv -Oz ' + '--ploidy 1 -o ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + 
                    patient_id + '/' + gene_name + '/bcftools.vcf.gz) > /dev/null 2>&1')
            else:
                _ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/gene.vcf.gz && rm ' + 
                    workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/gene.vcf.gz.tbi')
    _ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/all_reads.bam && rm ' + 
        workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/all_reads.bam.bai')

def CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome):
    chrom, pos, ref, alt = var.replace(':', '_').replace(' ', '_').replace('>', '_').split('_')
    targetDF = pd.read_csv(workdir + '/' + gene_panel.split('-')[0] + '/references/target_genes.bed', sep = '\t', header = None)
    target = targetDF[targetDF[4] == gene_name].values[0]
    command = bcftools + ' view -H -r ' + target[0] + ':' + str(target[1]+1) + '-' + str(target[2])
    longcallr_variants = pd.DataFrame([item.split() for item in subprocess.run(command + ' ' + workdir + '/' + gene_panel.split('-')[0] +
        '/' + patient_id + '/RNA/stripe/target_genes/longcallR/output.vcf.gz', shell = True, capture_output = True, text = True).stdout.splitlines()], 
        columns = list(range(10)))
    longcallr_variants = longcallr_variants[(longcallr_variants[3] != longcallr_variants[4])]
    longcallr_variants['Dense'] = CheckDensity(sorted(longcallr_variants[1].astype(int).tolist()))
    longcallr_variants = longcallr_variants[(longcallr_variants[1] == pos) & (longcallr_variants[3] == ref) & (longcallr_variants[4] == alt)]
    longcallr_variants['DP'] = [PullFeature(item, 'DP') for item in longcallr_variants[[8, 9]].values.tolist()]
    longcallr_variants['AF'] = [PullFeature(item, 'AF').split(',')[0] for item in longcallr_variants[[8, 9]].values.tolist()]
    longcallr_variants['GT'] = [PullFeature(item, 'GT').split(',')[0] for item in longcallr_variants[[8, 9]].values.tolist()]
    longcallr_variants = longcallr_variants[[0, 1, 3, 4, 6, 'Dense', 'DP', 'AF', 'GT']]
    longcallr_variants.columns = ['CHROM', 'POS', 'REF', 'ALT', 'LONGCALLR_FILTER', 'LONGCALLR_DENSE', 'LONGCALLR_DP', 'LONGCALLR_AF', 'LONGCALLR_GT']
    longcallr_variants.loc[:, 'LONGCALLR_AF'] = longcallr_variants['LONGCALLR_AF'].apply(lambda x: '{:,.2f}'.format(float(x)))
    clair3_variants = pd.DataFrame([item.split() for item in subprocess.run(command + ' ' + workdir + '/' + gene_panel.split('-')[0] +
        '/' + patient_id + '/RNA/stripe/target_genes/clair3/pileup.vcf.gz', shell = True, capture_output = True, text = True).stdout.splitlines()], 
        columns = list(range(10)))
    clair3_variants = clair3_variants[(clair3_variants[3] != clair3_variants[4]) & clair3_variants[6].isin({'PASS', 'LowQual'})]
    clair3_variants['Dense'] = CheckDensity(sorted(clair3_variants[1].astype(int).tolist()))
    clair3_variants = clair3_variants[(clair3_variants[1] == pos) & (clair3_variants[3] == ref) & (clair3_variants[4] == alt)]
    clair3_variants['DP'] = [PullFeature(item, 'DP') for item in clair3_variants[[8, 9]].values.tolist()]
    clair3_variants['AF'] = [PullFeature(item, 'AF').split(',')[0] for item in clair3_variants[[8, 9]].values.tolist()]
    clair3_variants['GT'] = [PullFeature(item, 'GT').split(',')[0] for item in clair3_variants[[8, 9]].values.tolist()]
    clair3_variants = clair3_variants[[0, 1, 3, 4, 6, 'Dense', 'DP', 'AF', 'GT']]
    clair3_variants.columns = ['CHROM', 'POS', 'REF', 'ALT', 'CLAIR3_FILTER', 'CLAIR3_DENSE', 'CLAIR3_DP', 'CLAIR3_AF', 'CLAIR3_GT']
    clair3_variants.loc[:, 'CLAIR3_AF'] = clair3_variants['CLAIR3_AF'].apply(lambda x: '{:,.2f}'.format(float(x)))
    dv_variants = pd.DataFrame([item.split() for item in subprocess.run(command + ' ' + workdir + '/' + gene_panel.split('-')[0] +
        '/' + patient_id + '/RNA/stripe/target_genes/deepvariant/PROBAND.vcf.gz', shell = True, capture_output = True, text = True).stdout.splitlines()], 
        columns = list(range(10)))
    dv_variants = dv_variants[(dv_variants[3] != dv_variants[4]) & dv_variants[6].isin({'PASS', 'LowQual'})]
    dv_variants['Dense'] = CheckDensity(sorted(dv_variants[1].astype(int).tolist()))
    dv_variants = dv_variants[(dv_variants[1] == pos) & (dv_variants[3] == ref) & (dv_variants[4] == alt)]
    dv_variants['DP'] = [PullFeature(item, 'DP') for item in dv_variants[[8, 9]].values.tolist()]
    dv_variants['VAF'] = [PullFeature(item, 'VAF').split(',')[0] for item in dv_variants[[8, 9]].values.tolist()]
    dv_variants['GT'] = [PullFeature(item, 'GT').split(',')[0] for item in dv_variants[[8, 9]].values.tolist()]
    dv_variants = dv_variants[[0, 1, 3, 4, 6, 'Dense', 'DP', 'VAF', 'GT']]
    dv_variants.columns = ['CHROM', 'POS', 'REF', 'ALT', 'DEEPVARIANT_FILTER', 'DEEPVARIANT_DENSE', 'DEEPVARIANT_DP', 'DEEPVARIANT_AF', 'DEEPVARIANT_GT']
    dv_variants.loc[:, 'DEEPVARIANT_AF'] = dv_variants['DEEPVARIANT_AF'].apply(lambda x: '{:,.2f}'.format(float(x)))
    callset = pd.DataFrame([[chrom, pos, ref, alt]], columns = ['CHROM', 'POS', 'REF', 'ALT']).merge(longcallr_variants,
        on = ['CHROM', 'POS', 'REF', 'ALT'], how = 'left').merge(clair3_variants, on = ['CHROM', 'POS', 'REF', 'ALT'], how = 'left').merge(
        dv_variants, on = ['CHROM', 'POS', 'REF', 'ALT'], how = 'left').fillna('.')
    if os.path.isfile(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/phasing_stats.txt'):
        phasing_variants = pd.DataFrame([item.split() for item in subprocess.run(bcftools + ' view -H ' + workdir + 
            '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/bcftools.vcf.gz', 
            shell = True, capture_output = True, text = True).stdout.splitlines()], columns = list(range(11)))
        phasing_variants = phasing_variants[phasing_variants[3] != phasing_variants[4]]
        phasing_variants['Dense'] = CheckDensity(sorted(phasing_variants[1].astype(int).tolist()))
        phasing_variants = phasing_variants[(phasing_variants[1] == pos) & (phasing_variants[3] == ref) & (phasing_variants[4] == alt)]
        phasing_variants['REF1'] = [PullFeature(item, 'AD').split(',')[0] for item in phasing_variants[[8, 9]].values.tolist()]
        phasing_variants['ALT1'] = [PullFeature(item, 'AD').split(',')[1] for item in phasing_variants[[8, 9]].values.tolist()]
        phasing_variants['GT1'] = [PullFeature(item, 'GT').split(',')[0] for item in phasing_variants[[8, 9]].values.tolist()]
        phasing_variants['REF2'] = [PullFeature(item, 'AD').split(',')[0] for item in phasing_variants[[8, 10]].values.tolist()]
        phasing_variants['ALT2'] = [PullFeature(item, 'AD').split(',')[1] for item in phasing_variants[[8, 10]].values.tolist()]
        phasing_variants['GT2'] = [PullFeature(item, 'GT').split(',')[0] for item in phasing_variants[[8, 10]].values.tolist()]
        phasing_variants['DP1'] = phasing_variants['REF1'].astype(int) + phasing_variants['ALT1'].astype(int)
        phasing_variants['DP2'] = phasing_variants['REF2'].astype(int) + phasing_variants['ALT2'].astype(int)
        phasing_variants['DP'] = (phasing_variants['DP1'] + phasing_variants['DP2']).astype(object)
        phasing_variants['AF'] = ((phasing_variants['ALT1'].astype(int) + phasing_variants['ALT2'].astype(int))/phasing_variants['DP']).apply(lambda x: '{:,.2f}'.format(float(x)))
        if phasing_variants.shape[0] > 0:
            phasing_variants['GT'] = phasing_variants['GT1'].replace('.', '0') + '|' + phasing_variants['GT2'].replace('.', '0')
            phasing_variants = phasing_variants[[0, 1, 3, 4, 6, 'Dense', 'DP', 'AF', 'GT']]
            phasing_variants.columns = ['CHROM', 'POS', 'REF', 'ALT', 'PHASING_FILTER', 'PHASING_DENSE', 'PHASING_DP', 'PHASING_AF', 'PHASING_GT']
        else:
            phasing_variants = pd.DataFrame(columns = ['CHROM', 'POS', 'REF', 'ALT', 'PHASING_FILTER', 'PHASING_DENSE', 'PHASING_DP', 'PHASING_AF', 'PHASING_GT'])
        callset = callset.merge(phasing_variants, on = ['CHROM', 'POS', 'REF', 'ALT'], how = 'left').fillna('.')
    else:
        callset = callset.merge(pd.DataFrame(columns = ['CHROM', 'POS', 'REF', 'ALT', 'PHASING_FILTER', 'PHASING_DENSE', 'PHASING_DP', 'PHASING_AF', 'PHASING_GT']), 
            on = ['CHROM', 'POS', 'REF', 'ALT'], how = 'left').fillna('.')
    _ = os.system('mkdir -p ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/CADD')
    pd.DataFrame([[item[0].replace('chr', ''), item[1], '.', item[2], item[3]] for item in callset.values]).to_csv(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + 
        patient_id + '/' + gene_name + '/CADD/input.vcf', sep = '\t', header = False, index = False)
    _ = os.system('(' + cadd + ' -o ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/CADD/output.tsv.gz ' +
        workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/CADD/input.vcf) > /dev/null 2>&1')
    if os.path.isfile(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/CADD/output.tsv.gz'):
        cadd_results = pd.read_csv(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/CADD/output.tsv.gz', header = None, 
            sep = '\t', compression = 'gzip', names = ['Chrom', 'Pos', 'Ref', 'Alt', 'RawScore', 'PHRED'], comment = '#')
        callset['CADD'] = ((callset['CHROM'] + '_' + callset['POS'].astype(str) + '_' + callset['REF'] + '_' + callset['ALT']).map(dict(zip(
            'chr' + cadd_results['Chrom'].astype(str) + '_' + cadd_results['Pos'].astype(str) + '_' + cadd_results['Ref'] + '_' + cadd_results['Alt'], 
            cadd_results['PHRED'])))).apply(lambda x: '{:,.1f}'.format(float(x)))
        _ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/CADD/input.vcf && rm ' + 
            workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/CADD/output.tsv.gz && rmdir ' + 
            workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/CADD')
    else:
        callset['CADD'] = 0
        _ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/CADD/input.vcf && rmdir ' + 
             workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/CADD')
    _ = os.system('mkdir -p ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/SpliceAI')
    _ = os.system('echo "##fileformat=VCFv4.2" > ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/SpliceAI/input.vcf')
    _ = os.system(bcftools + ' view -h ' + workdir + '/' + gene_panel.split('-')[0] + '/' + patient_id + '/RNA/stripe/target_genes/longcallR/output.vcf.gz | grep "##contig=" >> ' + 
        workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/SpliceAI/input.vcf')
    pd.DataFrame([[item[0], item[1], '.', item[2], item[3], '.', '.', '.'] for item in callset.values], columns = ['#CHROM', 'POS', 'ID',
        'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']).to_csv(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/SpliceAI/input.vcf', 
        mode = 'a', sep = '\t', index = False)
    _ = os.system('(spliceai -I ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/SpliceAI/input.vcf -O ' +
        workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/SpliceAI/output.vcf -R ' + genome + ' -A grch38) > /dev/null 2>&1')
    if os.path.isfile(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/SpliceAI/output.vcf'):
        spliceai_results = pd.read_csv(workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/SpliceAI/output.vcf', header = None, sep = '\t',
            names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'], comment = '#')
        spliceai_results['SPLICEAI'] = spliceai_results['INFO'].apply(lambda x: np.max([float(val) if val != '.' else 0 for val in x.split('|')[2:6]]) if 'SpliceAI' in x else np.nan)
        callset['SPLICEAI'] = ((callset['CHROM'] + '_' + callset['POS'].astype(str) + '_' + callset['REF'] + '_' + callset['ALT']).map(dict(zip(
            spliceai_results['CHROM'] + '_' + spliceai_results['POS'].astype(str) + '_' + spliceai_results['REF'] + '_' + spliceai_results['ALT'],
            spliceai_results['SPLICEAI'])))).apply(lambda x: '{:,.2f}'.format(float(x)))
        _ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/SpliceAI/input.vcf && rm ' + workdir + 
            '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/SpliceAI/output.vcf && rmdir ' + workdir + 
            '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/SpliceAI')
    else:
        callset['SPLICEAI'] = 0
        _ = os.system('rm ' + workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/SpliceAI/input.vcf && rmdir ' + workdir + 
            '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/' + patient_id + '/' + gene_name + '/SpliceAI')
    callset['GNOMAD_AF_GRPMAX_JOINT'], callset['GNOMAD_NHOMALT_GRPMAX_JOINT'] = zip(*AnnotateGnomad((callset['CHROM'] + '_' + callset['POS'].astype(str) + 
        '_' + callset['REF'] + '_' + callset['ALT']).values, workdir + '/' + gene_panel.split('-')[0] + '/references/gnomAD_v4/' + target[4] + '/gnomad_variants.tsv'))
    callset.loc[:, 'GNOMAD_AF_GRPMAX_JOINT'] = callset['GNOMAD_AF_GRPMAX_JOINT'].apply(lambda x: '{:,.3e}'.format(x))
    callset = callset.drop('GNOMAD_NHOMALT_GRPMAX_JOINT', axis = 1)
    callset['CLINVAR'] = AnnotateClinvar(callset[['CHROM', 'POS', 'REF', 'ALT']].values, clinvar)
    callset['HOMOPOLYMER'] = CheckHomopolymer(callset[['CHROM', 'POS']].values, genome)
    callset['END_DISTANCE_PCT'], callset['SKIP_FREQ'] = zip(*CalculateAlignmentStats(callset[['CHROM', 'POS', 'REF', 'ALT']].values, 
        workdir + '/' + gene_panel.split('-')[0] + '/' + patient_id + '/RNA/' + patient_id + '_TEQUILA.bam'))
    callset.loc[:, 'END_DISTANCE_PCT'] = callset['END_DISTANCE_PCT'].apply(lambda x: '{:,.3f}'.format(x))
    callset.loc[:, 'SKIP_FREQ'] = callset['SKIP_FREQ'].apply(lambda x: '{:,.3f}'.format(x))
    return list(callset.values[0])[4:]
        
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
cadd = '/scr1/users/wangr5/tools/CADD-scripts-1.6.post1/CADD.sh'
clinvar = '/scr1/users/wangr5/databases/ClinVar/clinvar_20240917.vcf.gz'
outfile = workdir + '/manuscript/Revisions/20241219/Supplementary_Tables/Table_S5/Table_S5.txt'

# Import functions used in Run_STRIPE.py
sys.path.append(workdir + '/scripts/main')
from Run_STRIPE import *

# Create an empty data frame with columns of interest
outDF = pd.DataFrame(columns = ['patient_id', 'gene_name', 'variant_grch38', 'variant_class', 'genotype'] + 
    [tool + '_' + field for tool in ['longcallR', 'clair3', 'deepvariant', 'phasing'] for field in 
    ['filter', 'dense', 'dp', 'af', 'gt']] + ['cadd', 'spliceai', 'gnomad_af_grpmax_joint', 'clinvar', 
    'homopolymer', 'end_distance_pct', 'skipping_frequency'])

# =====================================================================================================================
#                                                     CDG-114-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-114-1', 'CDG-466', 'ALG3'
variants, classes, genotypes = ['chr3:184243927 G>A'], ['Missense'], ['1/1']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                       Q1513
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1513', 'PMD-359', 'FBXL4'
variants, classes, genotypes = ['chr6:98874306 A>T'], ['Missense'], ['1|0']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                      FH-2209
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'FH-2209', 'PMD-359', 'FH'
variants, classes, genotypes = ['chr1:241519682 A>AG'], ['Frameshift duplication'], ['0|1']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                        JaWe
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'JaWe', 'CDG-466', 'MPDU1'
variants, classes, genotypes = ['chr17:7583881 G>T'], ['Nonsense'], ['1|0']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                    FCDGC-02003
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'FCDGC-02003', 'CDG-466', 'NGLY1'
variants, classes, genotypes = ['chr3:25729207 TTTGA>T'], ['Frameshift deletion'], ['1|0']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                     CDG-110-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-110-1', 'CDG-466', 'PIGN'
variants, classes, genotypes = ['chr18:62095902 C>T'], ['Missense'], ['1|0']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                     CDG-137-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-137-1', 'CDG-466', 'PIGN'
variants, classes, genotypes = ['chr18:62143337 A>C'], ['Missense'], ['0|1']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                     CDG-147-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-147-1', 'CDG-466', 'PIGN'
variants, classes, genotypes = ['chr18:62113277 C>T', 'chr18:62143306 C>T'], ['Missense', 'Synonymous'], ['1|0', '0|1']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                     CDG-183-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-183-1', 'CDG-466', 'PIGN'
variants, classes, genotypes = ['chr18:62106797 G>A'], ['Nonsense'], ['0|1']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                     CDG-184-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-184-1', 'CDG-466', 'PIGN'
variants, classes, genotypes = ['chr18:62045973 G>C', 'chr18:62143337 A>C'], ['Missense', 'Missense'], ['1|0', '0|1']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                     CDG-132-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-132-1', 'CDG-466', 'PIGQ'
variants, classes, genotypes = ['chr16:578911 TCTA>T'], ['Inframe deletion'], ['0|1']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                       11547
# =====================================================================================================================

patient_id, gene_panel, gene_name = '11547', 'CDG-466', 'PMM2'
variants, classes, genotypes = ['chr16:8811146 G>A', 'chr16:8811153 G>A'], ['Missense', 'Missense'], ['1|0', '0|1']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                     CDG-125-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-125-1', 'CDG-466', 'PMM2'
variants, classes, genotypes = ['chr16:8811088 C>A', 'chr16:8811153 G>A'], ['Missense', 'Missense'], ['1|0', '0|1']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                     CDG-152-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-152-1', 'CDG-466', 'PMM2'
variants, classes, genotypes = ['chr16:8806398 C>T'], ['Missense'], ['0|1']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                     CDG-115-1
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'CDG-115-1', 'CDG-466', 'SLC35A2'
variants, classes, genotypes = ['chrX:48909877 C>T'], ['Missense'], ['1|0']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                       08013
# =====================================================================================================================

patient_id, gene_panel, gene_name = '08013', 'CDG-466', 'SLC35C1'
variants, classes, genotypes = ['chr11:45806301 CCTT>C'], ['Inframe deletion'], ['1|0']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                       Q1819
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1819', 'PMD-359', 'TRIT1'
variants, classes, genotypes = ['chr1:39883470 G>A'], ['Nonsense'], ['0|1']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# =====================================================================================================================
#                                                       Q1642
# =====================================================================================================================

patient_id, gene_panel, gene_name = 'Q1642', 'PMD-359', 'TRMU'
variants, classes, genotypes = ['chr22:46356055 G>A'], ['Missense'], ['0|1']

PhaseGeneReads(workdir, patient_id, gene_panel, gene_name, bcftools, samtools, genome, minimap2)
for idx, var in enumerate(variants):
    results = CheckVariant(var, workdir, patient_id, gene_panel, gene_name, bcftools, samtools, cadd, clinvar, genome)
    outDF = pd.concat([outDF, pd.DataFrame([[patient_id, gene_name, var, classes[idx], genotypes[idx]] + results], 
        columns = outDF.columns)])

# Save outDF to outfile
outDF.to_csv(outfile, sep = '\t', index = False)
