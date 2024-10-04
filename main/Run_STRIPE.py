#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2024.09.30
Version: 1.0.0

A script designed to run the computational pipeline of STRIPE. Note that this script assumes that R is installed 
and available in your $PATH.
'''

# =====================================================================================================================
#                                                       PACKAGES 
# =====================================================================================================================

# Load required packages
import sys, argparse, os, re
import numpy as np
import pandas as pd
import pysam
import subprocess
import tabix
from pysam.libcalignmentfile import IteratorColumnRegion
from scipy.optimize import brentq
from scipy.special import digamma
from scipy.stats import betabinom, beta
pd.options.mode.chained_assignment = None

# =====================================================================================================================
#                                                   HELPER FUNCTIONS
# =====================================================================================================================

def pileup_truncated(bam, contig, start, stop):
    '''
    Function designed to only return read pileups at a specific genomic interval
    '''
    has_coord, rtid, rstart, rstop = bam.parse_region(contig, start, stop)
    yield from IteratorColumnRegion(bam, tid = rtid, start = rstart, stop = rstop, truncate = True)

def ParseCigarJunctions(cigarstring, chrom, start):
    '''
    Function designed to parse a CIGAR string and retrieve coordinates of splice junctions
    '''
    blocks = re.findall('(\d+)(\D+)', cigarstring)
    cumlen = np.cumsum([0] + [int(item[0]) if item[1] in {'M', 'D', 'N'} else 0 for item in blocks][:-1])

    return ['_'.join([chrom, str(cumlen[i] + int(start) + 1), str(cumlen[i+1] + int(start))])
        for i in range(len(blocks)) if blocks[i][1] == 'N']

def ExtractSpliceJunctions(infile, samtools, tmpfile, targets = None, region = None):
    '''
    Function designed to extract splice junctions from primary alignments (MAPQ >= 1) in infile
    '''
    if targets is None and region is None:
        _ = os.system(samtools + ' view -hb -F 256 -q 1 ' + infile + ' > ' + tmpfile + 
            ' && ' + samtools + ' index ' + tmpfile)
    elif targets is None:
        _ = os.system(samtools + ' view -hb -F 256 -q 1 ' + infile + ' ' + region + ' > ' + tmpfile +
            ' && ' + samtools + ' index ' + tmpfile)
    elif region is None:
        _ = os.system(samtools + ' view -hb -L ' + targets + ' -F 256 -q 1 ' + infile + 
            ' > ' + tmpfile + ' && ' + samtools + ' index ' + tmpfile)
    
    junctions = dict()
    for read in pysam.AlignmentFile(tmpfile, 'rb').fetch():
        for item in ParseCigarJunctions(read.cigarstring, read.reference_name, read.reference_start):
            junctions[item] = junctions.get(item, 0) + 1
    
    _ = os.system('rm ' + tmpfile + ' && rm ' + tmpfile + '.bai')
    junctions = pd.DataFrame(junctions.items(), columns = ['Junction', 'Read_Counts'])
    if junctions.shape[0] > 0:
        junctions[['Chrom', 'Start', 'End']] = junctions['Junction'].str.split('_', n = 2, expand = True)
        return junctions.sort_values(by = ['Chrom', 'Start', 'End'], ascending = True).reset_index(drop = True)[['Chrom', 'Start', 'End', 'Read_Counts']]
    else:
        return pd.DataFrame(columns = ['Chrom', 'Start', 'End', 'Read_Counts'])

def PullFeature(item, featureName):
    '''
    Function designed to pull out the value of featureName in item
    '''
    return item[1].split(':')[item[0].split(':').index(featureName)] if featureName in item[0].split(':') else np.nan

def CountAlleles(genotype):
    '''
    Function designed to count the number of alleles represented in a genotype string
    '''
    return sum(list(map(int, re.split(r'\D+', genotype))))

def PhaseByTransmission(proband, parent1, parent2):
    '''
    Function designed to phase variants in proband by transmission (if unambiguous)
    '''    
    trioDF = proband[[0, 1, 3, 4, 'DNA']].merge(parent1[[0, 1, 3, 4, 'P1']], on = [0, 1, 3, 4], how = 'left'
        ).merge(parent2[[0, 1, 3, 4, 'P2']], on = [0, 1, 3, 4], how = 'left').sort_values(by = [0, 1]).fillna('0/0')
    trioDF['DNA'] = trioDF.apply(lambda x: '1|0' if x['DNA'] == '0/1' and CountAlleles(x['P1']) > CountAlleles(x['P2']) else '0|1' 
        if x['DNA'] == '0/1' and CountAlleles(x['P1']) < CountAlleles(x['P2']) else x['DNA'], axis = 1)
    return trioDF[[0, 1, 3, 4, 'DNA']]

def MergeVariants(rna_variants, dna_variants):
    '''
    Function designed to merge rna_variants and dna_variants while considering phase information
    '''
    merge_variants = rna_variants.merge(dna_variants, on = [0, 1, 3, 4], how = 'outer').fillna('.')
    merge_variants[1] = merge_variants[1].astype(int)
    merge_variants = merge_variants.sort_values(by = [0, 1])
    merge_variants[1] = merge_variants[1].astype(object)
    phase_sets = merge_variants[merge_variants['PS'] != '.']['PS'].drop_duplicates().tolist()

    # Flip phased RNA genotypes if they are completely reversed with respect to phased DNA genotypes
    for phase_set in phase_sets:
        phaseDF = merge_variants[merge_variants['PS'] == phase_set][['RNA', 'DNA']].applymap(lambda x: 1 if x == '1|0' 
            else -1 if x == '0|1' else np.nan).dropna()
        
        if phaseDF.shape[0] > 0:
            if np.all(phaseDF.product(axis = 1) == -1):
                merge_variants['RNA'] = merge_variants.apply(lambda x: '1|0' if x['PS'] == phase_set and x['RNA'] == '0|1'
                    else '0|1' if x['PS'] == phase_set and x['RNA'] == '1|0' else x['RNA'], axis = 1)
                merge_variants['PS'] = merge_variants.apply(lambda x: 'DNA' if x['PS'] == phase_set else x['PS'], axis = 1)
            elif np.all(phaseDF.product(axis = 1) == 1):
                merge_variants['PS'] = merge_variants.apply(lambda x: 'DNA' if x['PS'] == phase_set else x['PS'], axis = 1)

    # Update RNA genotypes and phase sets for remaining variants in merge_variants
    for idx, row in merge_variants.iterrows():
        if row['PS'] == '.':
            if row['DNA'] in {'1|0', '0|1'}:
                merge_variants.loc[idx, 'RNA'], merge_variants.loc[idx, 'PS'] = row['DNA'], 'DNA'
            elif row['DNA'] == '0/1' or row['RNA'] == '0/1':
                merge_variants.loc[idx, 'RNA'], merge_variants.loc[idx, 'PS'] = '1|0', '_'.join([row[0], str(row[1]), row[3], row[4]])

    return merge_variants[[0, 1, 3, 4, 'RNA', 'PS']]

def CheckDensity(positions, window_size = 200, dense_count = 4):
    '''
    Function designed to check what values in positions are dense
    '''
    density = [False] * len(positions)
    for idx in range(len(positions) - dense_count + 1):
        if positions[idx+dense_count-1] - positions[idx] <= window_size:
            density[idx:(idx+dense_count)] = [True] * dense_count
    return density

def PrintVcf(bcftools, vcfdata, outfile):
    '''
    Function designed to print vcfdata into a VCF file
    '''
    headerString = ('##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' +
        '##FORMAT=<ID=PS,Number=Integer,Type=Integer,Description="Phase Set">\n#CHROM\tPOS\tID\tREF\tALT\t' +
        'QUAL\tFILTER\tINFO\tFORMAT\tSample\n')
    with open(outfile, 'w') as vcffile:
        _ = vcffile.write(headerString)
    vcfdata.to_csv(outfile, mode = 'a', sep = '\t', index = False, header = False)
    _ = os.system(bcftools + ' view -I ' + outfile + ' -Oz -o ' + outfile + '.gz --write-index="tbi" >/dev/null 2>&1 && rm ' + outfile)

def MakeHaplotypeFasta(samtools, bcftools, variants, genome, regionString, outdir):
    '''
    Function designed to construct haplotype-specific sequences for a target gene
    '''
    PrintVcf(bcftools, pd.DataFrame([[item[0], item[1], '.', item[2], item[3], '.', 'PASS', '.', 'GT', item[4]] for item in 
        variants[variants['RNA'] == '1|0'].values]), outdir + '/hap1.vcf')
    PrintVcf(bcftools, pd.DataFrame([[item[0], item[1], '.', item[2], item[3], '.', 'PASS', '.', 'GT', item[4]] for item in 
        variants[variants['RNA'] == '0|1'].values]), outdir + '/hap2.vcf')
    PrintVcf(bcftools, pd.DataFrame([[item[0], item[1], '.', item[2], item[3], '.', 'PASS', '.', 'GT', item[4]] for item in 
        variants.values]), outdir + '/gene.vcf')
    
    _ = os.system(samtools + ' faidx ' + genome + ' ' + regionString + ' | ' + bcftools + ' consensus -s - ' +
        outdir + '/hap1.vcf.gz -o ' + outdir + '/hap1.fa >/dev/null 2>&1')
    _ = os.system(samtools + ' faidx ' + genome + ' ' + regionString + ' | ' + bcftools + ' consensus -s - ' +
        outdir + '/hap2.vcf.gz -o ' + outdir + '/hap2.fa >/dev/null 2>&1')
    _ = os.system('cat <(cat <(echo ">hap1") <(tail -n +2 ' + outdir + '/hap1.fa' + ')) <(cat <(echo ">hap2") <(tail -n +2 ' +
        outdir + '/hap2.fa' + ')) > ' + outdir + '/gene.fa && ' + samtools + ' faidx ' + outdir + '/gene.fa && rm ' + 
        outdir + '/hap1.fa && rm ' + outdir + '/hap2.fa && rm ' + outdir + '/hap1.vcf.gz && rm ' + outdir + '/hap1.vcf.gz.tbi && rm ' +
        outdir + '/hap2.vcf.gz && rm ' + outdir + '/hap2.vcf.gz.tbi')

def HaplotypeAssignment(samfile):
    '''
    Function designed to retrieve haplotype assignments for each read in samfile
    '''
    assignments = dict()
    for read in pysam.AlignmentFile(samfile, 'r').fetch():
        if read.is_unmapped or read.mapping_quality == 0:
            assignments[read.query_name] = 'NA'
        else:
            assignments[read.query_name] = read.reference_name

    return assignments

def MeasureOverlap(a, b, c, d):
    '''
    Measures the amount of overlap between intervals [a, b] and [c, d]
    '''
    return max(min(b, d) - max(a, c), 0)

def PhasingQuality(bamfile, assignments, start, end):
    '''
    Function designed to assess phasing quality of haplotype read assignments
    '''
    weights = dict()
    for read in pysam.AlignmentFile(bamfile, 'rb').fetch():
        if read.query_name in assignments:
            weights[read.query_name] = MeasureOverlap(start, end, read.reference_start, read.reference_end)/(end - start)
    
    denom = np.sum([weights.get(rname, 0) for rname in assignments.keys()])
    if denom == 0:
        return np.nan
    else:
        return np.sum([weights.get(rname, 0) for rname in assignments.keys() if assignments[rname] != 'NA']) / denom

def WriteBam(samtools, reads, infile, outfile):
    '''
    Function designed to output an indexed BAM file containing the reads in reads
    '''
    inbam = pysam.AlignmentFile(infile, 'rb')
    outbam = pysam.AlignmentFile(outfile, 'wb', template = inbam)

    for read in inbam.fetch():
        if read.query_name in reads:
            outbam.write(read)
        
    outbam.close()
    _ = os.system(samtools + ' index ' + outfile)

def AnnotateGnomad(variants, gnomadfile):
    '''
    Function designed to return AF_GRPMAX_JOINT (maximum allele frequency across gnomAD populations) 
    and NHOMALT_GRPMAX_JOINT (maximum number of homozygous individuals across gnomAD populations) for variants
    '''
    gnomadDF = pd.read_csv(gnomadfile, sep = '\t', header = 0)
    gnomadDF['ID'] = gnomadDF['CHROM'] + '_' + gnomadDF['POS'].astype(str) + '_' + gnomadDF['REF'] + '_' + gnomadDF['ALT']
    return list(zip(pd.DataFrame(variants)[0].map(dict(zip(gnomadDF['ID'], gnomadDF['AF_GRPMAX_JOINT']))).fillna(0), 
        pd.DataFrame(variants)[0].map(dict(zip(gnomadDF['ID'], gnomadDF['NHOMALT_GRPMAX_JOINT']))).fillna(0)))

def AnnotateClinvar(variants, clinvar):
    '''
    Function designed to return CLINSIG for variants if present
    '''
    clinsig = []
    for variant in variants:
        clinsig.append(next(iter([next(iter([item.split('=')[1] for item in record[-1].split(';') if 'CLNSIG=' in item]), '.') for record in 
            tabix.open(clinvar).query(variant[0].replace('chr', ''), int(variant[1])-1, int(variant[1])) if 
            record[3] == variant[2] and record[4] == variant[3]]), '.'))
    return clinsig

def CheckHomopolymer(positions, genome):
    '''
    Function designed to check whether a variant position falls in a homopolymer region
    '''
    homopolymer = []
    for position in positions:
        query = pysam.FastaFile(genome).fetch(position[0], int(position[1]) - 5, int(position[1]) + 4).upper()
        ctr, idx, run, base = False, 0, query[0], query[0] 
        while idx < len(query) - 1:
            if query[idx+1] == base:
                run += query[idx+1]
                if len(run) >= 4:
                    ctr = True
            else:
                run, base = query[idx+1], query[idx+1]
            idx += 1
        homopolymer.append(ctr)
    return homopolymer

def CalculateAlignmentStats(variants, infile):
    '''
    Function designed to compute different alignment statistics for reads supporting a variant, including
    the percentage of variant-supporting reads in which the variant is found within 30 nt from the
    alignment ends and the percentage of reads covering the variant site that are reference skip (introns)
    '''
    values = []
    for variant in variants:
        distances, skipcount, totalcount = [], 0, 0
        for pileupcolumn in pileup_truncated(pysam.AlignmentFile(infile, 'rb'), variant[0], int(variant[1])-1, int(variant[1])):
            for pileupread in pileupcolumn.pileups:
                totalcount += 1
                if not pileupread.is_del and not pileupread.is_refskip:
                    adj = pileupread.alignment.cigartuples[0][1] if pileupread.alignment.cigartuples[0][0] == 4 else 0
                    if len(variant[2]) == len(variant[3]):
                        if pileupread.alignment.query_sequence[pileupread.query_position] == variant[3]:
                            distances.append(min(pileupread.query_position - adj, len(pileupread.alignment.query_alignment_sequence) - 
                                pileupread.query_position + 1))
                    elif len(variant[2]) > len(variant[3]):
                        if pileupread.indel == (len(variant[3]) - len(variant[2])):
                            distances.append(min(pileupread.query_position - adj, len(pileupread.alignment.query_alignment_sequence) - 
                                pileupread.query_position + 1))
                    else:
                        if pileupread.alignment.query_sequence[pileupread.query_position:(pileupread.query_position + len(variant[3]))] == variant[3]:
                            if pileupread.indel == (len(variant[3]) - len(variant[2])):
                                distances.append(min(pileupread.query_position - adj, len(pileupread.alignment.query_alignment_sequence) - 
                                    pileupread.query_position + 1))
                elif pileupread.is_refskip:
                    skipcount += 1
        
        end_dist_pct = sum(np.array(distances) <= 30)/len(distances) if len(distances) > 0 else np.nan
        skip_freq = skipcount/totalcount if totalcount > 0 else np.nan
        values.append((end_dist_pct, skip_freq))
    return values

def GetTotalCoverage(countDF):
    '''
    Function designed to compute total read coverage over splice sites for junctions in countDF
    '''
    coverageDF = countDF[['Chrom', 'Start', 'End']].copy()
    lssDF = countDF.drop('End', axis = 1).groupby(['Chrom', 'Start']).sum().reset_index()
    rssDF = countDF.drop('Start', axis = 1).groupby(['Chrom', 'End']).sum().reset_index()
    samples = [sampid for sampid in countDF.columns if sampid not in {'Chrom', 'Start', 'End'}]

    for clm in samples:
        appendDF = pd.DataFrame(((coverageDF['Chrom'] + '_' + coverageDF['Start'].astype(str)).map(dict(zip(
            lssDF['Chrom'] + '_' + lssDF['Start'].astype(str), lssDF[clm]))).fillna(0) +
            (coverageDF['Chrom'] + '_' + coverageDF['End'].astype(str)).map(dict(zip(
            rssDF['Chrom'] + '_' + rssDF['End'].astype(str), rssDF[clm]))).fillna(0) - countDF[clm]), columns = [clm])
        coverageDF = pd.concat([coverageDF, appendDF], axis = 1)
    
    return coverageDF

def MaskOutliers(x):
    '''
    Function designed to convert "boxplot" outlier values in x to NA
    '''
    if np.sum(~np.isnan(x)) >= 100:
        lower_cutoff = 2.5 * np.nanquantile(x, 0.25) - 1.5 * np.nanquantile(x, 0.75)
        upper_cutoff = 2.5 * np.nanquantile(x, 0.75) - 1.5 * np.nanquantile(x, 0.25)    
        return np.where(np.logical_or(x < lower_cutoff, x > upper_cutoff), np.nan, x)

    else:
        return x

def RescaleValues(x, old_lower, old_upper, new_lower, new_upper):
    '''
    Function designed to rescale values from the interval [old_lower, old_upper]
    to the interval [new_lower, new_upper]
    '''
    return (x - old_lower) * (new_upper - new_lower)/(old_upper - old_lower) + new_lower

def FitBetaDist(x, tol):
    '''
    Function designed to fit a beta distribution on values in x
    '''
    return (np.array(beta.fit(x[~np.isnan(x)], floc = 0, fscale = 1))[:2] if 
        np.nanmax(x) > np.nanmin(x) else np.array([np.nanmin(x), 1 - np.nanmin(x)])/tol)

def BetaBinomialTest(x, n, params):
    '''
    Function designed to return the probability of observing a value x or
    more extreme from a BetaBinomial distribution parameterized on n and params
    '''
    return np.clip(2 * min(betabinom.sf(x - 1, n, params[0], params[1]), betabinom.cdf(x, n, 
        params[0], params[1])) - betabinom.pmf(x, n, params[0], params[1]), 0, 1)

# =====================================================================================================================
#                                                    MAIN FUNCTIONS
# =====================================================================================================================

def QualityControlCheck(infile, targets, gencode, samtools, stringtie, threads, outdir):
    # Create output directory for quality control analysis
    _ = os.system('mkdir -p ' + outdir + '/quality_control')

    # Determine sample name from infile and path to working directory for Run_STRIPE.py
    sampleName = os.path.basename(infile).replace('.bam', '')
    stripeWD = os.path.dirname(os.path.realpath(__file__))

    # Calculate read mapping statistics for infile
    _ = os.system(samtools + ' flagstat --threads ' + str(threads) + ' -O tsv ' + infile + ' > ' + 
        outdir + '/quality_control/' + sampleName + '.mapping_stats.txt')

    # Run stringtie on infile 
    _ = os.system('(python ' + stripeWD + '/utils/Stringtie.py --infile ' + infile + ' --gencode ' + 
        gencode + ' --samtools ' + samtools + ' --stringtie ' + stringtie + ' --threads ' + str(threads) + 
        ' --outprefix ' + outdir + '/quality_control/stringtie/' + sampleName + ') > /dev/null 2>&1')
    
    # Compute on-target rate and identify top 2000 most abundant genes in infile
    _ = os.system('(Rscript ' + stripeWD + '/utils/OnTargetRate.R --gtffile ' + outdir + '/quality_control/stringtie/' + 
        sampleName + '.gtf ' + '--gene_abundance ' +  outdir + '/quality_control/stringtie/' + sampleName + 
        '.gene_abundance.tsv --targets ' + targets + ' --stats ' + outdir + '/quality_control/' + sampleName +
        '.mapping_stats.txt ' + '--gene_output ' + outdir + '/quality_control/' + sampleName + '.top_2000_genes.txt ' + 
        '--plot_pdf ' + outdir + '/quality_control/' + sampleName + '.top_2000_genes.pdf) > /dev/null 2>&1')
    
    # Extract splice junctions from filtered alignments
    ExtractSpliceJunctions(infile, samtools, outdir + '/quality_control/' + sampleName + '.tmp.bam', 
        targets = targets).to_csv(outdir + '/quality_control/' + sampleName + '.junctions.tsv', sep = '\t', index = False)

def PhaseTargets(infile, targets, genome, samtools, longcallr, bcftools, minimap2, threads, phase_indels, outdir,
    proband_vcf, parent1_vcf, parent2_vcf):
    # Create output directory for phasing analysis
    _ = os.system('mkdir -p ' + outdir + '/target_genes')

    # Determine path to working directory for Run_STRIPE.py
    stripeWD = os.path.dirname(os.path.realpath(__file__))

    # Run longcallR on infile
    _ = os.system('(python ' + stripeWD + '/utils/LongcallR.py --infile ' + infile + ' --longcallr ' +
        longcallr + ' --bcftools ' + bcftools + ' --genome ' + genome + ' --outprefix ' + outdir + 
        '/target_genes/longcallR/output --threads ' + str(threads) + ') > /dev/null 2>&1')
    
    # Phase reads mapping to each target gene
    for target in pd.read_csv(targets, sep = '\t', header = None).values:
        _ = os.system('mkdir -p ' + outdir + '/target_genes/' + target[4])

        # Construct shell command to view bi-allelic variants mapping to target gene
        command = bcftools + ' view -H -r ' + target[0] + ':' + str(target[1]+1) + '-' + str(target[2]) + ' -m2 -M2 '
        dna_command = command + '-v snps -f PASS ' if not phase_indels else command + '-f PASS '

        proband_variants = (pd.DataFrame([item.split() for item in subprocess.run(dna_command + proband_vcf,
            shell = True, capture_output = True, text = True).stdout.splitlines()], columns = list(range(10)))
            if proband_vcf is not None else pd.DataFrame(columns = list(range(10))))
        parent1_variants = (pd.DataFrame([item.split() for item in subprocess.run(dna_command + parent1_vcf,
            shell = True, capture_output = True, text = True).stdout.splitlines()], columns = list(range(10)))
            if parent1_vcf is not None else pd.DataFrame(columns = list(range(10))))
        parent2_variants = (pd.DataFrame([item.split() for item in subprocess.run(dna_command + parent2_vcf,
            shell = True, capture_output = True, text = True).stdout.splitlines()], columns = list(range(10)))
            if parent2_vcf is not None else pd.DataFrame(columns = list(range(10))))

        proband_variants['DNA'] = [PullFeature(item, 'GT') for item in proband_variants[[8, 9]].values.tolist()]
        parent1_variants['P1'] = [PullFeature(item, 'GT') for item in parent1_variants[[8, 9]].values.tolist()]
        parent2_variants['P2'] = [PullFeature(item, 'GT') for item in parent2_variants[[8, 9]].values.tolist()]

        # Phase proband variants by transmission (if unambiguous) and filter for heterozygous variants
        dna_variants = PhaseByTransmission(proband_variants, parent1_variants, parent2_variants)
        dna_variants = dna_variants[dna_variants['DNA'].isin({'0/1', '1|0', '0|1'})]

        # Filter longcallR VCF for bi-allelic SNPs mapping to target gene
        rna_variants = pd.DataFrame([item.split() for item in subprocess.run(command + outdir + 
            '/target_genes/longcallR/output.vcf.gz', shell = True, capture_output = True, 
            text = True).stdout.splitlines()], columns = list(range(10)))
        
        # Retreive variant genotype, read depth, allele frequency, and phase set
        rna_variants['RNA'] = [PullFeature(item, 'GT') for item in rna_variants[[8, 9]].values.tolist()]
        rna_variants['DP'] = [PullFeature(item, 'DP') for item in rna_variants[[8, 9]].values.tolist()]
        rna_variants['AF'] = [PullFeature(item, 'AF').split(',')[0] for item in rna_variants[[8, 9]].values.tolist()]
        rna_variants['PS'] = [PullFeature(item, 'PS') for item in rna_variants[[8, 9]].values.tolist()]

        # Annotate variants in rna_variants based on whether they are densely clustered
        rna_variants['Dense'] = CheckDensity(sorted(rna_variants[1].astype(int).tolist()))

        # Filter down rna_variants based on the following criteria:
        #   * Annotated as PASS and not dense
        #   * Annotated as LowQual or RnaEdit but have DP >= 100, AF >= 0.3, and not dense
        rna_variants = pd.concat([rna_variants[(rna_variants[6] == 'PASS') & ~rna_variants['Dense']], rna_variants[rna_variants[6].isin({'LowQual', 'RnaEdit'}) & 
            (rna_variants['DP'].astype(int) >= 100) & (rna_variants['AF'].astype(float) >= 0.3) & ~rna_variants['Dense']]])
        rna_variants = rna_variants[rna_variants['RNA'].isin({'0/1', '1|0', '0|1'})][[0, 1, 3, 4, 'RNA', 'PS']]
        rna_variants.loc[rna_variants['RNA'] == '0/1', 'PS'] = np.nan

        # Merge rna_variants with dna_variants
        merge_variants = MergeVariants(rna_variants, dna_variants)
        
        # Filter infile for primary/supplementary read alignments (MAPQ >= 1) mapping to the target gene
        _ = os.system(samtools + ' view -hb -F 256 -q 1 ' + infile + ' ' + target[0] + ':' + str(target[1]+1) + '-' + 
            str(target[2]) + ' > ' + outdir + '/target_genes/' + target[4] + '/all_reads.bam && ' + samtools +
            ' index ' + outdir + '/target_genes/' + target[4] + '/all_reads.bam')

        # Annotate variants in merge_variants with RNA coverage
        if merge_variants.shape[0] > 0:
            samfile = pysam.AlignmentFile(outdir + '/target_genes/' + target[4] + '/all_reads.bam', mode = 'rb')
            merge_variants['Coverage'] = merge_variants.apply(lambda x: np.sum(samfile.count_coverage(x[0], int(x[1])-1, int(x[1]))), axis = 1)
            merge_variants = merge_variants[merge_variants['Coverage'] >= 20] # Require at least 20 reads over each haplotype-informative variant

            if merge_variants.shape[0] > 0:
                # Identify phase set with the highest coverage
                bestPS = merge_variants[['PS', 'Coverage']].groupby('PS').max().reset_index().sort_values(by = 'Coverage', ascending = False).head(1)['PS'].item()
                merge_variants = merge_variants[merge_variants['PS'] == bestPS][[0, 1, 3, 4, 'RNA']]

                # Construct haplotype-specific sequences for the target gene
                MakeHaplotypeFasta(samtools, bcftools, merge_variants, genome, target[0] + ':' + str(target[1]+1) + '-' + str(target[2]), 
                    outdir + '/target_genes/' + target[4])
                
                # Convert all_reads.bam into a compressed FASTQ file
                _ = os.system('(' + samtools + ' fastq -n -@ ' + str(threads) + ' ' + outdir + '/target_genes/' + target[4] + 
                    '/all_reads.bam | gzip > ' + outdir + '/target_genes/' + target[4] + '/all_reads.fastq.gz) > /dev/null 2>&1')

                # Align all_reads.fastq.gz to gene.fa with minimap2 and determine haplotype assignments for each read
                _ = os.system('(' + minimap2 + ' -t ' + str(threads) + ' -ax splice --secondary=no ' + outdir + '/target_genes/' +
                    target[4] + '/gene.fa ' + outdir + '/target_genes/' + target[4] + '/all_reads.fastq.gz > ' + outdir + 
                    '/target_genes/' + target[4] + '/all_reads.sam' + ') > /dev/null 2>&1')
                _ = os.system('rm ' + outdir + '/target_genes/' + target[4] + '/all_reads.fastq.gz && rm ' + outdir + 
                    '/target_genes/' + target[4] + '/gene.fa && rm ' + outdir + '/target_genes/' + target[4] + '/gene.fa.fai')
                
                # Determine haplotype assignments for each read
                assignments = HaplotypeAssignment(outdir + '/target_genes/' + target[4] + '/all_reads.sam')
                _ = os.system('rm ' + outdir + '/target_genes/' + target[4] + '/all_reads.sam')

                if len({key for key in assignments if assignments[key] != 'NA'}) >= 100:
                    # Compute a phasing quality score for haplotype read assignments and report it along with other phasing statistics
                    qscore = PhasingQuality(outdir + '/target_genes/' + target[4] + '/all_reads.bam', assignments, target[1], target[2])
                    with open(outdir + '/target_genes/' + target[4] + '/phasing_stats.txt', 'w') as quality_file:
                        _ = quality_file.write('Phasing quality score\t%s\n' % format(qscore, '.4f'))
                        _ = quality_file.write('Read count (haplotype 1)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap1'}))
                        _ = quality_file.write('Read count (haplotype 2)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap2'}))
                        _ = quality_file.write('Read count (unassigned)\t%s\n' % len({key for key in assignments if assignments[key] == 'NA'}))
                    
                    # Save read phasing results to haplotype-specific BAM files
                    WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap1'}, 
                        outdir + '/target_genes/' + target[4] + '/all_reads.bam', outdir + '/target_genes/' + target[4] + '/hap1_reads.bam')
                    WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap2'}, 
                        outdir + '/target_genes/' + target[4] + '/all_reads.bam', outdir + '/target_genes/' + target[4] + '/hap2_reads.bam')
                
                else:
                    _ = os.system('rm ' + outdir + '/target_genes/' + target[4] + '/gene.vcf.gz && rm ' + outdir + '/target_genes/' + target[4] + '/gene.vcf.gz.tbi')
        
        # Remove intermediate BAM files
        _ = os.system('rm ' + outdir + '/target_genes/' + target[4] + '/all_reads.bam && rm ' + outdir + '/target_genes/' + target[4] + '/all_reads.bam.bai')

def CallVariants(infile, targets, genome, samtools, bcftools, longcallr, clair3, model, cadd, gnomad, clinvar, threads, outdir):
    # Determine path to working directory for Run_STRIPE.py
    stripeWD = os.path.dirname(os.path.realpath(__file__))

    # Run Clair3 and DeepVariant on infile
    if os.path.isdir(outdir + '/target_genes/clair3'):
        _ = os.system('rm -rf ' + outdir + '/target_genes/clair3')
    if os.path.isdir(outdir + '/target_genes/deepvariant'):
        _ = os.system('rm -rf ' + outdir + '/target_genes/deepvariant')

    _ = os.system('(python ' + stripeWD + '/utils/Clair3.py --clair3 ' + clair3 + ' --infile ' + infile +
        ' --genome ' + genome + ' --threads ' + str(threads) + ' --model ' + model + ' --outdir ' + outdir +
        '/target_genes/clair3) > /dev/null 2>&1')
    _ = os.system('(python ' + stripeWD + '/utils/DeepVariant.py --model_type ONT_R104 --genome ' + genome +
        ' --bam_proband ' + infile + ' --threads ' + str(threads) + ' --outdir ' + outdir + 
        '/target_genes/deepvariant) > /dev/null 2>&1')

    # Iterate over each target gene
    for target in pd.read_csv(targets, sep = '\t', header = None).values:
        _ = os.system('mkdir -p ' + outdir + '/target_genes/' + target[4] + '/variant_calling')

        # Remove output file if it already exists
        if os.path.isfile(outdir + '/target_genes/' + target[4] + '/variant_calling/merged_variants.tsv'):
            os.remove(outdir + '/target_genes/' + target[4] + '/variant_calling/merged_variants.tsv')

        # Construct shell command to view variant calls mapping to target gene for longcallR, Clair3, and DeepVariant
        command = bcftools + ' view -H -r ' + target[0] + ':' + str(target[1]+1) + '-' + str(target[2])

        # Retain bi-allelic variants from longcallR mapping to target gene that meet the following criteria:
        #   * Not located within a dense cluster of variant calls
        #   * DP >= 100 and AF >= 0.3
        longcallr_variants = pd.DataFrame([item.split() for item in subprocess.run(command + ' ' + outdir + 
            '/target_genes/longcallR/output.vcf.gz', shell = True, capture_output = True, text = True).stdout.splitlines()], 
            columns = list(range(10)))
        longcallr_variants = longcallr_variants[(longcallr_variants[3] != longcallr_variants[4])]
        longcallr_variants['Dense'] = CheckDensity(sorted(longcallr_variants[1].astype(int).tolist()))
        longcallr_variants['DP'] = [PullFeature(item, 'DP') for item in longcallr_variants[[8, 9]].values.tolist()]
        longcallr_variants['AF'] = [PullFeature(item, 'AF').split(',')[0] for item in longcallr_variants[[8, 9]].values.tolist()]
        longcallr_variants['GT'] = [PullFeature(item, 'GT').split(',')[0] for item in longcallr_variants[[8, 9]].values.tolist()]
        longcallr_variants = longcallr_variants[longcallr_variants[6].isin({'PASS', 'RnaEdit', 'LowQual'}) & ~longcallr_variants['Dense']
            & ~longcallr_variants[3].str.contains(',') & ~longcallr_variants[4].str.contains(',') & (longcallr_variants['DP'].astype(int) >= 100) &
            (longcallr_variants['AF'].astype(float) >= 0.3)][[0, 1, 3, 4, 'DP', 'AF', 'GT']]
        longcallr_variants.columns = ['CHROM', 'POS', 'REF', 'ALT', 'LONGCALLR_DP', 'LONGCALLR_AF', 'LONGCALLR_GT']
        longcallr_variants.loc[:, 'LONGCALLR_AF'] = longcallr_variants['LONGCALLR_AF'].apply(lambda x: '{:,.2f}'.format(float(x)))

        # Retain bi-allelic variants from Clair3 mapping to target gene that meet the following criteria:
        #   * Not located within a dense cluster of variant calls
        #   * DP >= 100 and AF >= 0.3
        clair3_variants = pd.DataFrame([item.split() for item in subprocess.run(command + ' ' + outdir + 
            '/target_genes/clair3/pileup.vcf.gz', shell = True, capture_output = True, text = True).stdout.splitlines()], 
            columns = list(range(10)))
        clair3_variants = clair3_variants[(clair3_variants[3] != clair3_variants[4]) & clair3_variants[6].isin({'PASS', 'LowQual'})]
        clair3_variants['Dense'] = CheckDensity(sorted(clair3_variants[1].astype(int).tolist()))
        clair3_variants['DP'] = [PullFeature(item, 'DP') for item in clair3_variants[[8, 9]].values.tolist()]
        clair3_variants['AF'] = [PullFeature(item, 'AF').split(',')[0] for item in clair3_variants[[8, 9]].values.tolist()]
        clair3_variants['GT'] = [PullFeature(item, 'GT').split(',')[0] for item in clair3_variants[[8, 9]].values.tolist()]
        clair3_variants = clair3_variants[~clair3_variants['Dense'] & ~clair3_variants[3].str.contains(',') & 
            ~clair3_variants[4].str.contains(',') & (clair3_variants['DP'].astype(int) >= 100) &
            (clair3_variants['AF'].astype(float) >= 0.3)][[0, 1, 3, 4, 'DP', 'AF', 'GT']]
        clair3_variants.columns = ['CHROM', 'POS', 'REF', 'ALT', 'CLAIR3_DP', 'CLAIR3_AF', 'CLAIR3_GT']
        clair3_variants.loc[:, 'CLAIR3_AF'] = clair3_variants['CLAIR3_AF'].apply(lambda x: '{:,.2f}'.format(float(x)))

        # Retain bi-allelic variants from DeepVariant mapping to target gene that meet the following criteria:
        #   * Not located within a dense cluster of variant calls
        #   * DP >= 100 and AF >= 0.3
        dv_variants = pd.DataFrame([item.split() for item in subprocess.run(command + ' ' + outdir + 
            '/target_genes/deepvariant/PROBAND.vcf.gz', shell = True, capture_output = True, text = True).stdout.splitlines()], 
            columns = list(range(10)))
        dv_variants = dv_variants[(dv_variants[3] != dv_variants[4]) & dv_variants[6].isin({'PASS', 'LowQual'})]
        dv_variants['Dense'] = CheckDensity(sorted(dv_variants[1].astype(int).tolist()))
        dv_variants['DP'] = [PullFeature(item, 'DP') for item in dv_variants[[8, 9]].values.tolist()]
        dv_variants['VAF'] = [PullFeature(item, 'VAF').split(',')[0] for item in dv_variants[[8, 9]].values.tolist()]
        dv_variants['GT'] = [PullFeature(item, 'GT').split(',')[0] for item in dv_variants[[8, 9]].values.tolist()]
        dv_variants = dv_variants[~dv_variants['Dense'] & ~dv_variants[3].str.contains(',') & 
            ~dv_variants[4].str.contains(',') & (dv_variants['DP'].astype(int) >= 100) &
            (dv_variants['VAF'].astype(float) >= 0.3)][[0, 1, 3, 4, 'DP', 'VAF', 'GT']]
        dv_variants.columns = ['CHROM', 'POS', 'REF', 'ALT', 'DEEPVARIANT_DP', 'DEEPVARIANT_AF', 'DEEPVARIANT_GT']
        dv_variants.loc[:, 'DEEPVARIANT_AF'] = dv_variants['DEEPVARIANT_AF'].apply(lambda x: '{:,.2f}'.format(float(x)))

        # Merge longcallr_variants, clair3_variants, and dv_variants into a unified callset
        callset = longcallr_variants.merge(clair3_variants, on = ['CHROM', 'POS', 'REF', 'ALT'], how = 'outer').merge(
            dv_variants, on = ['CHROM', 'POS', 'REF', 'ALT'], how = 'outer').sort_values(by = ['CHROM', 'POS']).fillna('.')

        # Check if phasing results exist for current gene
        if os.path.isfile(outdir + '/target_genes/' + target[4] + '/phasing_stats.txt'):
            # Perform haploid variant calling on the phased read sets using samtools/bcftools
            _ = os.system('(' + bcftools + ' mpileup -Ou -f ' + genome + ' -r ' + target[0] + ':' + str(target[1]+1) + '-' + 
                str(target[2]) + ' -B -Q 5 ' + outdir + '/target_genes/' + target[4] + '/hap1_reads.bam ' +
                outdir + '/target_genes/' + target[4] + '/hap2_reads.bam | ' + bcftools + ' call -mv -Oz ' +
                '--ploidy 1 -o ' + outdir + '/target_genes/' + target[4] + '/variant_calling/bcftools.vcf.gz) > /dev/null 2>&1')
            
            # Retain bi-allelic variants from bcftools mapping to target gene that meet the following criteria:
            #   * Not located within a dense cluster of variant calls
            #   * DP1 >= 30 or DP2 >= 30
            phasing_variants = pd.DataFrame([item.split() for item in subprocess.run(bcftools + ' view -H ' + outdir + 
                '/target_genes/' + target[4] + '/variant_calling/bcftools.vcf.gz', shell = True, capture_output = True, 
                text = True).stdout.splitlines()], columns = list(range(11)))
            phasing_variants = phasing_variants[phasing_variants[3] != phasing_variants[4]]
            phasing_variants['Dense'] = CheckDensity(sorted(phasing_variants[1].astype(int).tolist()))
            phasing_variants = phasing_variants[~phasing_variants['Dense'] & ~phasing_variants[3].str.contains(',') & ~phasing_variants[4].str.contains(',')]
            phasing_variants['REF1'] = [PullFeature(item, 'AD').split(',')[0] for item in phasing_variants[[8, 9]].values.tolist()]
            phasing_variants['ALT1'] = [PullFeature(item, 'AD').split(',')[1] for item in phasing_variants[[8, 9]].values.tolist()]
            phasing_variants['GT1'] = [PullFeature(item, 'GT').split(',')[0] for item in phasing_variants[[8, 9]].values.tolist()]
            phasing_variants['REF2'] = [PullFeature(item, 'AD').split(',')[0] for item in phasing_variants[[8, 10]].values.tolist()]
            phasing_variants['ALT2'] = [PullFeature(item, 'AD').split(',')[1] for item in phasing_variants[[8, 10]].values.tolist()]
            phasing_variants['GT2'] = [PullFeature(item, 'GT').split(',')[0] for item in phasing_variants[[8, 10]].values.tolist()]
            phasing_variants['DP1'] = phasing_variants['REF1'].astype(int) + phasing_variants['ALT1'].astype(int)
            phasing_variants['DP2'] = phasing_variants['REF2'].astype(int) + phasing_variants['ALT2'].astype(int)
            phasing_variants['DP'] = (phasing_variants['DP1'] + phasing_variants['DP2']).astype(object)
            phasing_variants = phasing_variants[(phasing_variants['DP1'] >= 30) | (phasing_variants['DP2'] >= 30)]
            phasing_variants['AF'] = ((phasing_variants['ALT1'].astype(int) + phasing_variants['ALT2'].astype(int))/phasing_variants['DP']).apply(lambda x: '{:,.2f}'.format(float(x)))
            
            if phasing_variants.shape[0] > 0:
                phasing_variants['GT'] = phasing_variants['GT1'].replace('.', '0') + '|' + phasing_variants['GT2'].replace('.', '0')
                phasing_variants = phasing_variants[[0, 1, 3, 4, 'DP', 'AF', 'GT']]
                phasing_variants.columns = ['CHROM', 'POS', 'REF', 'ALT', 'PHASING_DP', 'PHASING_AF', 'PHASING_GT']

            else:
                phasing_variants = pd.DataFrame(columns = ['CHROM', 'POS', 'REF', 'ALT', 'PHASING_DP', 'PHASING_AF', 'PHASING_GT'])

            # Merge phasing_variants with callset
            callset = callset.merge(phasing_variants, on = ['CHROM', 'POS', 'REF', 'ALT'], how = 'outer').sort_values(by = ['CHROM', 'POS']).fillna('.')
        
        else:
            callset = callset.merge(pd.DataFrame(columns = ['CHROM', 'POS', 'REF', 'ALT', 'PHASING_DP', 'PHASING_AF', 'PHASING_GT']), 
                on = ['CHROM', 'POS', 'REF', 'ALT'], how = 'outer').sort_values(by = ['CHROM', 'POS']).fillna('.')

        if callset.shape[0] > 0:
            # Score variants in callset using CADD
            _ = os.system('mkdir -p ' + outdir + '/target_genes/' + target[4] + '/variant_calling/CADD')
            pd.DataFrame([[item[0].replace('chr', ''), item[1], '.', item[2], item[3]] for item in callset.values]).to_csv(outdir + 
                '/target_genes/' + target[4] + '/variant_calling/CADD/input.vcf', sep = '\t', header = False, index = False)
            _ = os.system('(' + cadd + ' -o ' + outdir + '/target_genes/' + target[4] + '/variant_calling/CADD/output.tsv.gz ' + outdir + 
                '/target_genes/' + target[4] + '/variant_calling/CADD/input.vcf) > /dev/null 2>&1')
            if os.path.isfile(outdir + '/target_genes/' + target[4] + '/variant_calling/CADD/output.tsv.gz'):
                cadd_results = pd.read_csv(outdir + '/target_genes/' + target[4] + '/variant_calling/CADD/output.tsv.gz', header = None, sep = '\t', 
                    compression = 'gzip', names = ['Chrom', 'Pos', 'Ref', 'Alt', 'RawScore', 'PHRED'], comment = '#')
                callset['CADD'] = ((callset['CHROM'] + '_' + callset['POS'].astype(str) + '_' + callset['REF'] + '_' + callset['ALT']).map(dict(zip(
                    'chr' + cadd_results['Chrom'].astype(str) + '_' + cadd_results['Pos'].astype(str) + '_' + cadd_results['Ref'] + '_' + cadd_results['Alt'], 
                    cadd_results['PHRED'])))).apply(lambda x: '{:,.1f}'.format(float(x)))
                _ = os.system('rm ' + outdir + '/target_genes/' + target[4] + '/variant_calling/CADD/input.vcf && rm ' + outdir + '/target_genes/' + 
                    target[4] + '/variant_calling/CADD/output.tsv.gz && rmdir ' + outdir + '/target_genes/' + target[4] + '/variant_calling/CADD')
            else:
                callset['CADD'] = 0
                _ = os.system('rm ' + outdir + '/target_genes/' + target[4] + '/variant_calling/CADD/input.vcf && rmdir ' + outdir + '/target_genes/' + 
                    target[4] + '/variant_calling/CADD')

            # Score variants in callset using SpliceAI
            _ = os.system('mkdir -p ' + outdir + '/target_genes/' + target[4] + '/variant_calling/SpliceAI')
            _ = os.system('echo "##fileformat=VCFv4.2" > ' + outdir + '/target_genes/' + target[4] + '/variant_calling/SpliceAI/input.vcf')
            _ = os.system(bcftools + ' view -h ' + outdir + '/target_genes/longcallR/output.vcf.gz | grep "##contig=" >> ' + outdir + 
                '/target_genes/' + target[4] + '/variant_calling/SpliceAI/input.vcf')
            pd.DataFrame([[item[0], item[1], '.', item[2], item[3], '.', '.', '.'] for item in callset.values], columns = ['#CHROM', 'POS', 'ID',
                'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']).to_csv(outdir + '/target_genes/' + target[4] + '/variant_calling/SpliceAI/input.vcf', 
                mode = 'a', sep = '\t', index = False)
            _ = os.system('(spliceai -I ' + outdir + '/target_genes/' + target[4] + '/variant_calling/SpliceAI/input.vcf -O ' +
                outdir + '/target_genes/' + target[4] + '/variant_calling/SpliceAI/output.vcf -R ' + genome + ' -A grch38) > /dev/null 2>&1')
            if os.path.isfile(outdir + '/target_genes/' + target[4] + '/variant_calling/SpliceAI/output.vcf'):
                spliceai_results = pd.read_csv(outdir + '/target_genes/' + target[4] + '/variant_calling/SpliceAI/output.vcf', header = None, sep = '\t',
                    names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'], comment = '#')
                spliceai_results['SPLICEAI'] = spliceai_results['INFO'].apply(lambda x: np.max([float(val) if val != '.' else 0 for val in x.split('|')[2:6]]) if 'SpliceAI' in x else np.nan)
                callset['SPLICEAI'] = ((callset['CHROM'] + '_' + callset['POS'].astype(str) + '_' + callset['REF'] + '_' + callset['ALT']).map(dict(zip(
                    spliceai_results['CHROM'] + '_' + spliceai_results['POS'].astype(str) + '_' + spliceai_results['REF'] + '_' + spliceai_results['ALT'],
                    spliceai_results['SPLICEAI'])))).apply(lambda x: '{:,.2f}'.format(float(x)))
                _ = os.system('rm ' + outdir + '/target_genes/' + target[4] + '/variant_calling/SpliceAI/input.vcf && rm ' + outdir + '/target_genes/' + 
                    target[4] + '/variant_calling/SpliceAI/output.vcf && rmdir ' + outdir + '/target_genes/' + target[4] + '/variant_calling/SpliceAI')
            else:
                callset['SPLICEAI'] = 0
                _ = os.system('rm ' + outdir + '/target_genes/' + target[4] + '/variant_calling/SpliceAI/input.vcf && rmdir ' + outdir + '/target_genes/' + 
                    target[4] + '/variant_calling/SpliceAI')
            
            # Only keep variants with CADD > 15 or SpliceAI > 0.1
            callset = callset[(callset['CADD'].astype(float) > 15) | (callset['SPLICEAI'].astype(float) > 0.1)]

            if callset.shape[0] > 0:
                # Annotate variants with information from gnomAD (maximum allele frequency across populations and maximum number of
                # homozygous individuals across populations) and only keep variants where GNOMAD_AF_GRPMAX_JOINT < 1%
                callset['GNOMAD_AF_GRPMAX_JOINT'], callset['GNOMAD_NHOMALT_GRPMAX_JOINT'] = zip(*AnnotateGnomad((callset['CHROM'] + '_' + callset['POS'].astype(str) + 
                    '_' + callset['REF'] + '_' + callset['ALT']).values, gnomad + '/' + target[4] + '/gnomad_variants.tsv'))
                callset = callset[callset['GNOMAD_AF_GRPMAX_JOINT'] < 0.01]
                callset.loc[:, 'GNOMAD_AF_GRPMAX_JOINT'] = callset['GNOMAD_AF_GRPMAX_JOINT'].apply(lambda x: '{:,.3e}'.format(x))
                callset.loc[:, 'GNOMAD_NHOMALT_GRPMAX_JOINT'] = callset['GNOMAD_NHOMALT_GRPMAX_JOINT'].astype(object)

                # Annotate variants with information from ClinVar
                callset['CLINVAR'] = AnnotateClinvar(callset[['CHROM', 'POS', 'REF', 'ALT']].values, clinvar)

                # Check whether variant falls in homopolymer region
                callset['HOMOPOLYMER'] = CheckHomopolymer(callset[['CHROM', 'POS']].values, genome)
                
                if callset.shape[0] > 0:
                    # Compute alignment statistics for reads supporting the variant
                    callset['END_DISTANCE_PCT'], callset['SKIP_FREQ'] = zip(*CalculateAlignmentStats(callset[['CHROM', 'POS', 'REF', 'ALT']].values, infile))
                    callset.loc[:, 'END_DISTANCE_PCT'] = callset['END_DISTANCE_PCT'].apply(lambda x: '{:,.3f}'.format(x))
                    callset.loc[:, 'SKIP_FREQ'] = callset['SKIP_FREQ'].apply(lambda x: '{:,.3f}'.format(x))
                    
                    # Save callset to an output file
                    callset.to_csv(outdir + '/target_genes/' + target[4] + '/variant_calling/merged_variants.tsv', sep = '\t', index = False) 

def AberrantJunctions(infile, targets, samtools, gtex_splice, outdir):
    # Read in gtex_splice as a dataframe
    gtexDF = pd.read_csv(gtex_splice, sep = '\t', header = 0)
    gtexDF[['Chrom', 'Start', 'End']] = gtexDF['Name'].str.split('_', n = 2, expand = True)
    gtexDF = gtexDF[['Chrom', 'Start', 'End'] + list(gtexDF.columns[1:-3])]
    
    # Iterate over each target gene
    for target in pd.read_csv(targets, sep = '\t', header = None).values:
        _ = os.system('mkdir -p ' + outdir + '/target_genes/' + target[4] + '/aberrant_splicing')

        # Remove output files if they exist already
        for suffix in ['all', 'hap1', 'hap2']:
            currfile = outdir + '/target_genes/' + target[4] + '/aberrant_splicing/output_' + suffix + '.txt'
            if os.path.isfile(currfile):
                os.remove(currfile)
        
        # Retrieve junctions mapping to target gene from infile
        allJunctions = ExtractSpliceJunctions(infile, samtools, outdir + '/target_genes/' + target[4] + '/aberrant_splicing/all_reads.bam', 
            region = target[0] + ':' + str(target[1]+1) + '-' + str(target[2]))
        allJunctions.columns = ['Chrom', 'Start', 'End', 'All']

        if allJunctions.shape[0] > 0:
            # Check if phasing results exist for current gene
            if os.path.isfile(outdir + '/target_genes/' + target[4] + '/phasing_stats.txt'):
                # Retrieve junctions mapping to target gene from each haplotype-specific BAM
                hapJunctions1 = ExtractSpliceJunctions(outdir + '/target_genes/' + target[4] + '/hap1_reads.bam', samtools, outdir + '/target_genes/' + 
                    target[4] + '/aberrant_splicing/hap1_reads.bam', region = target[0] + ':' + str(target[1]+1) + '-' + str(target[2]))
                hapJunctions1.columns = ['Chrom', 'Start', 'End', 'Hap1']

                hapJunctions2 = ExtractSpliceJunctions(outdir + '/target_genes/' + target[4] + '/hap2_reads.bam', samtools, outdir + '/target_genes/' + 
                    target[4] + '/aberrant_splicing/hap2_reads.bam', region = target[0] + ':' + str(target[1]+1) + '-' + str(target[2]))
                hapJunctions2.columns = ['Chrom', 'Start', 'End', 'Hap2']

                # Merge hapJunctions1 and hapJunctions2 with allJunctions
                allJunctions = allJunctions.merge(hapJunctions1, on = ['Chrom', 'Start', 'End'], how = 'outer').merge(hapJunctions2, 
                    on = ['Chrom', 'Start', 'End'], how = 'outer').sort_values(by = ['Chrom', 'Start', 'End']).fillna(0)
            
            # Compute total read coverage for junctions in allJunctions
            allCoverage = GetTotalCoverage(allJunctions)

            # Only consider junctions where at least one input sample has at least 20 junction spanning reads
            countFilter = list(allJunctions[allJunctions.columns[3:]].max(axis = 1) >= 20)
            allJunctions, allCoverage = allJunctions[countFilter], allCoverage[countFilter]

            if sum(countFilter) > 0:
                # Retrieve read counts and total read coverage for junctions in allJunctions across GTEx samples
                ctrlCounts = gtexDF[(gtexDF['Chrom'] == target[0]) & (gtexDF['Start'].isin(set(allJunctions['Start'])) | 
                    gtexDF['End'].isin(set(allJunctions['End'])))].drop_duplicates().reset_index(drop = True)
                ctrlCounts = allJunctions[['Chrom', 'Start', 'End']].merge(ctrlCounts, on = ['Chrom', 'Start', 'End'], 
                    how = 'outer').sort_values(by = ['Chrom', 'Start', 'End']).fillna(0)
                ctrlCoverage = GetTotalCoverage(ctrlCounts)
                ctrlCounts = allJunctions[['Chrom', 'Start', 'End']].merge(ctrlCounts, on = ['Chrom', 'Start', 'End'], how = 'left')
                ctrlCoverage = allJunctions[['Chrom', 'Start', 'End']].merge(ctrlCoverage, on = ['Chrom', 'Start', 'End'], how = 'left')

                # Compute PSI values for each splice junction across GTEx samples (require a total coverage of at least 20 reads)
                ctrlPSI = pd.concat([allJunctions[['Chrom', 'Start', 'End']].reset_index(drop = True), 
                    (ctrlCounts.iloc[:, 3:]/ctrlCoverage.iloc[:, 3:].mask(ctrlCoverage.iloc[:, 3:] < 20)).reset_index(drop = True)], axis = 1)
                
                # Mask any PSI values in ctrlPSI that are "boxplot" outliers in each row
                ctrlPSI.iloc[:, 3:] = pd.DataFrame(ctrlPSI.iloc[:, 3:].apply(lambda x: MaskOutliers(np.array(x, dtype = float)), axis = 1).tolist())

                # Drop junctions from ctrlPSI, allJunctions, and allCoverage where you have less than 100 GTEx samples with defined PSI values
                na_filter = list((~ctrlPSI.iloc[:, 3:].isna()).sum(axis = 1) >= 100)
                ctrlPSI, allJunctions, allCoverage = ctrlPSI[na_filter], allJunctions[na_filter], allCoverage[na_filter]

                if sum(na_filter) > 0:
                    # Rescale PSI values in ctrlPSI from the interval [0, 1] to [1e-3, 1-1e-3]
                    ctrlPSI.iloc[:, 3:] = np.round(RescaleValues(ctrlPSI.iloc[:, 3:], 0, 1, 1e-3, 1-1e-3), 3)

                    # Fit beta distributions on junctions represented in ctrlPSI
                    parameters = list(ctrlPSI.iloc[:, 3:].apply(lambda x: FitBetaDist(np.array(x, dtype = float), 1e-3), axis = 1))

                    # Run a Beta-Binomial Test on junctions based on all reads and save to output folder
                    allDF = pd.DataFrame({
                        'Junction': list(allJunctions['Chrom'] + ':' + allJunctions['Start'].astype(str) + '-' + allJunctions['End'].astype(str)),
                        'Count': list(allJunctions['All'].astype(int).astype(object)), 'Coverage': list(allCoverage['All'].astype(int).astype(object)),
                        'Usage': list((allJunctions['All']/allCoverage['All'].mask(allCoverage['All'] < 20)).apply(lambda x: '{:,.3f}'.format(x))),
                        'Population': ['{:,.3f}'.format(item[0]/(item[0] + item[1])) for item in parameters],
                        'Shift': ['{:,.3f}'.format(allJunctions.iloc[idx]['All']/allCoverage['All'].mask(allCoverage['All'] < 20).iloc[idx] - 
                            parameters[idx][0]/np.sum(parameters[idx])) for idx in range(len(parameters))],
                        'PVal': ['{:,.3e}'.format(BetaBinomialTest(allJunctions.iloc[idx]['All'], allCoverage.iloc[idx]['All'],
                            parameters[idx])) for idx in range(len(parameters))]})
                    allDF[allDF['Count'].astype(int) >= 20].to_csv(outdir + '/target_genes/' + target[4] + '/aberrant_splicing/output_all.txt', 
                        sep = '\t', index = False, na_rep = 'NA')

                    # Run Beta-Binomial Tests on junctions based on haplotype-resolved read counts and save to output folder
                    if 'Hap1' in allJunctions.columns and 'Hap2' in allJunctions.columns:
                        hap1DF = pd.DataFrame({
                            'Junction': list(allJunctions['Chrom'] + ':' + allJunctions['Start'].astype(str) + '-' + allJunctions['End'].astype(str)),
                            'Count': list(allJunctions['Hap1'].astype(int).astype(object)), 'Coverage': list(allCoverage['Hap1'].astype(int).astype(object)),
                            'Usage': list((allJunctions['Hap1']/allCoverage['Hap1'].mask(allCoverage['Hap1'] < 20)).apply(lambda x: '{:,.3f}'.format(x))),
                            'Population': ['{:,.3f}'.format(item[0]/(item[0] + item[1])) for item in parameters],
                            'Shift': ['{:,.3f}'.format(allJunctions.iloc[idx]['Hap1']/allCoverage['Hap1'].mask(allCoverage['Hap1'] < 20).iloc[idx] - 
                                parameters[idx][0]/np.sum(parameters[idx])) for idx in range(len(parameters))],
                            'PVal': ['{:,.3e}'.format(BetaBinomialTest(allJunctions.iloc[idx]['Hap1'], allCoverage.iloc[idx]['Hap1'],
                                parameters[idx])) for idx in range(len(parameters))]})
                        hap1DF[hap1DF['Count'].astype(int) >= 20].to_csv(outdir + '/target_genes/' + target[4] + '/aberrant_splicing/output_hap1.txt', 
                            sep = '\t', index = False, na_rep = 'NA')
                        
                        hap2DF = pd.DataFrame({
                            'Junction': list(allJunctions['Chrom'] + ':' + allJunctions['Start'].astype(str) + '-' + allJunctions['End'].astype(str)),
                            'Count': list(allJunctions['Hap2'].astype(int).astype(object)), 'Coverage': list(allCoverage['Hap2'].astype(int).astype(object)),
                            'Usage': list((allJunctions['Hap2']/allCoverage['Hap2'].mask(allCoverage['Hap2'] < 20)).apply(lambda x: '{:,.3f}'.format(x))),
                            'Population': ['{:,.3f}'.format(item[0]/(item[0] + item[1])) for item in parameters],
                            'Shift': ['{:,.3f}'.format(allJunctions.iloc[idx]['Hap2']/allCoverage['Hap2'].mask(allCoverage['Hap2'] < 20).iloc[idx] - 
                                parameters[idx][0]/np.sum(parameters[idx])) for idx in range(len(parameters))],
                            'PVal': ['{:,.3e}'.format(BetaBinomialTest(allJunctions.iloc[idx]['Hap2'], allCoverage.iloc[idx]['Hap2'],
                                parameters[idx])) for idx in range(len(parameters))]})
                        hap2DF[hap2DF['Count'].astype(int) >= 20].to_csv(outdir + '/target_genes/' + target[4] + '/aberrant_splicing/output_hap2.txt', 
                            sep = '\t', index = False, na_rep = 'NA')

def HaplotypeDosage(targets, gtex_haplotype, outdir):
    # Read in gtex_haplotype as a dataframe
    gtexDF = pd.read_csv(gtex_haplotype, sep = '\t', header = 0)
    gtexDF.loc[:, 'name'] = gtexDF['name'].str.split('.').str[0]

    # Iterate over each target gene
    for target in pd.read_csv(targets, sep = '\t', header = None).values:
        _ = os.system('mkdir -p ' + outdir + '/target_genes/' + target[4] + '/haplotype_dosage')

        # Remove output file if it exists already
        if os.path.isfile(outdir + '/target_genes/' + target[4] + '/haplotype_dosage/output.txt'):
            os.remove(outdir + '/target_genes/' + target[4] + '/haplotype_dosage/output.txt')

        # Check if phasing results exist for current gene
        if os.path.isfile(outdir + '/target_genes/' + target[4] + '/phasing_stats.txt'):
            # Read in phasing_stats.txt and pull out phasing score and read counts for haplotype 1 and 2
            phasingInfo = list(pd.read_csv(outdir + '/target_genes/' + target[4] + '/phasing_stats.txt',
                sep = '\t', header = None).iloc[:3, 1])
            
            # Check if total number of haplotype-assigned reads is at least 100
            if phasingInfo[1] + phasingInfo[2] >= 100:
                # Retrieve haplotype read counts for current gene across GTEx samples
                ctrlDF = gtexDF[gtexDF['name'] == target[3]].iloc[:, 1:]
                if ctrlDF.shape[0] > 0:
                    ctrlCounts = [tuple(map(int, item.split('|'))) for item in ctrlDF.values[0]]

                    # Compute haplotype expression ratios for GTEx samples (require a total count of at least 20)
                    ctrlRatios = [min(item)/sum(item) if sum(item) >= 20 else np.nan for item in ctrlCounts]

                    # Mask any haplotype expression ratios that are "boxplot" outliers
                    ctrlRatios = MaskOutliers(np.array(ctrlRatios))
                    
                    if np.sum(~np.isnan(ctrlRatios)) >= 100:
                        # Rescale haplotype expression ratios from the interval [0, 1] to [1e-3, 1-1e-3]
                        ctrlRatios = np.round(RescaleValues(ctrlRatios, 0, 1, 1e-3, 1-1e-3), 3)

                        # Fit ctrlRatios to a Beta Distribution with a fixed mean of 0.5
                        myFunc = lambda x: 2*(digamma(x)-digamma(2*x)) - np.nanmean(np.log(ctrlRatios*(1-ctrlRatios)))

                        if np.sign(myFunc(1)) != np.sign(myFunc(999)):
                            fit = brentq(myFunc, 1, 999, full_output = True, disp = False)
                            
                            if fit[1].converged:
                                # Run Beta Binomial Test on haplotype read counts and save to output folder
                                pval = BetaBinomialTest(min(phasingInfo[1], phasingInfo[2]), phasingInfo[1] + phasingInfo[2], [fit[0]] * 2)
                                with open(outdir + '/target_genes/' + target[4] + '/haplotype_dosage/output.txt', 'w') as outfile:
                                    _ = outfile.write('Read count (haplotype 1)\t%s\n' % int(phasingInfo[1]))
                                    _ = outfile.write('Read count (haplotype 2)\t%s\n' % int(phasingInfo[2]))
                                    _ = outfile.write('Phasing quality score\t%s\n' % format(phasingInfo[0], '.4f'))
                                    _ = outfile.write('Haplotype dosage (p-value)\t%s\n' % format(pval, '.3e'))

def main():
    message = 'Runs computational pipeline for STRIPE'
    parser = argparse.ArgumentParser(description = message)

    # Add arguments
    parser.add_argument('--infile', metavar = '/path/to/input/TEQUILA/bam', required = True,
        help = 'path to input file (BAM) with TEQUILA-seq read alignments')
    parser.add_argument('--targets', metavar = '/path/to/target/gene/BED', required = True,
        help = 'path to BED file with coordinates of target genes')
    parser.add_argument('--genome', metavar = '/path/to/reference/genome/FASTA', required = True,
        help = 'path to FASTA file with reference genome sequence')
    parser.add_argument('--gencode', metavar = '/path/to/gencode/GTF', required = True,
        help = 'path to GTF file with gene annotations from GENCODE')
    parser.add_argument('--samtools', metavar = '/path/to/samtools', required = True,
        help = 'path to samtools executable')
    parser.add_argument('--stringtie', metavar = '/path/to/stringtie', required = True,
        help = 'path to stringtie executable')
    parser.add_argument('--longcallr', metavar = '/path/to/longcallr', required = True,
        help = 'path to longcallR executable')
    parser.add_argument('--clair3', metavar = '/path/to/clair3/sif', required = True,
        help = 'path to clair3 singularity image file')
    parser.add_argument('--model', metavar = '/path/to/clair3/model', required = True,
        help = 'path to clair3 model file')
    parser.add_argument('--cadd', metavar = '/path/to/CADD', required = True,
        help = 'path to CADD.sh script')
    parser.add_argument('--gnomad', metavar = '/path/to/gnomad/database', required = True,
        help = 'path to gnomAD database (from running Build_STRIPE.py)')
    parser.add_argument('--clinvar', metavar = '/path/to/clinvar/vcf', required = True,
        help = 'path to VCF file with ClinVar variants')
    parser.add_argument('--bcftools', metavar = '/path/to/bcftools', required = True,
        help = 'path to bcftools executable')
    parser.add_argument('--minimap2', metavar = '/path/to/minimap2', required = True,
        help = 'path to minimap2 executable')
    parser.add_argument('--gtex_splice', metavar = '/path/to/gtex/splice/database', required = True,
        help = 'path to database of splice junction read counts for GTEx controls (from running Build_STRIPE.py)')
    parser.add_argument('--gtex_haplotype', metavar = '/path/to/gtex/haplotype/database', required = True,
        help = 'path to database of haplotype expression for GTEx controls (from running Build_STRIPE.py)')
    parser.add_argument('--threads', metavar = '<number of threads>', required = False, 
        help = 'number of threads', type = int, default = 1)
    parser.add_argument('--proband_vcf', metavar = '/path/to/proband/vcf', required = False,
        help = 'path to VCF file for proband (based on prior WES/WGS)')
    parser.add_argument('--parent1_vcf', metavar = '/path/to/parent1/vcf', required = False,
        help = 'path to VCF file for parent 1 (based on prior WES/WGS)')
    parser.add_argument('--parent2_vcf', metavar = '/path/to/parent2/vcf', required = False,
        help = 'path to VCF file for parent 2 (based on prior WES/WGS)')
    parser.add_argument('--only_qc', action = 'store_true', help = 'only perform quality control check for STRIPE')
    parser.add_argument('--phase_indels', action = 'store_true', help = 'use indels from WES/WGS data to phase reads')
    parser.add_argument('--outdir', metavar = '/path/to/output/directory', required = True,
        help = 'path to output directory for STRIPE results')
    
    # Parse command line arguments
    args = parser.parse_args()

    # Create output directory for STRIPE
    _ = os.system('mkdir -p ' + args.outdir)

    # Perform quality control analysis of input TEQUILA-seq data
    QualityControlCheck(args.infile, args.targets, args.gencode, args.samtools, args.stringtie, 
        args.threads, args.outdir)
    
    if not args.only_qc:    
        # Phase TEQUILA-seq reads for each target gene
        PhaseTargets(args.infile, args.targets, args.genome, args.samtools, args.longcallr,
            args.bcftools, args.minimap2, args.threads, args.phase_indels, args.outdir, args.proband_vcf,
            args.parent1_vcf, args.parent2_vcf)
        
        # Call rare deleterious sequence variants from TEQUILA-seq reads for each target gene
        CallVariants(args.infile, args.targets, args.genome, args.samtools, args.bcftools, args.longcallr, 
            args.clair3, args.model, args.cadd, args.gnomad, args.clinvar, args.threads, args.outdir)
        
        # Detect aberrant splice junctions from TEQUILA-seq reads for each target gene
        AberrantJunctions(args.infile, args.targets, args.samtools, args.gtex_splice, args.outdir)

        # Assess whether target genes show unusually strong haplotype dosage imbalance
        HaplotypeDosage(args.targets, args.gtex_haplotype, args.outdir)

if __name__ == '__main__':
    main()
