#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.05
# Figure 6a

# (a) Sashimi plot visualization of haplotype-resolved MTFMT splicing patterns in individual PMD-P23.

# We used the following Python code to phase reads for MTFMT. Analysis of unphased TEQUILA-seq reads for individual
# PMD-P23 uncovered two splicing defects: exon 4 skipping and inclusion of a 53 bp pseudoexon between exons 6 and 7.
# We adopted the following "pseudo" phasing strategy to determine whether these two splicing defects are in trans
# (unfortunately, this individual does not appear to have any other heterozygous exonic variants that can be used to
# inform read phasing): 
#   1. Split up the reads into the following groups/haplotypes: reads that feature canonical splicing of intron 6 and
#      reads that feature inclusion of the 53 bp pseudoexon. 
#   2. If exon 4 skipping is only observed in reads that feature canonical splicing of intron 6, then the two splicing
#      defects are indeed in trans

import sys

# Establish working directories and file paths
workdir = '/mnt/isilon/lin_lab_share/STRIPE'
samtools = '/scr1/users/wangr5/tools/samtools-1.21/samtools'

# Import functions used in Run_STRIPE.py
sys.path.append(workdir + '/scripts/main')
from Run_STRIPE import *

# Read in mapping file linking patient IDs to their original IDs
sampleMap = pd.read_csv(workdir + '/manuscript/Revisions/20250228/sample.map', sep = '\t', header = None)
sampleMap = dict(zip(sampleMap[1], sampleMap[0]))

patient_id, gene_panel, gene_name = sampleMap['PMD-P23'], 'PMD-359', 'MTFMT'

targetDF = pd.read_csv(workdir + '/' + gene_panel.split('-')[0] + '/references/target_genes.bed', sep = '\t', header = None)
target, assignments = targetDF[targetDF[4] == gene_name].values[0], dict()
_ = os.system(samtools + ' view -hb -F 256 -q 1 ' + workdir + '/' + gene_panel.split('-')[0] +
    '/' + patient_id + '/RNA/' + patient_id + '_TEQUILA.bam' + ' ' + target[0] + ':' + str(target[1]+1) + 
    '-' + str(target[2]) + ' > ' + workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam && ' + 
    samtools + ' index ' + workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam')
samfile = pysam.AlignmentFile(workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam', 'rb')

for read in samfile.fetch():
    blocks = re.findall('(\d+)(\D+)', read.cigarstring)
    cumlen = np.cumsum([0] + [int(item[0]) if item[1] in {'M', 'D', 'N'} else 0 for item in blocks][:-1])
    blocks = [(blocks[i][0], blocks[i][1], cumlen[i] + int(read.reference_start) + 1) for i in range(len(blocks))]
    junctions = [(cumlen[i] + int(read.reference_start) + 1, cumlen[i+1] + int(read.reference_start)) for i in range(len(blocks)) if blocks[i][1] == 'N']
    if (65006192, 65016435) in junctions:
        assignments[read.query_name] = 'hap1'
    elif (65006192, 65008900) in junctions and (65008954, 65016435) in junctions:
        assignments[read.query_name] = 'hap2'
    else:
        assignments[read.query_name] = 'NA'

qscore = PhasingQuality(workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam', assignments, target[1], target[2])
with open(workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/phasing_stats.txt', 'w') as quality_file:
    _ = quality_file.write('Phasing quality score\t%s\n' % format(qscore, '.4f'))
    _ = quality_file.write('Read count (haplotype 1)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap1'}))
    _ = quality_file.write('Read count (haplotype 2)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap2'}))
    _ = quality_file.write('Read count (unassigned)\t%s\n' % len({key for key in assignments if assignments[key] == 'NA'}))

WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap1'}, workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam', 
    workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/hap1_reads.bam')

WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap2'}, workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam', 
    workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/hap2_reads.bam')

WriteBam(samtools, {key for key in assignments if assignments[key] == 'NA'}, workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam', 
    workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/unassigned_reads.bam')

_ = os.system('rm ' + workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam && rm ' + 
    workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam.bai')

# Our initial round of phasing confirmed that the two splicing defects are in trans. Furthermore, our phasing results revealed
# a missense variant in exon 4, c.626C>T, p.Ser209Leu, (among transcripts that do not include the pseudoexon in intron 6) that 
# is curated as pathogenic and is known to induce exon 4 skipping. Notably, some of the reads that do not include the pseudoexon in 
# intron 6 also show the reference allele at this missense variant position, suggesting that the pseudoexon is not always 
# included in transcripts of MTFMT. In light of this new evidence, we proposed the following second pass "pseudo" phasing
# strategy: 
#   * Haplotype 1: Reads that skip exon 4 or include exon 4 but have the c.626C>T, p.Ser209Leu variant 
#   * Haplotype 2: Reads that include exon 4 but do not have the c.626C>T, p.Ser209Leu variant

import sys

# Establish working directories and file paths
workdir = '/mnt/isilon/lin_lab_share/STRIPE'
samtools = '/scr1/users/wangr5/tools/samtools-1.21/samtools'

# Import functions used in Run_STRIPE.py
sys.path.append(workdir + '/scripts/main')
from Run_STRIPE import *

# Read in mapping file linking patient IDs to their original IDs
sampleMap = pd.read_csv(workdir + '/manuscript/Revisions/20250228/sample.map', sep = '\t', header = None)
sampleMap = dict(zip(sampleMap[1], sampleMap[0]))

patient_id, gene_panel, gene_name = sampleMap['PMD-P23'], 'PMD-359', 'MTFMT'

targetDF = pd.read_csv(workdir + '/' + gene_panel.split('-')[0] + '/references/target_genes.bed', sep = '\t', header = None)
target, assignments = targetDF[targetDF[4] == gene_name].values[0], dict()
_ = os.system(samtools + ' view -hb -F 256 -q 1 ' + workdir + '/' + gene_panel.split('-')[0] +
    '/' + patient_id + '/RNA/' + patient_id + '_TEQUILA.bam' + ' ' + target[0] + ':' + str(target[1]+1) + 
    '-' + str(target[2]) + ' > ' + workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam && ' + 
    samtools + ' index ' + workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam')
samfile = pysam.AlignmentFile(workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam', 'rb')

for read in samfile.fetch():
    blocks = re.findall('(\d+)(\D+)', read.cigarstring)
    cumlen = np.cumsum([0] + [int(item[0]) if item[1] in {'M', 'D', 'N'} else 0 for item in blocks][:-1])
    blocks = [(blocks[i][0], blocks[i][1], cumlen[i] + int(read.reference_start) + 1) for i in range(len(blocks))]
    junctions = [(cumlen[i] + int(read.reference_start) + 1, cumlen[i+1] + int(read.reference_start)) for i in range(len(blocks)) if blocks[i][1] == 'N']
    if (65020273, 65023671) in junctions:
        assignments[read.query_name] = 'hap1'
    else:
        assignments[read.query_name] = 'NA'

for pileupcolumn in pileup_truncated(samfile, 'chr15', 65021532, 65021533):
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            if pileupread.alignment.query_sequence[pileupread.query_position] == 'A':
                if assignments.get(pileupread.alignment.query_name, 'NA') == 'NA':
                    assignments[pileupread.alignment.query_name] = 'hap1'
            elif pileupread.alignment.query_sequence[pileupread.query_position] == 'G':
                if assignments.get(pileupread.alignment.query_name, 'NA') == 'NA':
                    assignments[pileupread.alignment.query_name] = 'hap2'

qscore = PhasingQuality(workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam', assignments, target[1], target[2])
with open(workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/phasing_stats.txt', 'w') as quality_file:
    _ = quality_file.write('Phasing quality score\t%s\n' % format(qscore, '.4f'))
    _ = quality_file.write('Read count (haplotype 1)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap1'}))
    _ = quality_file.write('Read count (haplotype 2)\t%s\n' % len({key for key in assignments if assignments[key] == 'hap2'}))
    _ = quality_file.write('Read count (unassigned)\t%s\n' % len({key for key in assignments if assignments[key] == 'NA'}))

WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap1'}, workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam', 
    workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/hap1_reads.bam')

WriteBam(samtools, {key for key in assignments if assignments[key] == 'hap2'}, workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam', 
    workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/hap2_reads.bam')

WriteBam(samtools, {key for key in assignments if assignments[key] == 'NA'}, workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam', 
    workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/unassigned_reads.bam')

_ = os.system('rm ' + workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam && rm ' + 
    workdir + '/manuscript/Revisions/20250228/Main_Figures/Figure_6/all_reads.bam.bai')

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_6/Figure_6a.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_6/tmp"

# Retrieve annotations for the canonical transcript of MTFMT
awk '{if($5 == "MTFMT"){print $4;}}' "$WORKDIR/PMD/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_6/tmp/input.gtf"

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_6/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_6/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_6/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_6/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved MTFMT splicing patterns in individual PMD-P23
# Only visualize splice junctions supported by at least 100 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_6/tmp/input_bam.txt" \
    -c "chr15:65003056-65029639" -o "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_6/Figure_6a" -M 100 \
    -g "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_6/tmp/input.gtf" --shrink -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_6/tmp"
