#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.03
# Supplementary Figure 10a

# (a) Sashimi plot visualization of haplotype-resolved SUCLG1 splicing patterns in individual PMD-P11. Haplotype 2 carries 
# a VUS (NM_003849.4:c.98-9A>G) located 9 bp upstream of the intron 1 acceptor splice site.

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S10/Figure_S10a.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S10/tmp"

# Retrieve annotations for the canonical transcript of SUCLG1
awk '{if($5 == "SUCLG1"){print $4;}}' "$WORKDIR/PMD/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S10/tmp/input.gtf"

# Pull out original sample ID for PMD-P11
SAMPID=$(grep -w "PMD-P11" "$WORKDIR/manuscript/Revisions/20250228/sample.map" | cut -f1)

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/PMD/$SAMPID/RNA/stripe/target_genes/SUCLG1/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S10/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/PMD/$SAMPID/RNA/stripe/target_genes/SUCLG1/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S10/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved SUCLG1 splicing patterns in individual PMD-P11
# Only visualize splice junctions supported by at least 200 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S10/tmp/input_bam.txt" \
    -c "chr2:84443284-84459280" -o "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S10/Figure_S10a" -M 200 \
    -g "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S10/tmp/input.gtf" --shrink -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S10/tmp"