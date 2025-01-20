#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.01.09
# Supplementary Figure 10a

# (a) Sashimi plot visualization of haplotype-resolved OPA1 splicing patterns in individual Q2032s1. Haplotype 2 carries 
# a VUS (NM_130837.3:c.1377+5del) that disrupts the +5 position of the intron 14 donor splice site region.

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S10/Figure_S10a.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S10/tmp"

# Retrieve annotations for the canonical transcript of OPA1
awk '{if($5 == "OPA1"){print $4;}}' "$WORKDIR/PMD/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S10/tmp/input.gtf"

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/PMD/Q2032s1/RNA/stripe/target_genes/OPA1/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S10/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/PMD/Q2032s1/RNA/stripe/target_genes/OPA1/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S10/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved OPA1 splicing patterns in individual Q2032s1
# Only visualize splice junctions supported by at least 200 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S10/tmp/input_bam.txt" \
    -c "chr3:193642765-193643627" -o "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S10/Figure_S10a" -M 200 \
    -g "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S10/tmp/input.gtf" --shrink -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S10/tmp"