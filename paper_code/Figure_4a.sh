#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.05
# Figure 4a

# (a) Sashimi plot visualization of haplotype-resolved PIGQ splicing patterns in individual CDG-132-1. Haplotype 1 carries
# a heterozygous variant (NM_004204.5: c.942+1G>A) that disrupts the donor splice site of intron 4.

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Main_Figures/Figure_4/Figure_4a.pdf"
mkdir -p "$WORKDIR/manuscript/Main_Figures/Figure_4/tmp"

# Retrieve annotations for the canonical transcript of PIGQ
awk '{if($5 == "PIGQ"){print $4;}}' "$WORKDIR/CDG/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Main_Figures/Figure_4/tmp/input.gtf"
REGION=$(awk '{if($3 == "transcript"){printf("%s:%s-%s", $1, $4, $5)}}' "$WORKDIR/manuscript/Main_Figures/Figure_4/tmp/input.gtf")

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "CDG-132-1 (haplotype 1)\t$WORKDIR/CDG/CDG-132-1/RNA/stripe/target_genes/PIGQ/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Main_Figures/Figure_4/tmp/input_bam.txt"
echo -e "CDG-132-1 (haplotype 2)\t$WORKDIR/CDG/CDG-132-1/RNA/stripe/target_genes/PIGQ/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Main_Figures/Figure_4/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved PIGQ splicing patterns in individual CDG-132-1
# Only visualize splice junctions supported by at least 50 reads
(python "/scr1/users/wangr5/tools/ggsashimi.py" -b "$WORKDIR/manuscript/Main_Figures/Figure_4/tmp/input_bam.txt" \
    -c "$REGION" -o "$WORKDIR/manuscript/Main_Figures/Figure_4/Figure_4a" -M 50 -g "$WORKDIR/manuscript/Main_Figures/Figure_4/tmp/input.gtf" \
    -C 3 --alpha 1 --height 1 --ann-height 0.5 --width 4 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Main_Figures/Figure_4/tmp"
