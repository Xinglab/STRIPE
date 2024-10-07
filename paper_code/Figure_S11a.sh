#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.07
# Supplementary Figure 11a

# (a) Sashimi plot visualization of haplotype-resolved PIGN splicing patterns in individual CDG-147-1. Haplotype 1 
# carries a synonymous variant (NM_176787.5: c.963G>A, p.Gln321=) that disrupts the donor splice site region of intron 11. 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Supplementary_Figures/Figure_S11/Figure_S11a.pdf"
mkdir -p "$WORKDIR/manuscript/Supplementary_Figures/Figure_S11/tmp"

# Retrieve annotations for the canonical transcript of PIGN
awk '{if($5 == "PIGN"){print $4;}}' "$WORKDIR/CDG/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Supplementary_Figures/Figure_S11/tmp/input.gtf"

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "CDG-147-1 (haplotype 1)\t$WORKDIR/CDG/CDG-147-1/RNA/stripe/target_genes/PIGN/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Supplementary_Figures/Figure_S11/tmp/input_bam.txt"
echo -e "CDG-147-1 (haplotype 2)\t$WORKDIR/CDG/CDG-147-1/RNA/stripe/target_genes/PIGN/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Supplementary_Figures/Figure_S11/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved PIGN splicing patterns in individual CDG-147-1
# Only visualize splice junctions supported by at least 50 reads
(python "/scr1/users/wangr5/tools/ggsashimi.py" -b "$WORKDIR/manuscript/Supplementary_Figures/Figure_S11/tmp/input_bam.txt" \
    -c "chr18:62138000-62148000" -o "$WORKDIR/manuscript/Supplementary_Figures/Figure_S11/Figure_S11a" -M 50 -g "$WORKDIR/manuscript/Supplementary_Figures/Figure_S11/tmp/input.gtf" \
    -C 3 --alpha 1 --height 1 --ann-height 0.5 --width 3.75 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Supplementary_Figures/Figure_S11/tmp"
