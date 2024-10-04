#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.04
# Figure 2c

# (c) Sashimi plot visualization of haplotype-resolved EARS2 splicing patterns in individual BS2-1. Haplotype 1 carries 
# a 25.9-kb structural deletion that removes a genomic region containing exons 3 to 7. 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Main_Figures/Figure_2/Figure_2c.pdf"
mkdir -p "$WORKDIR/manuscript/Main_Figures/Figure_2/tmp"

# Retrieve annotations for the canonical transcript of EARS2
awk '{if($5 == "EARS2"){print $4;}}' "$WORKDIR/PMD/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Main_Figures/Figure_2/tmp/input.gtf"
REGION=$(awk '{if($3 == "transcript"){printf("%s:%s-%s", $1, $4, $5)}}' "$WORKDIR/manuscript/Main_Figures/Figure_2/tmp/input.gtf")

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "BS2-1 (haplotype 1)\t$WORKDIR/PMD/BS2-1/RNA/stripe/target_genes/EARS2/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Main_Figures/Figure_2/tmp/input_bam.txt"
echo -e "BS2-1 (haplotype 2)\t$WORKDIR/PMD/BS2-1/RNA/stripe/target_genes/EARS2/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Main_Figures/Figure_2/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved EARS2 splicing patterns in individual BS2-1
# Only visualize splice junctions supported by at least 200 reads
(python "/scr1/users/wangr5/tools/ggsashimi.py" -b "$WORKDIR/manuscript/Main_Figures/Figure_2/tmp/input_bam.txt" \
    -c "$REGION" -o "$WORKDIR/manuscript/Main_Figures/Figure_2/Figure_2c" -M 200 -g "$WORKDIR/manuscript/Main_Figures/Figure_2/tmp/input.gtf" \
    --shrink -C 3 --alpha 1 --height 1 --ann-height 0.5 --width 4 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Main_Figures/Figure_2/tmp"
