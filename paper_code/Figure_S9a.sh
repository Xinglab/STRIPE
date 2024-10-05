#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.05
# Supplementary Figure 9a

# (a) Sashimi plot visualization of haplotype-resolved NUBPL splicing patterns in individual Q1687. Haplotype 1 carries
# a missense variant (NM_025152.3: c.166G>A, p.Gly56Arg) on exon 2.

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Supplementary_Figures/Figure_S9/Figure_S9a.pdf"
mkdir -p "$WORKDIR/manuscript/Supplementary_Figures/Figure_S9/tmp"

# Retrieve annotations for the canonical transcript of NUBPL
awk '{if($5 == "NUBPL"){print $4;}}' "$WORKDIR/PMD/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Supplementary_Figures/Figure_S9/tmp/input.gtf"
REGION=$(awk '{if($3 == "transcript"){printf("%s:%s-%s", $1, $4, $5)}}' "$WORKDIR/manuscript/Supplementary_Figures/Figure_S9/tmp/input.gtf")

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Q1687 (haplotype 1)\t$WORKDIR/PMD/Q1687/RNA/stripe/target_genes/NUBPL/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Supplementary_Figures/Figure_S9/tmp/input_bam.txt"
echo -e "Q1687 (haplotype 2)\t$WORKDIR/PMD/Q1687/RNA/stripe/target_genes/NUBPL/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Supplementary_Figures/Figure_S9/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved NUBPL splicing patterns in individual Q1687
# Only visualize splice junctions supported by at least 80 reads
(python "/scr1/users/wangr5/tools/ggsashimi.py" -b "$WORKDIR/manuscript/Supplementary_Figures/Figure_S9/tmp/input_bam.txt" \
    -c "$REGION" -o "$WORKDIR/manuscript/Supplementary_Figures/Figure_S9/Figure_S9a" -M 80 -g "$WORKDIR/manuscript/Supplementary_Figures/Figure_S9/tmp/input.gtf" \
    --shrink -C 3 --alpha 1 --height 1 --ann-height 0.5 --width 4 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Supplementary_Figures/Figure_S9/tmp"
