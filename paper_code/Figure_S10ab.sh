#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.05
# Supplementary Figure 10a,b

# (a, b) Sashimi plot visualization of haplotype-resolved PIGN splicing patterns in individuals (a) AnJa and (b) 
# CDG-161-1. Both individuals are heterozygous for a variant (NM_176787.5: c.1434+5G>A) that disrupts the donor 
# splice site region of intron 16.

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE1="$WORKDIR/manuscript/Supplementary_Figures/Figure_S10/Figure_S10a.pdf"
OUTFILE2="$WORKDIR/manuscript/Supplementary_Figures/Figure_S10/Figure_S10b.pdf"
mkdir -p "$WORKDIR/manuscript/Supplementary_Figures/Figure_S10/tmp"

# Retrieve annotations for the canonical transcript of PIGN
awk '{if($5 == "PIGN"){print $4;}}' "$WORKDIR/CDG/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Supplementary_Figures/Figure_S10/tmp/input.gtf"

# Construct text files containing file paths to haplotype-resolved BAM files
echo -e "AnJa (haplotype 1)\t/mnt/isilon/lin_lab_share/TEQUILA-Dx/CDG/rna_seq/AnJa/examples/PIGN/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Supplementary_Figures/Figure_S10/tmp/input_bam1.txt"
echo -e "AnJa (haplotype 2)\t/mnt/isilon/lin_lab_share/TEQUILA-Dx/CDG/rna_seq/AnJa/examples/PIGN/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Supplementary_Figures/Figure_S10/tmp/input_bam1.txt"

echo -e "CDG-161-1 (haplotype 1)\t$WORKDIR/CDG/CDG-161-1/RNA/stripe/target_genes/PIGN/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Supplementary_Figures/Figure_S10/tmp/input_bam2.txt"
echo -e "CDG-161-1 (haplotype 2)\t$WORKDIR/CDG/CDG-161-1/RNA/stripe/target_genes/PIGN/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Supplementary_Figures/Figure_S10/tmp/input_bam2.txt"

# Use ggsashimi.py to visualize haplotype-resolved PIGN splicing patterns in individuals AnJa and CDG-161-1
# Only visualize splice junctions supported by at least 20 reads
(python "/scr1/users/wangr5/tools/ggsashimi.py" -b "$WORKDIR/manuscript/Supplementary_Figures/Figure_S10/tmp/input_bam1.txt" \
    -c "chr18:62105000-62115000" -o "$WORKDIR/manuscript/Supplementary_Figures/Figure_S10/Figure_S10a" -M 20 -g "$WORKDIR/manuscript/Supplementary_Figures/Figure_S10/tmp/input.gtf" \
    -C 3 --alpha 1 --height 1 --ann-height 0.5 --width 3.75 --base-size 6) > /dev/null 2>&1

(python "/scr1/users/wangr5/tools/ggsashimi.py" -b "$WORKDIR/manuscript/Supplementary_Figures/Figure_S10/tmp/input_bam2.txt" \
    -c "chr18:62105000-62115000" -o "$WORKDIR/manuscript/Supplementary_Figures/Figure_S10/Figure_S10b" -M 20 -g "$WORKDIR/manuscript/Supplementary_Figures/Figure_S10/tmp/input.gtf" \
    -C 3 --alpha 1 --height 1 --ann-height 0.5 --width 3.75 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Supplementary_Figures/Figure_S10/tmp"
