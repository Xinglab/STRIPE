#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.01.06
# Supplementary Figure 13c,d

# (c) Sashimi plot visualization of haplotype-resolved PIGN splicing patterns in individual AnJa. Haplotype 2 carries a 
# variant (NM_176787.5:c.1434+5G>A) that maps to the +5 position of the intron 16 donor splice site region. (d) Same as 
# in (c) but for haplotype 1 in individual CDG-161-1. 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE1="$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/Figure_S13c.pdf"
OUTFILE2="$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/Figure_S13d.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp"

# Retrieve annotations for the canonical transcript of PIGN
awk '{if($5 == "PIGN"){print $4;}}' "$WORKDIR/CDG/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp/input.gtf"

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/AnJa/PIGN/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/AnJa/PIGN/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved PIGN splicing patterns in individual AnJa
# Only visualize splice junctions supported by at least 50 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp/input_bam.txt" \
    -c "chr18:62109834-62114639" -o "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/Figure_S13c" -M 50 \
    -g "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp/input.gtf" -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/CDG-161-1/PIGN/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/CDG-161-1/PIGN/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved PIGN splicing patterns in individual CDG-161-1
# Only visualize splice junctions supported by at least 50 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp/input_bam.txt" \
    -c "chr18:62109834-62114639" -o "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/Figure_S13d" -M 50 \
    -g "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp/input.gtf" -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp"