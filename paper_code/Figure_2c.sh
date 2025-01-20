#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.01.04
# Figure 2c

# (c) Sashimi plot visualization of haplotype-resolved PMM2 splicing patterns in individual CDG-152-1. Haplotype 2 carries 
# a 15.4-kb structural deletion that removes a genomic region containing exons 3 to 7. 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_2/Figure_2c.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_2/tmp"

# Retrieve annotations for the canonical transcript of PMM2
awk '{if($5 == "PMM2"){print $4;}}' "$WORKDIR/CDG/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_2/tmp/input.gtf"

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/CDG-152-1/PMM2/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_2/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/CDG-152-1/PMM2/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_2/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved PMM2 splicing patterns in individual CDG-152-1
# Only visualize splice junctions supported by at least 1000 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_2/tmp/input_bam.txt" \
    -c "chr16:8797839-8849325" -o "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_2/Figure_2c" -M 1000 \
    -g "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_2/tmp/input.gtf" --shrink -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_2/tmp"