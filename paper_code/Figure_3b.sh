#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.03
# Figure 3b

# (b) Sashimi plot visualization of ATP5MK splicing patterns in individual PMD-P01, who is homozygous for a variant 
# (NM_001206427.2:c.87+1G>C) that disrupts the donor splice site of intron 3, and an unaffected healthy control (PMD-C01). 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_3/Figure_3b.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp"

# Retrieve annotations for the canonical transcript of ATP5MK
awk '{if($5 == "ATP5MK"){print $4;}}' "$WORKDIR/PMD/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/input.gtf"

# Pull out original sample IDs for PMD-P01 and PMD-C01
SAMPID1=$(grep -w "PMD-P01" "$WORKDIR/manuscript/Revisions/20250228/sample.map" | cut -f1)
SAMPID2=$(grep -w "PMD-C01" "$WORKDIR/manuscript/Revisions/20250228/sample.map" | cut -f1)

# Construct text file containing file paths to TEQUILA-seq BAM files
echo -e "Control\t$WORKDIR/PMD/${SAMPID2}/RNA/${SAMPID2}_TEQUILA.bam\t1" > "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/input_bam.txt"
echo -e "Patient\t$WORKDIR/PMD/${SAMPID1}/RNA/${SAMPID1}_TEQUILA.bam\t2" >> "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/input_bam.txt"

# Use ggsashimi.py to visualize ATP5MK splicing patterns in individual PMD-P01
# Only visualize splice junctions supported by at least 130 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/input_bam.txt" \
    -c "chr10:103389050-103396475" -o "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_3/Figure_3b" -M 130 \
    -g "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp/input.gtf" --shrink -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_3/tmp"