#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.01.06
# Figure 3b

# (b) Sashimi plot visualization of ATP5MK splicing patterns in individual Q1663, who is homozygous for a mutation 
# (NM_001206427.2: c.87+1G>C) that disrupts the donor splice site of intron 3, and an unaffected healthy control (E1877).

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_3/Figure_3b.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_3/tmp"

# Retrieve annotations for the canonical transcript of ATP5MK
awk '{if($5 == "ATP5MK"){print $4;}}' "$WORKDIR/PMD/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_3/tmp/input.gtf"

# Construct text file containing file paths to TEQUILA-seq BAM files
echo -e "Control\t$WORKDIR/PMD/E1877/RNA/E1877_TEQUILA.bam\t1" > "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_3/tmp/input_bam.txt"
echo -e "Patient\t$WORKDIR/PMD/Q1663/RNA/Q1663_TEQUILA.bam\t2" >> "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_3/tmp/input_bam.txt"

# Use ggsashimi.py to visualize ATP5MK splicing patterns in individual Q1663
# Only visualize splice junctions supported by at least 130 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_3/tmp/input_bam.txt" \
    -c "chr10:103389050-103396475" -o "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_3/Figure_3b" -M 130 \
    -g "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_3/tmp/input.gtf" --shrink -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_3/tmp"