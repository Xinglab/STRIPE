#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.05
# Figure 3b

# (b) Sashimi plot visualization of ATP5MK splicing patterns in individual Q1663, who carries a homozygous mutation 
# (NM_001206427.2: c.87+1G>C) that disrupts the donor splice site of intron 3.

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Main_Figures/Figure_3/Figure_3b.pdf"
mkdir -p "$WORKDIR/manuscript/Main_Figures/Figure_3/tmp"

# Retrieve annotations for the canonical transcript of ATP5MK
awk '{if($5 == "ATP5MK"){print $4;}}' "$WORKDIR/PMD/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Main_Figures/Figure_3/tmp/input.gtf"
REGION=$(awk '{if($3 == "transcript"){printf("%s:%s-%s", $1, $4, $5)}}' "$WORKDIR/manuscript/Main_Figures/Figure_3/tmp/input.gtf")

# Construct text file containing file path to TEQUILA-seq BAM file
echo -e "Q1663\t$WORKDIR/PMD/Q1663/RNA/Q1663_TEQUILA.bam\t1" > "$WORKDIR/manuscript/Main_Figures/Figure_3/tmp/input_bam.txt"

# Use ggsashimi.py to visualize ATP5MK splicing patterns in individual Q1663
# Only visualize splice junctions supported by at least 100 reads
(python "/scr1/users/wangr5/tools/ggsashimi.py" -b "$WORKDIR/manuscript/Main_Figures/Figure_3/tmp/input_bam.txt" \
    -c "$REGION" -o "$WORKDIR/manuscript/Main_Figures/Figure_3/Figure_3b" -M 100 -g "$WORKDIR/manuscript/Main_Figures/Figure_3/tmp/input.gtf" \
    --shrink -C 3 --alpha 1 --height 1 --ann-height 0.5 --width 4 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Main_Figures/Figure_3/tmp"
