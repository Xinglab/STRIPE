#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.01.06
# Figure 4a

# (a) Sashimi plot visualization of haplotype-resolved PIGQ splicing patterns in individual CDG-132-1. Haplotype 1 
# carries a variant (NM_004204.5: c.942+1G>A) that disrupts the donor splice site of intron 4. 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_4/Figure_4a.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_4/tmp"

# Retrieve annotations for the canonical transcript of PIGQ
awk '{if($5 == "PIGQ"){print $4;}}' "$WORKDIR/CDG/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_4/tmp/input.gtf"

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/CDG-132-1/PIGQ/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_4/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/CDG-132-1/PIGQ/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_4/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved PIGQ splicing patterns in individual CDG-132-1
# Only visualize splice junctions supported by at least 50 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_4/tmp/input_bam.txt" \
    -c "chr16:574066-578938" -o "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_4/Figure_4a" -M 50 \
    -g "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_4/tmp/input.gtf" -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_4/tmp"