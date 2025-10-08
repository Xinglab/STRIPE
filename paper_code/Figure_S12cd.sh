#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.03
# Supplementary Figure 12c,d

# (c) Sashimi plot visualization of haplotype-resolved splicing patterns around PIGN exons 15 to 17 in individual CDG-P02. 
# Haplotype 2 carries a variant (NM_176787.5:c.1434+5G>A) that maps to the +5 position of the intron 16 donor splice site. 
# (d) Same as in (c) but for haplotype 1 in individual CDG-P06. 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE1="$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/Figure_S12c.pdf"
OUTFILE2="$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/Figure_S12d.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp"

# Retrieve annotations for the canonical transcript of PIGN
awk '{if($5 == "PIGN"){print $4;}}' "$WORKDIR/CDG/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input.gtf"

# Pull out original sample ID for CDG-P02
SAMPID=$(grep -w "CDG-P02" "$WORKDIR/manuscript/Revisions/20250228/sample.map" | cut -f1)

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/manuscript/Revisions/20250228/Supplementary_Tables/Table_S6/$SAMPID/PIGN/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/manuscript/Revisions/20250228/Supplementary_Tables/Table_S6/$SAMPID/PIGN/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved PIGN splicing patterns in individual CDG-P02
# Only visualize splice junctions supported by at least 50 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input_bam.txt" \
    -c "chr18:62109834-62114639" -o "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/Figure_S12c" -M 50 \
    -g "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input.gtf" -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

# Pull out original sample ID for CDG-P06
SAMPID=$(grep -w "CDG-P06" "$WORKDIR/manuscript/Revisions/20250228/sample.map" | cut -f1)

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/manuscript/Revisions/20250228/Supplementary_Tables/Table_S6/$SAMPID/PIGN/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/manuscript/Revisions/20250228/Supplementary_Tables/Table_S6/$SAMPID/PIGN/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved splicing patterns around PIGN exons 15 to 17 in individual CDG-P06
# Only visualize splice junctions supported by at least 50 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input_bam.txt" \
    -c "chr18:62109834-62114639" -o "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/Figure_S12d" -M 50 \
    -g "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input.gtf" -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp"
