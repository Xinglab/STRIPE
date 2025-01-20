#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.01.06
# Supplementary Figure 13a

# (a) Sashimi plot visualization of haplotype-resolved PIGN splicing patterns in individual CDG-147-1. Haplotype 2 carries 
# a synonymous variant (NM_176787.5:c.963G>A, p.Gln321=) on exon 11 that maps to the -1 position of the intron 11 donor 
# splice site region. 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/Figure_S13a.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp"

# Retrieve annotations for the canonical transcript of PIGN
awk '{if($5 == "PIGN"){print $4;}}' "$WORKDIR/CDG/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp/input.gtf"

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/CDG-147-1/PIGN/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/CDG-147-1/PIGN/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved PIGN splicing patterns in individual CDG-147-1
# Only visualize splice junctions supported by at least 25 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp/input_bam.txt" \
    -c "chr18:62140420-62147101" -o "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/Figure_S13a" -M 25 \
    -g "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp/input.gtf" -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S13/tmp"