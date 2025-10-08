#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.03
# Supplementary Figure 12a

# (a) Sashimi plot visualization of haplotype-resolved splicing patterns around PIGN exons 9 to 12 in individual CDG-P05. 
# Haplotype 2 carries a synonymous variant (NM_176787.5:c.963G>A, p.Gln321=) on exon 11 that maps to the -1 position of 
# the intron 11 donor splice site. 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/Figure_S12a.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp"

# Retrieve annotations for the canonical transcript of PIGN
awk '{if($5 == "PIGN"){print $4;}}' "$WORKDIR/CDG/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input.gtf"

# Pull out original sample ID for CDG-P05
SAMPID=$(grep -w "CDG-P05" "$WORKDIR/manuscript/Revisions/20250228/sample.map" | cut -f1)

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/manuscript/Revisions/20250228/Supplementary_Tables/Table_S6/$SAMPID/PIGN/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/manuscript/Revisions/20250228/Supplementary_Tables/Table_S6/$SAMPID/PIGN/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved splicing patterns around PIGN exons 9 to 12 in individual CDG-P05
# Only visualize splice junctions supported by at least 25 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input_bam.txt" \
    -c "chr18:62140420-62147101" -o "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/Figure_S12a" -M 25 \
    -g "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input.gtf" -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp"
