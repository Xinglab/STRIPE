#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.03
# Supplementary Figure 7a

# (a) Sashimi plot visualization of haplotype-resolved EARS2 splicing patterns in individual PMD-P03. Haplotype 1 carries 
# a 25.9-kb structural deletion that removes a genomic region containing exons 3 to 7. 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S7/Figure_S7a.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S7/tmp"

# Retrieve annotations for the canonical transcript of EARS2
awk '{if($5 == "EARS2"){print $4;}}' "$WORKDIR/PMD/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S7/tmp/input.gtf"

# Pull out original sample ID for PMD-P03
SAMPID=$(grep -w "PMD-P03" "$WORKDIR/manuscript/Revisions/20250228/sample.map" | cut -f1)

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/manuscript/Revisions/20250228/Supplementary_Tables/Table_S6/$SAMPID/EARS2/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S7/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/manuscript/Revisions/20250228/Supplementary_Tables/Table_S6/$SAMPID/EARS2/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S7/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved EARS2 splicing patterns in individual PMD-P03
# Only visualize splice junctions supported by at least 200 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S7/tmp/input_bam.txt" \
    -c "chr16:23523954-23557350" -o "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S7/Figure_S7a" -M 200 \
    -g "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S7/tmp/input.gtf" --shrink -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S7/tmp"