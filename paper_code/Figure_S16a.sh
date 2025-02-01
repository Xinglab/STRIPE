#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.01.18
# Supplementary Figure 16a

# Sashimi plot visualization of haplotype-resolved NUBPL splicing patterns in individual Q1687. Haplotype 2 carries a 
# variant (NM_025152.3:c.693+1G>A) that disrupts the donor splice site of intron 8 and activates a downstream cryptic 
# polyadenylation site. 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/Figure_S16a.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/tmp"

# Retrieve annotations for the canonical transcript of NUBPL
awk '{if($5 == "NUBPL"){print $4;}}' "$WORKDIR/PMD/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/tmp/input.gtf"

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/Q1687/NUBPL/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/Q1687/NUBPL/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved NUBPL splicing patterns in individual Q1687
# Only visualize splice junctions supported by at least 80 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/tmp/input_bam.txt" \
    -c "chr14:31561404-31859418" -o "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/Figure_S16a" -M 80 \
    -g "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/tmp/input.gtf" --shrink -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S16/tmp"
