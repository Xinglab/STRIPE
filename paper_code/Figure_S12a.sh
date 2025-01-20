#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.01.09
# Supplementary Figure 12a

# (a) Sashimi plot visualization of haplotype-resolved TRIT1 splicing patterns in individual Q1819. Haplotype 1 carries 
# a VUS (NM_017646.6:c.702A>G, p.Ala234=) that disrupts the -2 position of the intron 5 donor splice site region. 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S12/Figure_S12a.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S12/tmp"

# Retrieve annotations for the canonical transcript of TRIT1
awk '{if($5 == "TRIT1"){print $4;}}' "$WORKDIR/PMD/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S12/tmp/input.gtf"

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/Q1819/TRIT1/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S12/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/Q1819/TRIT1/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S12/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved TRIT1 splicing patterns in individual Q1819
# Only visualize splice junctions supported by at least 110 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S12/tmp/input_bam.txt" \
    -c "chr1:39847548-39883504" -o "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S12/Figure_S12a" -M 110 \
    -g "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S12/tmp/input.gtf" --shrink -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20241219/Supplementary_Figures/Figure_S12/tmp"