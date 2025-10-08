#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.03
# Supplementary Figure 12b

# (b) Sashimi plot visualization of haplotype-resolved splicing patterns around NUBPL exons 7 to 9 in individual PMD-P07. 
# Haplotype 2 carries a variant (NM_025152.3:c.693+1G>A) that disrupts the donor splice site of intron 8. 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/Figure_S12b.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp"

# Retrieve annotations for the canonical transcript of NUBPL
awk '{if($5 == "NUBPL"){print $4;}}' "$WORKDIR/PMD/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input.gtf"

# Pull out original sample ID for PMD-P07
SAMPID=$(grep -w "PMD-P07" "$WORKDIR/manuscript/Revisions/20250228/sample.map" | cut -f1)

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/manuscript/Revisions/20250228/Supplementary_Tables/Table_S6/$SAMPID/NUBPL/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/manuscript/Revisions/20250228/Supplementary_Tables/Table_S6/$SAMPID/NUBPL/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved splicing patterns around NUBPL exons 7 to 9 in individual PMD-P07
# Only visualize splice junctions supported by at least 50 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input_bam.txt" \
    -c "chr14:31787780-31846591" -o "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/Figure_S12b" -M 50 \
    -g "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp/input.gtf" -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S12/tmp"
