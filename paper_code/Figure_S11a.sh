#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.03
# Supplementary Figure 11a

# (a) Sashimi plot visualization of haplotype-resolved splicing patterns around TRIT1 exons 1 to 7 in individual PMD-P12. 
# Haplotype 1 carries a synonymous VUS (NM_017646.6:c.702A>G, p.Ala234=) that maps to the -2 position of the intron 5 donor 
# splice site. 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/Figure_S11a.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/tmp"

# Retrieve annotations for the canonical transcript of TRIT1
awk '{if($5 == "TRIT1"){print $4;}}' "$WORKDIR/PMD/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/tmp/input.gtf"

# Pull out original sample ID for PMD-P12
SAMPID=$(grep -w "PMD-P12" "$WORKDIR/manuscript/Revisions/20250228/sample.map" | cut -f1)

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/PMD/$SAMPID/RNA/stripe/target_genes/TRIT1/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/PMD/$SAMPID/RNA/stripe/target_genes/TRIT1/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved splicing patterns around TRIT1 exons 1 to 7 in individual PMD-P12
# Only visualize splice junctions supported by at least 110 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/tmp/input_bam.txt" \
    -c "chr1:39847548-39883504" -o "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/Figure_S11a" -M 110 \
    -g "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/tmp/input.gtf" --shrink -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S11/tmp"
