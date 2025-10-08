#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.03
# Supplementary Figure 8a

# (a) Sashimi plot visualization of haplotype-resolved splicing patterns around PMM2 exons 4 to 7 in individual CDG-P10. 
# Haplotype 1 carries a missense variant (NM_000303.3:c.415G>A, p.Glu139Lys) on exon 5 that is known to induce exon 
# skipping by disrupting a splicing enhancer sequence. 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S8/Figure_S8a.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S8/tmp"

# Retrieve annotations for the canonical transcript of PMM2
awk '{if($5 == "PMM2"){print $4;}}' "$WORKDIR/CDG/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S8/tmp/input.gtf"

# Pull out original sample ID for CDG-P10
SAMPID=$(grep -w "CDG-P10" "$WORKDIR/manuscript/Revisions/20250228/sample.map" | cut -f1)

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/manuscript/Revisions/20250228/Supplementary_Tables/Table_S6/$SAMPID/PMM2/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S8/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/manuscript/Revisions/20250228/Supplementary_Tables/Table_S6/$SAMPID/PMM2/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S8/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved splicing patterns around PMM2 exons 4 to 7 in individual CDG-P10
# Only visualize splice junctions supported by at least 400 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S8/tmp/input_bam.txt" \
    -c "chr16:8806316-8813106" -o "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S8/Figure_S8a" -M 400 \
    -g "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S8/tmp/input.gtf" --shrink -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20250228/Supplementary_Figures/Figure_S8/tmp"
