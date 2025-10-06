#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.03.04
# Figure 5a

# (a) Sashimi plot visualization of haplotype-resolved splicing patterns around NGLY1 exons 9 to 12 in individual CDG-P16. 
# Haplotype 1 (maternal) carries a frameshift deletion variant on exon 10 (NM_018297.4:c.1533_1536del, p.Asn511LysfsTer51).

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_5/Figure_5a.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_5/tmp"

# Retrieve annotations for the canonical transcript of NGLY1
awk '{if($5 == "NGLY1"){print $4;}}' "$WORKDIR/CDG/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_5/tmp/input.gtf"

# Pull out original sample ID for CDG-P16
SAMPID=$(grep -w "CDG-P16" "$WORKDIR/manuscript/Revisions/20250228/sample.map" | cut -f1)

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/CDG/$SAMPID/RNA/stripe/target_genes/NGLY1/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_5/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/CDG/$SAMPID/RNA/stripe/target_genes/NGLY1/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_5/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved splicing patterns around NGLY1 exons 9 to 12 in individual CDG-P16
# Only visualize splice junctions supported by at least 200 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_5/tmp/input_bam.txt" \
    -c "chr3:25,718,944-25732483" -o "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_5/Figure_5a" -M 200 \
    -g "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_5/tmp/input.gtf" --shrink -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20250228/Main_Figures/Figure_5/tmp"


