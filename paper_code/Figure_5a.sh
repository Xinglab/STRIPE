#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.09
# Figure 5a

# Sashimi plot visualization of haplotype-resolved NGLY1 splicing patterns in individual FCDGC-02003. Haplotype 1 (maternal) 
# carries a frameshift deletion variant on exon 10 (NM_018297.4: c.1533_1536del, p.Asn511LysfsTer51). 

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Main_Figures/Figure_5/Figure_5a.pdf"
mkdir -p "$WORKDIR/manuscript/Main_Figures/Figure_5/tmp"

# Retrieve annotations for the canonical transcript of NGLY1
awk '{if($5 == "NGLY1"){print $4;}}' "$WORKDIR/CDG/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Main_Figures/Figure_5/tmp/input.gtf"
REGION=$(awk '{if($3 == "transcript"){printf("%s:%s-%s", $1, $4, $5)}}' "$WORKDIR/manuscript/Main_Figures/Figure_5/tmp/input.gtf")

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "FCDGC-02003 (haplotype 1)\t$WORKDIR/CDG/FCDGC-02003/RNA/stripe/target_genes/NGLY1/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Main_Figures/Figure_5/tmp/input_bam.txt"
echo -e "FCDGC-02003 (haplotype 2)\t$WORKDIR/CDG/FCDGC-02003/RNA/stripe/target_genes/NGLY1/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Main_Figures/Figure_5/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved NGLY1 splicing patterns in individual FCDGC-02003
# Only visualize splice junctions supported by at least 250 reads
(python "/scr1/users/wangr5/tools/ggsashimi.py" -b "$WORKDIR/manuscript/Main_Figures/Figure_5/tmp/input_bam.txt" \
    -c "$REGION" -o "$WORKDIR/manuscript/Main_Figures/Figure_5/Figure_5a" -M 250 -g "$WORKDIR/manuscript/Main_Figures/Figure_5/tmp/input.gtf" \
    --shrink -C 3 --alpha 1 --height 1 --ann-height 0.5 --width 4 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Main_Figures/Figure_5/tmp"
