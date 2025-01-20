#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2025.01.17
# Figure 5a

# Sashimi plot visualization of haplotype-resolved NGLY1 splicing patterns in individual FCDGC-02003. Haplotype 1 (maternal) 
# carries a frameshift deletion variant on exon 10 (NM_018297.4:c.1533_1536del, p.Asn511LysfsTer51).

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
OUTFILE="$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_5/Figure_5a.pdf"
mkdir -p "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_5/tmp"

# Retrieve annotations for the canonical transcript of NGLY1
awk '{if($5 == "NGLY1"){print $4;}}' "$WORKDIR/CDG/references/target_genes.bed" | \
    grep -f - "/scr1/users/wangr5/references/gencode.v45.annotation.gtf" | grep "Ensembl_canonical" \
    > "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_5/tmp/input.gtf"

# Construct text file containing file paths to haplotype-resolved BAM files
echo -e "Haplotype 1\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/FCDGC-02003/NGLY1/hap1_reads.bam\t1" \
    > "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_5/tmp/input_bam.txt"
echo -e "Haplotype 2\t$WORKDIR/manuscript/Revisions/20241219/Supplementary_Tables/Table_S6/FCDGC-02003/NGLY1/hap2_reads.bam\t2" \
    >> "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_5/tmp/input_bam.txt"

# Use ggsashimi.py to visualize haplotype-resolved NGLY1 splicing patterns in individual FCDGC-02003
# Only visualize splice junctions supported by at least 200 reads
(singularity run -B "$WORKDIR" docker://guigolab/ggsashimi:v1.1.5 -b "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_5/tmp/input_bam.txt" \
    -c "chr3:25,718,944-25732483" -o "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_5/Figure_5a" -M 200 \
    -g "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_5/tmp/input.gtf" --shrink -C 3 --alpha 1 --height 1 \
    --ann-height 0.5 --width 3.5 --base-size 6) > /dev/null 2>&1

rm -rf "$WORKDIR/manuscript/Revisions/20241219/Main_Figures/Figure_5/tmp"
