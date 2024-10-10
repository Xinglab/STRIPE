#!/bin/bash

# Author: Robert Wang (Xing Lab)
# Date: 2024.10.10

# Generate a multiple sequence alignment between the canonical transcript sequence for ALG1 and those of its 8 known 
# pseudogenes (https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=56052)

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
WORKDIR="/mnt/isilon/lin_lab_share/STRIPE"
GENOME="/scr1/users/wangr5/references/GRCh38.primary_assembly.genome.fa"
GENCODE="/scr1/users/wangr5/references/gencode.v45.annotation.gtf"
OUTDIR="$WORKDIR/manuscript/Main_Figures/Figure_6"

mkdir -p "$OUTDIR/ALG1"

# Pull out sequences of canonical transcripts for ALG1 and its 8 annotated pseudogenes
GENE_ARR=("ALG1" "TSSC2" "ALG1L6P" "ALG1L3P" "ALG1L5P" "ALG1L9P" "ALG1L12P" "ALG1L7P" "ALG1L8P")
TXID_ARR=("ENST00000262374" "ENST00000450217" "ENST00000515248" "ENST00000503331"
    "ENST00000482043" "ENST00000532875" "ENST00000515046" "ENST00000502317" "ENST00000533887")

for i in ${!GENE_ARR[@]}; do
    gffread -w "$OUTDIR/ALG1/${GENE_ARR[$i]}.fa" -W -g "$GENOME" <(grep "${TXID_ARR[$i]}" "$GENCODE")
done

# Perform a multiple sequence alignment of all sequences using MUSCLE (v5.1)
muscle -align <(printf "%s\n" "${GENE_ARR[@]}" | sed "s|^|$OUTDIR/ALG1/|" | sed "s/$/.fa/" | xargs cat) \
    -output "$OUTDIR/ALG1/sequences.afa"

# Remove new lines from individual sequences in output file from MUSCLE
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' \
    "$OUTDIR/ALG1/sequences.afa" | paste - - | sort -k1,1V | tr '\t' '\n' > "$OUTDIR/ALG1/sequences.fa"
rm "$OUTDIR/ALG1/sequences.afa"
