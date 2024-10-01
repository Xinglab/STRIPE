## STRIPE (Sequencing Targeted RNAs Identifies Pathogenic Events)

STRIPE is a versatile workflow that enables deep sequencing of full-length transcripts for any disease-specific gene panel such that a wide range of clinically informative readouts, including transcript aberrations and sequence variants, can be detected at haplotype-level resolution. 

### Dependencies

**Tools:**
* Bcftools (version 1.21)
* CADD (version 1.6)
* LongcallR (version 0.1.0)
* Minimap2 (version 2.26)
* Python (version 3.8.18)
    * numpy (version 1.24.2)
    * pandas (version 2.0.3)
    * pysam (version 0.21.0)
    * pytabix (version 0.1)
    * tensorflow (version 2.13.1)
* R (version 4.3.1)
    * argparse (version 2.2.2)
    * dplyr (version 1.1.4)
    * ggplot2 (version 3.4.4)
    * tidyr (version 1.3.1)
* Samtools (version 1.21)
* Singularity (version 3.11.1)
    * Clair3 (version 1.0.10)
    * DeepVariant (version 1.6.1)
    * DeepTrio (version 1.6.1)
* SpliceAI (version 1.3.1)
* Stringtie (version 2.2.3)

**Databases:**
* GTEx v8:
    * [Gene TPM Matrix](https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz)
    * [Splice Junction Read Count Matrix](https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz)
    * [Haplotype Expression Matrix](https://storage.googleapis.com/adult-gtex/haplotype-expression/v8/haplotype-expression-matrices/phASER_GTEx_v8_matrix.txt.gz)
    * [Sample Annotations](https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt)
* gnomAD v4.1.0:
    * Please download the VCF files (chr1-22,X,Y) under "Genomes" [here](https://gnomad.broadinstitute.org/downloads).
* ClinVar:
    * Please download the latest VCF file of ClinVar variants [here](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/).

### Usage:

**Building a STRIPE database**

```
python Build_STRIPE.py [-h] --gene_list /path/to/target/gene/list --anno_gtf /path/to/gene/annotation/gtf
    [--gtex_files /path/to/GTEx/files/list] [--tissue_list /path/to/GTEx/tissues/file]
    [--gnomad_files /path/to/gnomAD/files/list] [--only_bed] --outdir /path/to/output/folder
```

|Parameter|Description|
|-|-|
|`--gene_list`|Path to a single-column text file with a list of Ensembl gene IDs for target genes|
|`--anno_gtf`|Path to a GTF file with reference gene annotations from GENCODE (download [here](https://www.gencodegenes.org/human/releases.html))|
|`--gtex_files`|(Optional) Path to a single-column text file containing file paths to GTEx v8 files (see Databases)|
|`--tissue_list`|(Optional) Path to a two-column text file in which column 1 has the official names of GTEx tissues (see [here](https://gtexportal.org/home/tissueSummaryPage)) and column 2 has the desired labels for each tissue|
|`--gnomad_files`|(Optional) Path to a single-column text file containing file paths to all gnomAD VCF files (chr1-22,X,Y)|
|`--only_bed`|(Optional) Forces `Build_STRIPE.py` to only construct a BED file for input target genes|
|`--outdir`|Path to output directory for STRIPE database files|

`Build_STRIPE.py` will generate the following intermediate folders/files in `--outdir` that are needed for running STRIPE:
* `target_genes.bed`: A BED6 file containing the genomic coordinates of input target genes (based on reference gene annotations from `--anno_gtf`)
* `GTEx_v8`: A folder containing subfolders for each tissue label in `--tissue_list`. Each subfolder contains the following files:
    * `target_genes.gene_tpm.txt`: A matrix file containing normalized abundances (TPM) for target genes across GTEx samples of the given tissue type
    * `target_genes.junction_counts.txt`: A matrix file containing read counts for splice junctions mapping to target genes across GTEx samples of the given tissue type
    * `target_genes.haplotype_expression.txt`: A matrix file containing haplotype-specific read counts for target genes across GTEx samples of the given tissue type
* `gnomAD_v4`: A folder containing subfolders for each target gene. Each subfolder contains a file `gnomad_variants.tsv`, which describes all known gnomAD variants mapping to the target gene as follows:
    * `CHROM`, `POS`, `REF`, `ALT`: Chromosome, position, reference allele, alternate allele
    * `AF_GRPMAX_JOINT`: Maximum allele frequency across gnomAD populations
    * `NHOMALT_GRPMAX_JOINT`: Maximum number of homozygous individuals across gnomAD populations

**Running STRIPE**

```
python Run_STRIPE.py [-h] --infile /path/to/input/TEQUILA/bam --targets /path/to/target/gene/BED --genome /path/to/reference/genome/FASTA
    --gencode /path/to/gencode/GTF --samtools /path/to/samtools --stringtie /path/to/stringtie --longcallr /path/to/longcallr
    --clair3 /path/to/clair3/sif --model /path/to/clair3/model --cadd /path/to/CADD --gnomad /path/to/gnomad/database
    --clinvar /path/to/clinvar/vcf --bcftools /path/to/bcftools --minimap2 /path/to/minimap2 --gtex_splice /path/to/gtex/splice/database
    --gtex_haplotype /path/to/gtex/haplotype/database [--threads <number of threads>] [--proband_vcf /path/to/proband/vcf]
    [--parent1_vcf /path/to/parent1/vcf] [--parent2_vcf /path/to/parent2/vcf] [--only_qc] [--phase_indels] --outdir /path/to/output/directory
```

|Parameter|Description|
|-|-|
|`--infile`|Path to BAM file with TEQUILA-seq read alignments|
|`--targets`|Path to BED6 file containing genomic coordinates of input target genes (output of running `Build_STRIPE.py`)|
|`--genome`|Path to FASTA file with reference genome sequence|
|`--gencode`|Path to GTF file with reference gene annotations from GENCODE (please use the same file as the one used in running `Build_STRIPE.py`)|
|`--samtools`|Path to `samtools` executable (see [here](https://github.com/samtools/samtools/releases))|
|`--stringtie`|Path to `stringtie` executable (see [here](https://github.com/gpertea/stringtie/releases))|
|`--longcallr`|Path to `longcallr` executable (see [here](https://github.com/huangnengCSU/longcallR/releases))|
|`--clair3`|Path to SIF file for Clair3. Run `singularity pull` on the appropriate SIF file from [here](https://hub.docker.com/r/hkubal/clair3/tags)|
|`--model`|Path to folder for Clair3 model. Clair3 models based on the latest nanopore sequencing chemistries can be downloaded from [here](https://github.com/nanoporetech/rerio/tree/master/clair3_models)|
|`--cadd`|Path to `CADD.sh` script (download and follow installation instructions [here](https://github.com/kircherlab/CADD-scripts/archive/v1.6.post1.zip))|
|`--gnomad`|Path to `gnomAD_v4` output folder generated by `Build_STRIPE.py`|
|`--clinvar`|Path to VCF file containing the latest release of ClinVar variants (see [here](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/))|
|`--bcftools`|Path to `bcftools` executable (see [here](https://github.com/samtools/bcftools/releases))|
|`--minimap2`|Path to `minimap2` executable (see [here](https://github.com/lh3/minimap2/releases))|
|`--gtex_splice`|Path to the appropriate `target_genes.junction_counts.txt` file generated by `Build_STRIPE.py`|
|`--gtex_haplotype`|Path to the appropriate `target_genes.haplotype_expression.txt` file generated by `Build_STRIPE.py`|
|`--threads`|(Optional) Number of threads (default: 1)|
|`--proband_vcf`|(Optional) Path to a VCF file with small variants called from prior exome or genome sequencing data for the proband|
|`--parent1_vcf`|(Optional) Path to a VCF file with small variants called from prior exome or genome sequencing data for the proband's parent|
|`--parent2_vcf`|(Optional) Path to a VCF file with small variants called from prior exome or genome sequencing data for the proband's parent|
|`--only_qc`|(Optional) Forces `Run_STRIPE.py` to only perform a quality control check on the input BAM file (see below for more details)|
|`--phase_indels`|(Optional) Allows `Run_STRIPE.py` to use indels called from prior exome or genome sequencing data to inform phasing of TEQUILA-seq reads|
|`--outdir`|Path to output directory for STRIPE|

`Run_STRIPE.py` consists of the following steps:
1. **Quality control analyses:** Generates the following files:
    * `[SAMPLE].mapping_stats.txt`: Read mapping statistics and the overall on-target rate for `SAMPLE`
    * `[SAMPLE].top_2000_genes.txt`: A three-column text file describing the 2000 most abundant genes in `SAMPLE` as follows:
        * Column 1: Gene ID
        * Column 2: Indicator for whether the given gene is a target gene or not
        * Column 3: Estimated abundance for the gene (TPM)
    * `[SAMPLE].top_2000_genes.pdf`: Plot showing the abundances (TPM) of the 2000 most abundant genes in `SAMPLE`, in which each gene is colored based on whether it is targeted or not. The overall on-target rate for `SAMPLE` is also indicated on the plot.
    * `[SAMPLE].junctions.tsv`: A four-column text file describing splice junctions discovered across all target genes in `SAMPLE` as follows:
        * Columns 1-3: (Chrom, Start, End) for each discovered splice junction
        * Column 4: Read count for each discovered splice junction
    * `stringtie/[SAMPLE].gtf`: A GTF file describing all genes and transcripts discovered by `stringtie` in `SAMPLE`
    * `stringtie/[SAMPLE].gene_abundance.tsv`: A tab-separated file describing all genes discovered by `stringtie` in `SAMPLE` and their abundances

2. **Phase TEQUILA-seq reads mapping to target genes:** For each target gene, STRIPE will generate the following files (if TEQUILA-seq reads mapping to the gene can be phased)
    * `gene.vcf.gz`: A compressed and indexed VCF file containing phased variants used for guiding assignment of TEQUILA-seq reads to individual haplotypes
    * `hap1_reads.bam`: A BAM file containing alignments of TEQUILA-seq reads assigned to haplotype 1
    * `hap2_reads.bam`: A BAM file containing alignments of TEQUILA-seq reads assigned to haplotype 2
    * `phasing_stats.txt`: A text file summarizing the results of phasing TEQUILA-seq reads, including:
        * Number of reads assigned to each haplotype and the number of reads that could not be assigned to a haplotype
        * A phasing "quality score", which is computed as follows. Suppose we have a total of $N$ reads mapping to some gene with start and end genomic coordinates $(S, T)$. The quality score, $Q$, is computed as:
            $$Q = \frac{\sum_{k=1}^{N}w_kz_k}{\sum_{k=1}^{N}w_k}$$
            where $w_k = \max(\min(t_k, T) - \max(s_k, S), 0)/(T-S)$ and $(s_k, t_k)$ denotes the start and end genomic coordinates for the given read alignment. Here, a higher value of $Q$ indicates that more full-length reads mapping to the gene were able to be phased. 

3. **Identify rare, deleterious variants in each target gene:** For each target gene, STRIPE will generate a merged callset of variants from `longcallR`, `clair3`, `deepvariant`, and `bcftools`. The merged callset will then be filtered for variants meeting the following criteria: (i) CADD score > 15 or SpliceAI score > 0.1, (ii) maximum allele frequency across gnomAD populations < 1%. The filtered merged callset, `merged_variants.tsv`, will contain the following fields:
    * `CHROM`, `POS`, `REF`, `ALT`: Chromosome, position, reference allele, alternate allele
    * `LONGCALLR_DP`, `LONGCALLR_AF`, `LONGCALLR_GT`: Read depth, allele frequency, and genotype from `longcallR` for given variant
    * `CLAIR3_DP`, `CLAIR3_AF`, `CLAIR3_GT`: Read depth, allele frequency, and genotype from `clair3` for given variant
    * `DEEPVARIANT_DP`, `DEEPVARIANT_AF`, `DEEPVARIANT_GT`: Read depth, allele frequency, and genotype from `deepvariant` for given variant
    * `PHASING_DP`, `PHASING_AF`, `PHASING_GT`: Read depth, allele frequency, and genotype from `bcftools` for given variant
        * Here, we ran `bcftools` on haplotype-assigned TEQUILA-seq read sets and constructed a phased genotype
    * `CADD`: Functional deleteriousness score from CADD (verison 1.6)
    * `SPLICEAI`: Likelihood of impacting splicing from SpliceAI
    * `GNOMAD_AF_GRPMAX_JOINT`: Maximum allele frequency across gnomAD populations
    * `GNOMAD_NHOMALT_GRPMAX_JOINT`: Maximum number of homozygous individuals across gnomAD populations
    * `CLINVAR`: ClinVar annotations for variant
    * `HOMOPOLYMER`: Indicator for whether the variant lies in a homopolymer region (useful for flagging false positive variant calls due to basecalling artifacts)
    * `END_DISTANCE_PCT`: Percentage of reads supporting the variant in which the variant is found within 30 bp of the alignment ends (useful for flagging false positive variant calls due to alignment artifacts)
    * `SKIP_FREQ`: Percentage of reads covering the variant locus in which the variant position is intronic (useful for flagging false positive variant calls due to local mis-alignments around splice sites)

4. **Detect aberrant splice junctions in each target gene:** For each target gene, STRIPE will extract unique splice junctions from all TEQUILA-seq reads (as well as haplotype-specific reads) mapping to the gene locus and evaluate whether usage frequencies for these splice junctions are unusually higher than expected compared to usage frequencies observed among tissue-matched GTEx controls. Each output file will contain the following fields:
    * `Junction`: Splice junction coordinates
    * `Count`: Number of reads supporting the splice junction
    * `Coverage`: Total number of reads covering the corresponding splice sites
    * `Usage`: Usage frequency of the splice junction in sample
    * `Population`: Average usage frequency of the splice junction across GTEx controls
    * `Shift`: Difference in usage frequencies between the sample and GTEx controls
    * `PVal`: P value for testing whether the sample has a usage frequency that is consistent with those observed in GTEx controls

5. **Find target genes showing strong haplotype dosage imbalance:** For each target gene, STRIPE will also assess whether the gene shows unusually strong haplotype dosage imbalance relative to haplotype expression ratios observed among tissue-matched GTEx controls. The output file will contain the following information:
    * Number of reads assigned to haplotypes 1 and 2
    * Quality score for read phasing results 
    * P value for testing whether the haplotype expression ratio observed in the sample is consistent with those observed in GTEx controls
