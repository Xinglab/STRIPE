#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2024.09.22
Version: 1.0.0

A script designed to prepare intermediate files that are needed when running STRIPE.
'''

# =====================================================================================================================
#                                                       PACKAGES 
# =====================================================================================================================

# Load required packages
import sys, argparse, os
import pandas as pd
import tabix

# =====================================================================================================================
#                                                   HELPER FUNCTIONS
# =====================================================================================================================

def PullFeature(infoString, featureName, format, suffix = True):
    '''
    Function designed to pull out the value for featureName in infoString.
    Version number on value will be removed if suffix is set to False.
    '''
    if format == 'GENCODE':
        return next(iter([item.split('"')[1] if suffix else item.split('"')[1].split('.')[0] 
            for item in infoString.split(';') if featureName in item]), 'NA')
    elif format == 'gnomAD':
        return next(iter([item.split('=')[1] for item in infoString.split(';') if featureName in item]), '0')
    else:
        return 'NA'

# =====================================================================================================================
#                                                    MAIN FUNCTIONS
# =====================================================================================================================

def MakeGeneBed(genelist, annotations, outdir):
    # Filter annotations for genes in genelist
    _ = os.system('awk \'$3 == "gene"\' ' + annotations + ' | grep -f <(cut -d \'.\' -f1 ' + genelist + 
        ') - > ' + outdir + '/tmp.gtf')
    
    # Read in filtered GTF file as a Pandas dataframe and pull out gene name and ID (w/o suffix)
    annoDF = pd.read_csv(outdir + '/tmp.gtf', sep = '\t', header = None)
    annoDF['Gene_Name'] = annoDF[8].apply(PullFeature, featureName = 'gene_name', format = 'GENCODE')
    annoDF['Gene_ID'] = annoDF[8].apply(PullFeature, featureName = 'gene_id', format = 'GENCODE', suffix = False)

    # Determine which genes in genelist are not represented in annoDF
    genes = [item.split('.')[0] for item in pd.read_csv(genelist, header = None)[0]]
    missing = sorted(list(set(genes) - set(annoDF['Gene_ID'])))
    if len(missing) > 0:
        with open(outdir + '/missing_genes.txt', 'w') as f:
            for item in missing:
                f.write('%s\n' % item)

    # Format annoDF as a BED6 object, then save into outdir
    annoDF[3] = annoDF[3] - 1
    annoDF = annoDF[[0, 3, 4, 'Gene_ID', 'Gene_Name', 6]]
    annoDF.columns = range(annoDF.shape[1])
    annoDF.to_csv(outdir + '/target_genes.bed', sep = '\t', header = False, index = False)

    # Remove intermediate files and return annoDF
    _ = os.system('rm ' + outdir + '/tmp.gtf')
    return annoDF

def FilterGtexMatrix(infile, genelist, tissueDF, metadata, mode, outdir):
    # Construct a shell command for filtering a given GTEx data matrix
    if mode == 'gene_tpm':
        command = 'zcat ' + infile + """ | tail -n +3 | awk '{if(FNR == NR){split($1, gidArr, "."); 
            arr[gidArr[1]] = 1;} else{if($1 == "Name"){print;} else{split($1, gidArr, ".");
            if(gidArr[1] in arr){print;}}}}' """
    
    elif mode == 'junction_counts':
        command = 'zcat ' + infile + """ | tail -n +3 | awk '{if(FNR == NR){split($1, gidArr, "."); 
            arr[gidArr[1]] = 1;} else{if($2 == "Description"){print;} else{split($2, gidArr, ".");
            if(gidArr[1] in arr){print;}}}}' """

    elif mode == 'haplotype_expression':
        command = 'zcat ' + infile + """ | awk '{if(FNR == NR){split($1, gidArr, "."); 
            arr[gidArr[1]] = 1;} else{if($2 == "name"){print;} else{split($2, gidArr, ".");
            if(gidArr[1] in arr){print;}}}}' """
    
    _ = os.system(command + genelist + ' - | gzip > ' + outdir + '/GTEx_v8/tmp.txt.gz')

    # Read in filtered GTEx matrix as a Pandas dataframe
    inDF = pd.read_csv(outdir + '/GTEx_v8/tmp.txt.gz', sep = '\t', compression = 'gzip', header = 0)

    # Filter inDF for samples corresponding to each tissue group
    for tissue in tissueDF.values:
        outfile = outdir + '/GTEx_v8/' + tissue[1] + '/target_genes.' + mode + '.txt'
        inDF[[clm for clm in inDF.columns if clm in {'Name', 'name'} | set(metadata['SAMPID'][metadata['SMTSD'] 
            == tissue[0]])]].to_csv(outfile, sep = '\t', index = False)
    
    # Remove intermediate files
    _ = os.system('rm ' + outdir + '/GTEx_v8/tmp.txt.gz')

def FilterGnomad(vcffile, chrom, start, end, outfile):
    # Pull out chromosome, position, reference allele, alternate allele, maximum allele frequency across gnomAD populations,
    # and maximum number of homozygous individuals across gnomAD populations for variants in vcffile falling within the specified region:
    resDF = pd.DataFrame([[record[0], record[1], record[3], record[4], PullFeature(record[7], 'AF_grpmax_joint=', 'gnomAD'), 
        PullFeature(record[7], 'nhomalt_grpmax_joint=', 'gnomAD')] for record in vcffile.query(chrom, start, end)],
        columns = ['CHROM', 'POS', 'REF', 'ALT', 'AF_GRPMAX_JOINT', 'NHOMALT_GRPMAX_JOINT'])
    
    # Save resDF to outfile
    resDF.to_csv(outfile, sep = '\t', index = False)

def main():
    message = 'Prepares intermediate files that are needed when running STRIPE'
    parser = argparse.ArgumentParser(description = message)

    # Add arguments
    parser.add_argument('--gene_list', metavar = '/path/to/target/gene/list', required = True,
        help = 'path to list of Ensembl IDs for target genes')
    parser.add_argument('--anno_gtf', metavar = '/path/to/gene/annotation/gtf', required = True,
        help = 'path to GTF file with gene annotations from GENCODE')
    parser.add_argument('--gtex_files', metavar = '/path/to/GTEx/files/list', required = False,
        help = 'path to list of file paths for GTEx database files')
    parser.add_argument('--tissue_list', metavar = '/path/to/GTEx/tissues/file', required = False,
        help = 'path to tsv file with two columns: (i) GTEx tissue name, (ii) tissue label')
    parser.add_argument('--gnomad_files', metavar = '/path/to/gnomAD/files/list', required = False,
        help = 'path to list of file paths for gnomAD database files')
    parser.add_argument('--only_bed', action = 'store_true', help = 'only make BED file for input gene list')
    parser.add_argument('--outdir', metavar = '/path/to/output/folder', required = True,
        help = 'path to output folder for STRIPE intermediate files')
    
    # Parse command line arguments
    args = parser.parse_args()

    # Create outdir (if it does not exist yet)
    _ = os.system('mkdir -p ' + args.outdir)

    # Construct a BED file for genes in gene_list using anno_gtf
    geneBed = MakeGeneBed(args.gene_list, args.anno_gtf, args.outdir)

    if not args.only_bed and args.gtex_files is not None and args.tissue_list is not None:
        # Create a folder in outdir to store filtered GTEx v8 database files
        _ = os.system('mkdir -p ' + args.outdir + '/GTEx_v8')

        # Read in tissue_list and make subfolders for each tissue
        tissueDF = pd.read_csv(args.tissue_list, sep = '\t', header = None, names = ['Tissue', 'Label'])
        for label in tissueDF['Label']:
            _ = os.system('mkdir -p ' + args.outdir + '/GTEx_v8/' + label)

        # Read in gtex_files and retrieve basenames for each input file
        gtexDF = pd.read_csv(args.gtex_files, header = None, names = ['Path'])
        gtexDF['Basename'] = gtexDF['Path'].apply(os.path.basename)
    
        # Check if GTEx sample annotations were provided in gtex_files
        if 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt' in set(gtexDF['Basename']):
            # Read in sample annotations as a dataframe
            metadata = pd.read_csv(gtexDF[gtexDF['Basename'] == 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt']['Path'].item(), sep = '\t', header = 0)

            # Filter GTEx gene TPM matrix for genes in genelist
            if 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz' in set(gtexDF['Basename']):
                FilterGtexMatrix(gtexDF[gtexDF['Basename'] == 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz']['Path'].item(),
                    args.gene_list, tissueDF, metadata, 'gene_tpm', args.outdir)
            
            # Filter GTEx splice junction read count matrix for genes in genelist
            if 'GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz' in set(gtexDF['Basename']):
                FilterGtexMatrix(gtexDF[gtexDF['Basename'] == 'GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz']['Path'].item(), 
                    args.gene_list, tissueDF, metadata, 'junction_counts', args.outdir)
            
            # Filter GTEx haplotype expression matrix for genes in genelist
            if 'phASER_GTEx_v8_matrix.txt.gz' in set(gtexDF['Basename']):
                FilterGtexMatrix(gtexDF[gtexDF['Basename'] == 'phASER_GTEx_v8_matrix.txt.gz']['Path'].item(), 
                    args.gene_list, tissueDF, metadata, 'haplotype_expression', args.outdir)
    
    if not args.only_bed and args.gnomad_files is not None:
        # Create a folder in outdir to store filtered gnomAD v4 database files
        _ = os.system('mkdir -p ' + args.outdir + '/gnomAD_v4')

        # Read in gnomadlist and extract chromosomes for each input file
        gnomadDF = pd.read_csv(args.gnomad_files, header = None, names = ['Path'])
        gnomadDF['Chrom'] = gnomadDF['Path'].apply(lambda x: os.path.basename(x).split('.')[5])

        for row in gnomadDF.values:
            # Extract target genes mapping to current chromosome
            for gene in geneBed[geneBed[0] == row[1]].values:
                _ = os.system('mkdir -p ' + args.outdir + '/gnomAD_v4/' + gene[4])

                # Filter gnomAD VCF file for variants mapping to current gene
                FilterGnomad(tabix.open(row[0]), row[1], gene[1], gene[2], args.outdir + '/gnomAD_v4/' + gene[4] + '/gnomad_variants.tsv')

if __name__ == '__main__':
    main()
