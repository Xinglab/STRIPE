#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2024.09.23
Version: 1.0.0

A script designed to run LongcallR on alignments of ONT cDNA sequencing reads. 
'''

# =====================================================================================================================
#                                                       PACKAGES 
# =====================================================================================================================

# Load required packages
import sys, argparse, os

# =====================================================================================================================
#                                                    MAIN FUNCTIONS
# =====================================================================================================================

def main():
    message = 'Runs LongcallR on alignments of ONT cDNA sequencing reads'
    parser = argparse.ArgumentParser(description = message)

    # Add arguments
    parser.add_argument('--infile', metavar = '/path/to/input/BAM', required = True,
        help = 'path to input file (BAM) containing alignments of ONT cDNA sequencing reads')
    parser.add_argument('--longcallr', metavar = '/path/to/longcallR', required = True,
        help = 'path to longcallR executable')
    parser.add_argument('--bcftools', metavar = '/path/to/bcftools', required = True,
        help = 'path to bcftools executable')
    parser.add_argument('--genome', metavar = '/path/to/reference/genome/FASTA', required = True,
        help = 'path to FASTA file containing reference genome sequence')
    parser.add_argument('--outprefix', metavar = '/path/to/output/file/prefix', required = True,
        help = 'path to output file prefix')
    parser.add_argument('--threads', metavar = '<number of threads>', required = False, 
        help = 'number of threads', type = int, default = 1)

    # Parse command line arguments
    args = parser.parse_args()

    # Create outdir (if it does not exist)
    _ = os.system('mkdir -p ' + os.path.dirname(args.outprefix))

    # Run longcallR on infile
    _ = os.system(args.longcallr + ' --bam-path ' + args.infile + ' --ref-path ' + args.genome +
        ' --output ' + args.outprefix + ' --threads ' + str(args.threads) + ' --platform ont ' +
        '--preset ont-cdna --no-bam-output')
    
    # Sort, compress, and index the VCF file from longcallR using bcftools
    _ = os.system(args.bcftools + ' sort ' + args.outprefix + '.vcf -Oz -o ' + args.outprefix +
        '.vcf.gz --write-index="tbi"')
    _ = os.system('rm ' + args.outprefix + '.vcf')

if __name__ == '__main__':
    main()
