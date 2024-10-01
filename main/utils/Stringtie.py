#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2024.09.21
Version: 1.0.0

A script designed to run stringtie on long-read RNA sequencing alignments.
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
    message = 'Runs Stringtie on long-read RNA sequencing alignments'
    parser = argparse.ArgumentParser(description = message)

    # Add arguments
    parser.add_argument('--infile', metavar = '/path/to/input/bam', required = True,
        help = 'path to input file (BAM) with long-read RNA sequencing alignments')
    parser.add_argument('--gencode', metavar = '/path/to/gencode/GTF', required = True,
        help = 'path to GTF file with gene annotations from GENCODE')
    parser.add_argument('--samtools', metavar = '/path/to/samtools', required = True,
        help = 'path to samtools executable')
    parser.add_argument('--stringtie', metavar = '/path/to/stringtie', required = True,
        help = 'path to stringtie executable')
    parser.add_argument('--threads', metavar = '<number of threads>', required = False, 
        help = 'number of threads', type = int, default = 1)
    parser.add_argument('--outprefix', metavar = '/path/to/output/file/prefix', required = True,
        help = 'path to output file prefix')
    
    # Parse command line arguments
    args = parser.parse_args()

    # Create output directory (if it does not exist yet)
    _ = os.system('mkdir -p ' + os.path.dirname(args.outprefix))

    # Filter infile for primary read alignments with MAPQ >= 1
    _ = os.system(args.samtools + ' view -hb -F 256 -q 1 ' + args.infile + ' > ' + args.outprefix + '.filtered.bam')
    _ = os.system(args.samtools + ' index ' + args.outprefix + '.filtered.bam')

    # Run stringtie on filtered read alignments
    _ = os.system(args.stringtie + ' ' + args.outprefix + '.filtered.bam' + ' -G ' + args.gencode + ' -o ' + 
        args.outprefix + '.gtf -p ' + str(args.threads) + ' -L -s 5 -c 5 -u -M 0 -A ' + args.outprefix + '.gene_abundance.tsv')

    # Remove intermediate files
    _ = os.system('rm ' + args.outprefix + '.filtered.bam')
    _ = os.system('rm ' + args.outprefix + '.filtered.bam.bai')

if __name__ == '__main__':
    main()
