#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2024.09.25
Version: 1.0.0

A script designed to run Clair3 on input sequencing alignment files
Note that this script assumes that Singularity is installed and available in your $PATH
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
    message = 'Runs Clair3 on input sequencing alignment files'
    parser = argparse.ArgumentParser(description = message)

    # Add arguments
    parser.add_argument('--clair3', metavar = '/path/to/clair3/sif', required = True,
        help = 'path to Clair3 singularity image file')
    parser.add_argument('--infile', metavar = '/path/to/input/BAM', required = True,
        help = 'path to input read alignments (BAM)')
    parser.add_argument('--genome', metavar = '/path/to/reference/genome/FASTA', required = True,
        help = 'path to FASTA file containing reference genome sequence')
    parser.add_argument('--threads', metavar = '<number of threads>', required = False, 
        help = 'number of threads', type = int, default = 1)
    parser.add_argument('--model', metavar = '/path/to/clair3/model', required = True,
        help = 'path to Clair3 model file')
    parser.add_argument('--outdir', metavar = '/path/to/outdir', required = True,
        help = 'path to output directory')

    # Parse command line arguments
    args = parser.parse_args()

    # Create outdir (if it does not exist)
    _ = os.system('mkdir -p ' + args.outdir)

    # Construct bind path using outdir and folders for infile, genome, and model
    bind_path = ','.join([args.outdir, os.path.dirname(args.infile), os.path.dirname(args.genome), 
        os.path.dirname(args.model)])
    
    _ = os.system('singularity exec -B ' + bind_path + ' ' + args.clair3 + ' /opt/bin/run_clair3.sh ' +
        '--bam_fn=' + args.infile + ' --ref_fn=' + args.genome + ' --threads=' + str(args.threads) +
        ' --platform=ont --pileup_only --model_path=' + args.model + ' --output=' + args.outdir)

if __name__ == '__main__':
    main()
