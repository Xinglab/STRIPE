#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2024.09.21
Version: 1.0.0

A script designed to run DeepVariant on input sequencing alignment files
Note that this script assumes the following: 
    * Singularity is installed and available in your $PATH
    * Docker image files for DeepVariant (version 1.6.1) or DeepTrio (version 1.6.1) have been pulled
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
    message = 'Runs DeepVariant on input sequencing alignment files'
    parser = argparse.ArgumentParser(description = message)

    # Add arguments
    parser.add_argument('--run_mode', required = False, default = 'standard',
        help = 'run mode for DeepVariant (standard or trio)')
    parser.add_argument('--model_type', required = True,
        help = 'model type for DeepVariant (WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA)')
    parser.add_argument('--genome', metavar = '/path/to/reference/genome/FASTA', required = True,
        help = 'path to FASTA file containing reference genome sequence')
    parser.add_argument('--bam_proband', metavar = '/path/to/input/BAM/proband', required = True,
        help = 'path to input read alignments (BAM) for proband')
    parser.add_argument('--bam_parent1', metavar = '/path/to/input/BAM/parent1', required = False,
        help = 'path to input read alignments (BAM) for parent 1')
    parser.add_argument('--bam_parent2', metavar = '/path/to/input/BAM/parent2', required = False,
        help = 'path to input read alignments (BAM) for parent 2')
    parser.add_argument('--label_proband', required = False, default = 'PROBAND',
        help = 'sample label for proband')
    parser.add_argument('--label_parent1', required = False, default = 'PROBAND_P1',
        help = 'sample label for parent 1')
    parser.add_argument('--label_parent2', required = False, default = 'PROBAND_P2',
        help = 'sample label for parent 2')
    parser.add_argument('--threads', metavar = '<number of threads>', required = False, 
        help = 'number of threads', type = int, default = 1)
    parser.add_argument('--outdir', metavar = '/path/to/outdir', required = True,
        help = 'path to output directory')

    # Parse command line arguments
    args = parser.parse_args()

    # Create outdir (if it does not exist) and a temporary directory in outdir
    _ = os.system('mkdir -p ' + args.outdir)
    _ = os.system('mkdir -p ' + args.outdir + '/tmp')

    # Construct bind path using folders for genome and outdir
    bind_path = ','.join(['/usr/lib/locale/:/usr/lib/locale/', os.path.dirname(args.genome), args.outdir])

    if args.run_mode == 'standard':
        # Add folder for bam_proband to bind_path
        bind_path = ','.join([bind_path, os.path.dirname(args.bam_proband)])

        _ = os.system('singularity run -B ' + bind_path + ' docker://google/deepvariant:1.6.1 ' +
            '/opt/deepvariant/bin/run_deepvariant --model_type=' + args.model_type + 
            ' --ref=' + args.genome + ' --reads=' + args.bam_proband + ' --sample_name ' +
            args.label_proband + ' --output_vcf=' + args.outdir + '/' + args.label_proband + '.vcf.gz' +
            ' --output_gvcf=' + args.outdir + '/' + args.label_proband + '.g.vcf.gz ' +
            '--intermediate_results_dir ' + args.outdir + '/tmp --num_shards=' + str(args.threads))

    elif args.run_mode == 'trio':
        if args.bam_parent1 is not None and args.bam_parent2 is not None:
            # Add folders for bam_proband, bam_parent1, and bam_parent2 to bind_path
            bind_path = ','.join([bind_path, os.path.dirname(args.bam_proband), os.path.dirname(args.bam_parent1),
                os.path.dirname(args.bam_parent2)])
            
            _ = os.system('singularity run -B ' + bind_path + ' docker://google/deepvariant:deeptrio-1.6.1 ' +
                '/opt/deepvariant/bin/deeptrio/run_deeptrio --model_type=' + args.model_type + 
                ' --ref=' + args.genome + ' --reads_child=' + args.bam_proband + ' --reads_parent1=' +
                args.bam_parent1 + ' --reads_parent2=' + args.bam_parent2 + ' --sample_name_child ' +
                args.label_proband + ' --sample_name_parent1 ' + args.label_parent1 + ' --sample_name_parent2 ' +
                args.label_parent2 + ' --output_vcf_child ' + args.outdir + '/' + args.label_proband + '.vcf.gz' +
                ' --output_vcf_parent1 ' + args.outdir + '/' + args.label_parent1 + '.vcf.gz' +
                ' --output_vcf_parent2 ' + args.outdir + '/' + args.label_parent2 + '.vcf.gz' +
                ' --output_gvcf_child ' + args.outdir + '/' + args.label_proband + '.g.vcf.gz' +
                ' --output_gvcf_parent1 ' + args.outdir + '/' + args.label_parent1 + '.g.vcf.gz' +
                ' --output_gvcf_parent2 ' + args.outdir + '/' + args.label_parent2 + '.g.vcf.gz' +
                ' --intermediate_results_dir ' + args.outdir + '/tmp --num_shards ' + str(args.threads))

if __name__ == '__main__':
    main()
