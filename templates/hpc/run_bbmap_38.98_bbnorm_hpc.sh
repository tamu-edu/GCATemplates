#!/bin/bash
#SBATCH --job-name=bbnorm           # job name
#SBATCH --time=1:00:00              # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=12          # CPUs (threads) per command
#SBATCH --mem=64G                   # total memory per node
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

######### SYNOPSIS #########
# This script will first error correct and then normalize fastq files to a targeted coverage

module purge
module load GCC/11.2.0 BBMap/38.98

<<README
    - BBNorm manual: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbnorm-guide
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
pe_1='DATA_DIR/e_coli/seqs/SRR23746297_1.fastq.gz'
pe_2='DATA_DIR/e_coli/seqs/SRR23746297_2.fastq.gz'

######## PARAMETERS ########
# tadpole.sh params
mode='correct'
correction_kmer=63

# bbnorm.sh params
target_coverage=75
tmpdir=$TMPDIR          # change this if $TMPDIR is not defined on compute nodes

########## OUTPUTS #########
# error corrected files
output_ecc_file1='pe_1_ecc.fastq.gz'
output_ecc_file2='pe_2_ecc.fastq.gz'

# normalized error corrected files
output_norm_file1='pe_1_norm.fastq.gz'
output_norm_file2='pe_2_norm.fastq.gz'
histogram_outfile='norm_hist.txt'

################################### COMMANDS ###################################
# error correct with tadpole first
tadpole.sh in=$pe_1 in2=$pe_2 out=$output_ecc_file1 out2=$output_ecc_file2 mode=$mode k=$correction_kmer

# normalize error corrected reads to $target_coverage
bbnorm.sh in=$output_ecc_file1 in2=$output_ecc_file2 out=$output_norm_file1 out2=$output_norm_file2 hist=$histogram_outfile target=$target_coverage tmpdir=$tmpdir

################################################################################

<<CITATIONS
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - BBMap: https://sourceforge.net/projects/bbmap
CITATIONS

#1:C:M
