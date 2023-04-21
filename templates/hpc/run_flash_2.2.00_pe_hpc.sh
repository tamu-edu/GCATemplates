#!/bin/bash
#SBATCH --job-name=flash2           # job name
#SBATCH --time=1:00:00              # max job run time
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=3           # CPUs (threads) per command
#SBATCH --mem=2G                    # total memory
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

######### SYNOPSIS #########
# This script will merge overlapping paired end short reads into a single sequence

module purge
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 FLASH/2.2.00

<<README
    - FLASH homepage: https://github.com/dstreett/FLASH2
    - FLASH (Fast Length Adjustment of SHort reads) is an accurate and fast tool to merge
        paired-end reads that were generated from DNA fragments whose
        lengths are shorter than twice the length of reads.
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
pe1_1='DATA_DIR/e_coli/seqs/SRR23746297_1.fastq.gz'
pe1_2='DATA_DIR/e_coli/seqs/SRR23746297_2.fastq.gz'

######## PARAMETERS ########
threads=$SLURM_CPUS_PER_TASK
min_overlap=10                  # default 10

avg_read_length=250             # default 100
avg_frag_length=400             # default 180
stdev_frag_length=50            # default 20

########## OUTPUTS #########
output_prefix='SRR1056110'

################################### COMMANDS ###################################

flash2 -z -t $threads -m $min_overlap -r $avg_read_length -f $avg_frag_length -s $stdev_frag_length -o $output_prefix $pe1_1 $pe1_2

################################################################################

<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - FLASH:
        FLASH: Fast length adjustment of short reads to improve genome assemblies. T. Magoc and S. Salzberg.
        Bioinformatics 27:21 (2011), 2957-63.
CITATION

#1:3:m
