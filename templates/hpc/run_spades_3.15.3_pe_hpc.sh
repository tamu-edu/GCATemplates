#!/bin/bash
#SBATCH --job-name=spades           # job name
#SBATCH --time=1:00:00              # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=12          # CPUs (threads) per command
#SBATCH --mem=64G                   # total memory per node
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

######### SYNOPSIS #########
# This script will assemble Illumina paired end reads using the spades assembler

module purge
module load GCC/10.2.0 SPAdes/3.15.3

<<README
    - SPAdes manual: http://spades.bioinf.spbau.ru/release3.5.0/manual.html
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
pe_1='DATA_DIR/e_coli/seqs/SRR23746297_1.fastq.gz'
pe_2='DATA_DIR/e_coli/seqs/SRR23746297_2.fastq.gz'

######## PARAMETERS ########
threads=$SLURM_CPUS_PER_TASK
max_memory=$SLURM_MEM_PER_NODE

########## OUTPUTS #########
output_dir='SRR23746297_out'

################################### COMMANDS ###################################

spades.py -1 $pe_1 -2 $pe_2 --threads $threads --memory $max_memory -o $output_dir

################################################################################

<<CITATIONS
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - SPAdes:
        Bankevich A., et al. SPAdes: A New Genome Assembly Algorithm 
        and Its Applications to Single-Cell Sequencing.
        J Comput Biol. 2012 May; 19(5): 455–477. doi:  10.1089/cmb.2012.0021
CITATIONS

#1:C:M
