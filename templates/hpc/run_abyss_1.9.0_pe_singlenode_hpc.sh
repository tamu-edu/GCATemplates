#!/bin/bash
#SBATCH --job-name=abyss            # job name
#SBATCH --time=1:00:00              # max job run time: dd-hh:mm:ss
#SBATCH --nodes=1                   # total compute nodes
#SBATCH --ntasks-per-node=12        # tasks (commands) per compute node
#SBATCH --cpus-per-task=1           # CPUs (threads) per command
#SBATCH --mem=64G                   # total memory per node
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

######### SYNOPSIS #########
# This script will assemble Illumina paired end reads using the ABySS assembler on a single compute node

module purge
module load GCC/8.3.0 OpenMPI/3.1.4 ABySS/2.1.5

<<README
    - ABySS Manual: https://github.com/bcgsc/abyss#abyss
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
pe1_1='DATA_DIR/e_coli/seqs/ERR10828770_1.fastq.gz'
pe1_2='DATA_DIR/e_coli/seqs/ERR10828770_2.fastq.gz'

######## PARAMETERS ########
mpi_threads=$(( SLURM_NTASKS_PER_NODE * SLURM_NNODES ))
serial_threads=$SLURM_CPUS_PER_TASK
kmer=75                             # max kmer available is 128
min_pairs4scaffolding=5             # indicates 5 mate pairs needed to join contigs

########## OUTPUTS #########
prefix='build_abyss_1pe'

################################### COMMANDS ###################################

abyss-pe j=$serial_threads np=$mpi_threads k=$kmer n=$min_pairs4scaffolding name=$prefix lib='lib1' lib1="$pe1_1 $pe1_2"

################################################################################

<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - ABySS:
        Simpson, J. T., Wong, K., Jackman, S. D., Schein, J. E., Jones, S. J., & Birol, I. (2009).
        ABySS: a parallel assembler for short read sequence data. Genome research, 19(6), 1117-1123.
CITATION

#C:1:M
