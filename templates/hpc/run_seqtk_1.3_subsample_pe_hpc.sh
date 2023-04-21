#!/bin/bash
#SBATCH --job-name=seqtk            # job name
#SBATCH --time=1:00:00              # max job run time
#SBATCH --ntasks-per-node=2         # tasks (commands) per compute node
#SBATCH --cpus-per-task=1           # CPUs (threads) per command
#SBATCH --mem=2G                    # total memory
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

######### SYNOPSIS #########
# This script will subsample paired end short reads

module purge
module load GCC/8.3.0 seqtk/1.3

<<README
    - Seqtk manual: https://github.com/lh3/seqtk
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
pe_1='DATA_DIR/e_coli/seqs/SRR23746297_1.fastq.gz'
pe_2='DATA_DIR/e_coli/seqs/SRR23746297_2.fastq.gz'

######## PARAMETERS ########
subsample_fraction='0.1'        # keep 10%
random_seed=12                  # use the same random seed to maintain pairing

########## OUTPUTS #########
pe_1_out="pe_1_subsample_${subsample_fraction}.fastq.gz"
pe_2_out="pe_2_subsample_${subsample_fraction}.fastq.gz"

################################### COMMANDS ###################################
# the wait command will wait until the two seqtk background process are complete before ending the job script

seqtk sample -s$random_seed $pe_1 $subsample_fraction | gzip > $pe_1_out &
seqtk sample -s$random_seed $pe_2 $subsample_fraction | gzip > $pe_2_out &
wait

################################################################################

<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - Seqtk: https://github.com/lh3/seqtk
CITATION

#2:1:m
