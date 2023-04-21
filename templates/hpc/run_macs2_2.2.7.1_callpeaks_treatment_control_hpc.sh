#!/bin/bash
#SBATCH --job-name=macs2            # job name
#SBATCH --time=1:00:00              # max job run time
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=1           # CPUs (threads) per command
#SBATCH --mem=2G                    # total memory
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

######### SYNOPSIS #########
# This script will call peaks in chip-seq data and create a peak model image .pdf file

module purge
module load GCC/11.2.0 OpenMPI/4.1.1 MACS2/2.2.7.1
module load R/4.1.2

<<README
    - MACS: Model-based Analysis for ChIP-Sequencing
    - MACS homepage: http://liulab.dfci.harvard.edu/MACS
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
control_file='DATA_DIR/chipseq/UW_K562_H3K4me3_Control.bam'
treatment_file='DATA_DIR/chipseq/UW_K562_H3K4me3.bam'

######## PARAMETERS ########
genome_size='hs'        # shortcuts:'hs' for human, 'mm' for mouse, 'ce' for C. elegans and 'dm' for fruitfly, Default:hs
tempdir="$TMPDIR"       # use a different directory if $TMPDIR does not exist on your cluster

########## OUTPUTS #########
output_dir='out_macs2'
output_prefix='UW_K562_H3K4me3'

################################### COMMANDS ###################################

macs2 callpeak -g $genome_size --tempdir $tempdir -c $control_file -t $treatment_file -n $output_prefix --outdir $output_dir

# create pdf for peak model
Rscript $output_dir/${output_prefix}_model.r

################################################################################

<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - MACS:
        Zhang et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biol (2008) vol. 9 (9) pp. R137.
CITATION

#1:1:m
