#!/bin/bash
#SBATCH --job-name=freebayes        # job name
#SBATCH --time=1:00:00              # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=12          # CPUs (threads) per command
#SBATCH --mem=64G                   # total memory per node
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

######### SYNOPSIS #########
# This script will call variants using a reference genome and aligned .bam file of Illumina reads

module purge
module load GCC/11.2.0  OpenMPI/4.1.1 freebayes/1.3.6-R-4.1.2

<<README
    - FreeBayes manual: https://github.com/ekg/freebayes
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
reference_fasta='DATA_DIR/e_coli/ref/GCF_000005845.2_ASM584v2_genomic.fna'
alignments_bam='DATA_DIR/e_coli/ref/aln/SRR23746297_sorted.bam'

######## PARAMETERS ########
sample_name='SRR23746297'

########## OUTPUTS #########
output_vcf="${sample_name}_freebayes_out.vcf"

################################### COMMANDS ###################################

freebayes --fasta-reference $reference_fasta --bam $alignments_bam > $output_vcf

################################################################################

<<CITATIONS
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - FreeBayes:
        Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing.
        arXiv preprint arXiv:1207.3907 [q-bio.GN] 2012
CITATIONS

#1:C:M
