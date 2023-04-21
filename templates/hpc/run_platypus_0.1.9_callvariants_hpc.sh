#!/bin/bash
#SBATCH --job-name=platypus         # job name
#SBATCH --time=1:00:00              # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=12          # CPUs (threads) per command
#SBATCH --mem=64G                   # total memory per node
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

######### SYNOPSIS #########
# This script will call variants using the Platypus tool

module purge
module load GCC/10.2.0  OpenMPI/4.0.5  Platypus/0.8.1-Python-2.7.18

<<README
    - Platypus manual: https://www.rdm.ox.ac.uk/research/lunter-group/lunter-group/platypus-documentation
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
reference_genome='DATA_DIR/e_coli/ref/GCF_000005845.2_ASM584v2_genomic.fna'
alignments_bam='DATA_DIR/e_coli/ref/aln/SRR23746297_sorted.bam'

######## PARAMETERS ########
sample_name='SRR23746297'

########## OUTPUTS #########
output_vcf="${sample_name}_platypus_out.vcf"

################################### COMMANDS ###################################

Platypus.py callVariants --nCPU=$SLURM_CPUS_PER_TASK --bamFiles=$alignments_bam --refFile=$reference_genome --output=$output_vcf

################################################################################

<<CITATIONS
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - Platypus:
            Andy Rimmer, Hang Phan, Iain Mathieson, Zamin Iqbal, Stephen R. F. Twigg, WGS500 Consortium, Andrew O. M. Wilkie,
            Gil McVean, Gerton Lunter. Integrating mapping-, assembly- and haplotype-based approaches for calling variants
            in clinical sequencing applications. Nature Genetics (2014) doi:10.1038/ng.3036
CITATIONS

#1:C:M
