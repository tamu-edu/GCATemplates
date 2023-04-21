#!/bin/bash
#SBATCH --job-name=trinity          # job name
#SBATCH --time=1:00:00              # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=12          # CPUs (threads) per command
#SBATCH --mem=64G                   # total memory per node
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

######### SYNOPSIS #########
# This script will assemble Illumina RNA-seq paired end transcriptomic sequence reads

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 Trinity/2.12.0-Python-3.8.2

<<README
    - Trinity manual: https://github.com/trinityrnaseq/trinityrnaseq/wiki
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
pe_1='DATA_DIR/e_coli/seqs/ERR10828770_1.fastq.gz'
pe_2='DATA_DIR/e_coli/seqs/ERR10828770_2.fastq.gz'

######## PARAMETERS ########
tmpdir=$TMPDIR          # change this if $TMPDIR is not defined on compute nodes
seqType='fq'            # fa, fq
max_memory="$(( $SLURM_MEM_PER_NODE*95/100000 ))G"   # use 95% of SBATCH configured memory; max memory to use by Trinity where limiting can be enabled.
threads=$SLURM_CPUS_PER_TASK
inchworm_cpus=6         # recommended max of 6

########## OUTPUTS #########
# output files are saved to $tmpdir and then copied to pwd after Trinity completes

################################### COMMANDS ###################################
# all files are saved to the compute node disk so they don't count against your file quota
Trinity --seqType $seqType --max_memory $max_memory --left $pe_1 --right $pe_2 --CPU $threads --no_version_check --inchworm_cpu $inchworm_cpus --output $tmpdir/trinity_out

# when Trinity is complete, the following will copy the results files from the $tmpdir to the working directory
cp $tmpdir/trinity_out/Trinity.fasta.gene_trans_map ./
cp $tmpdir/trinity_out/Trinity.fasta ./

################################################################################

<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - Trinity citation:
        Full-length transcriptome assembly from RNA-Seq data without a reference genome.
        Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L,
        Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F,
        Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A.
        Nature Biotechnology 29, 644–652 (2011).
CITATION

#1:C:M
