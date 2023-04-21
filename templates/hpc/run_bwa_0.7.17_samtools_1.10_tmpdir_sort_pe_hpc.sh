#!/bin/bash
#SBATCH --job-name=bwa_samtools     # job name
#SBATCH --time=1:00:00              # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=12          # CPUs (threads) per command
#SBATCH --mem=64G                   # total memory per node
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

module purge
module load GCC/9.3.0 BWA/0.7.17 SAMtools/1.10

######### SYNOPSIS #########
# aligns to a reference; uses $TMPDIR for sorting; output is a sorted bam file

<<README
    - BWA manual: http://bio-bwa.sourceforge.net/bwa.shtml
    - SAMtools manual: http://samtools.github.io/hts-specs/SAMv1.pdf
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
pe1_1='DATA_DIR/e_coli/seqs/SRR23746297_1.fastq.gz'
pe1_2='DATA_DIR/e_coli/seqs/SRR23746297_2.fastq.gz'

ref_genome='DATA_DIR/e_coli/ref/GCF_000005845.2_ASM584v2_genomic.fna'

######## PARAMETERS ########
tmpdir=$TMPDIR          # change this if $TMPDIR is not defined on compute nodes
cpus=$SLURM_CPUS_PER_TASK
readgroup='ecoli_sra'
library='pe'
sample='SRR23746297'
platform='ILLUMINA'     # ILLUMINA, CAPILLARY, LS454, SOLID, HELICOS, IONTORRENT, ONT, PACBIO

########## OUTPUTS #########
output_bam="${sample}_sorted.bam"

################################### COMMANDS ###################################
# NOTE index genome only if not using already indexed genome
if [ ! -f ${ref_genome}.bwt ]; then
  bwa index $ref_genome
fi

bwa mem -M -t $cpus -R "@RG\tID:$readgroup\tLB:$library\tSM:$sample\tPL:$platform" $ref_genome $pe1_1 $pe1_2 -o $tmpdir/${sample}_aln.sam

samtools sort $tmpdir/${sample}_aln.sam -o $output_bam -m 1G -@ $cpus -T $tmpdir/tmp4sort

################################################################################

<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - BWA:
        Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25, 1754-1760.

    - SAMtools:
        Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project
        Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9.
CITATION

#1:C:M
