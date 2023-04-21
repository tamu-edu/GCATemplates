#!/bin/bash
#SBATCH --job-name=hisat2           # job name
#SBATCH --time=1:00:00              # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=12          # CPUs (threads) per command
#SBATCH --mem=64G                   # total memory per node
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file

######### SYNOPSIS #########
# This template script aligns paired end reads and sorts the output into a bam file

<<README
    - HISAT2 manual: http://ccb.jhu.edu/software/hisat2/manual.shtml
README

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 HISAT2/2.2.1
module load Python/3.8.2 SAMtools/1.10

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
pe_1='DATA_DIR/e_coli/seqs/ERR10828770_1.fastq.gz'
pe_2='DATA_DIR/e_coli/seqs/ERR10828770_2.fastq.gz'
ref_genome='DATA_DIR/e_coli/ref/GCF_000005845.2_ASM584v2_genomic.fna'

######## PARAMETERS ########
tmpdir=$TMPDIR          # change this if $TMPDIR is not defined on compute nodes
threads=$SLURM_CPUS_PER_TASK
# read group information
id='ecoli_sra'
library='sra_pe'
platform='ILLUMINA'
sample='ERR10828770'

########## OUTPUTS #########
output_bam="${sample}_pe_aln.bam"

################################### COMMANDS ###################################
# NOTE index genome only if not using already indexed genome
if [ ! -f ${ref_genome}.1.ht2 ]; then
 hisat2-build $ref_genome $ref_genome
fi

# If you will be using cufflinks downstream, run hisat2 with the --dta-cufflinks option instead of --dta
hisat2 --dta -p $threads --rg-id "$id" --rg "LB:$library" --rg "SM:$sample" --rg "PL:$platform" -x $ref_genome -q -1 $pe_1 -2 $pe_2 -S $tmpdir/out.sam

samtools view -bS $tmpdir/out.sam | samtools sort -m 1G -@ $threads - -T $tmpdir/$sample -o $output_bam

################################################################################

<<CITATION
    - Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

    - HISAT2:
        Kim D, Langmead B and Salzberg SL. HISAT: a fast spliced aligner with low memory requirements.
        Nature Methods 12, 357–360 (2015). doi:10.1038/nmeth.3317
    -SAMtools:
	Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis R, Durbin R, and 1000 Genome Project
	Data Processing Subgroup. The Sequence alignment/map (SAM) format and SAMtools. 
	Bioinformatics 25, 2078-2079 (2021). doi:10.1093/bioinformatics/btp352 

CITATION

#1:C:M
