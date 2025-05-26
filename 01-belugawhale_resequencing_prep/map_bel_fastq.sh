#!/bin/bash

# Script from Claudio Muller, modified by E.D. to add notes on beluga whale specifics. #Using Bowtie2 program to map samples (Langmead & Salzberg, 2012)

# Call this script with 2 fastq file2 for reads 1 and 2.
#
# e.g. sbatch [thisscript.sh](http://thisscript.sh) myfile.R1.fastq  myfile.R2.fastq
#
# myfile.R1.fastq is referenced by the variable $1
# myfile.R2.fastq is referenced by the variable $2
# $3 is the base name of your bowtie2 index

#SBATCH -J Bowtie2
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 0-12:00 
#SBATCH --mem=24000
#SBATCH -o bt2.%x.out
#SBATCH -e bt2.%x.err

# Load modules
module load StdEnv/2023
module load bowtie2/2.5.2
module load samtools/1.18

# Index reference genome with bowtie2-build before running mapping step. Beluga reference genome ID is GCA_029941455.3 (on NCBI). The reference genome used in this script was previously named "S_20_00703_MOBELS_0202.22xii22.FINAL_scaffolds.fsa" because it was shared/used before it was officially published as GCA_029941455.3.
$BT2_HOME/bowtie2-build $BT2_HOME/reference/S_20_00703_MOBELS_0202.22xii22.FINAL_scaffolds.fsa S_20_00703_MOBELS_0202.22xii22.FINAL_scaffolds

# Run mapping step. This set up is designed to run one sample at a time. 
bowtie2 -x $3 -1 $1 -2 $2 -X 1000 --fr --no-mixed --no-discordant -p 8 | samtools view -b -S - | samtools sort - -o ${1}.bam

# Index bam file.
samtools index ${1}.bam
