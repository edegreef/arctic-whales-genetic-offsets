#!/bin/bash

#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=trim_fastq
#SBATCH --output=%x-%j.out

# Trim raw sequence files with Trimmomatic program (Bolger et al., 2014)

# Load modules
module load nixpkgs/16.09 trimmomatic/0.36

# Use TruSeq3-PE-2.fa to account for read-through, obtained at https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE-2.fa
# Path for TruSeq3-PE-2.fa file needs to be updated on line 24.

# Make trimmed folder for all the outputs
mkdir trimmed
	
# Run trimmomatic	
ls ./*_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz$//' | parallel --jobs $SLURM_CPUS_PER_TASK 'java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE {}_R1_001.fastq.gz {}_R2_001.fastq.gz trimmed/{}_R1_paired.fastq.gz trimmed/{}_R1_unpaired.fastq.gz trimmed/{}_R2_paired.fastq.gz trimmed/{}_R2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36'