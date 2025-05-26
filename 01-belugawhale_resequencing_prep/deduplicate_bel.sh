#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=300GB
#SBATCH --job-name=deduplicateRG
#SBATCH --output=%x-%j.out

# Removing duplicate reads and adding read group information with Picard program (Broad Institute 2019)

# Load modules
module load nixpkgs/16.09 picard/2.20.6

# Noting here need to be in directory that has sorted.bam files, or use cd to navigate there in this script. 

# Also noting beluga bam files at this stage were named in this format: CS_Pagnirtung_1022_1_R1.fastq.gz.bam. 

# Remove duplicate reads for all samples in current folder. # This turns it into example: CS_Pagnirtung_1022_1_R1.fastq.gz.deDup.bam
ls *.bam | sed 's/.bam$//' | parallel --jobs 4 'java -Xmx50000m -jar $EBROOTPICARD/picard.jar MarkDuplicates I={}.bam O={}.deDup.bam M={}_deDupMetrics.txt REMOVE_DUPLICATES=true USE_JDK_DEFLATER=true USE_JDK_INFLATER=true'

# Add read group information to each deduplicated bam file. # This turns it into example: CS_Pagnirtung_1022_1_R1.fastq.gz.deDupRG.bam
ls .*deDup.bam | sed 's/.deDup.bam$//' | parallel --jobs 4 'java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I={}.deDup.bam O={}.deDupRG.bam RGID={} RGPL=illumina RGSM={} RGLB=lib1 RGPU=unit1'

# Recommend keeping "deDupRG.bam" as suffix (e.g., CS_Pagnirtung_1022_1_R1.fastq.gz.deDupRG.bam). 
# Though if want to keep same name as what was used previously, can rename all bam files with lines 27-34 to match VCF header (e.g., "no_dups_CS_Pagnirtung_1022_1_R1.fastq.gz.bam")

## Remove last 11 characters (deDupRG.bam) for all file names in folder
#for i in *; do mv "$i" "${i%????????????}"; done

## add "bam" to each file name in folder
#for i in *; do mv "$i" "$i.bam"; done

## add "no_dups_" to the beginning of each file name
#for i in *; do mv "$i" "no_dups_$i"; done

# Load samtools to index final bam files
module load gcc/8.3.0 samtools/1.9

# Index new bam files (adjust "*.deDupRG.bam" suffix if using different suffix)
ls *.deDupRG.bam | sed 's/.deDupRG.bam$//' | parallel --jobs 4 'samtools index {}.deDupRG.bam'
