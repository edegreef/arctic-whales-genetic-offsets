#!/bin/bash

#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=300M
#SBATCH --job-name=edit_IDs
#SBATCH --output=%x-%j.out

# Quick script for editing labels in the vcf, shortening sample ID to remove unncessary characters/paths in sample IDs
# Using vcftools (Danecek et al., 2011) and bcftools (Danecek et al., 2021)

# Load modules
module load StdEnv/2023 vcftools/0.1.16 gcc/12.3 bcftools/1.19

# Zip and index
bgzip 146Beluga_raw.vcf
tabix -p vcf 146Beluga_raw.vcf.gz

# Use "vcf_sample_id.txt" to modify header in vcf file to change sample ID to the shortened ID list
bcftools reheader -s vcf_sample_id.txt 146Beluga_raw.vcf.gz -o 146Beluga_raw.ID.vcf.gz
tabix -p vcf 146Beluga_raw.ID.vcf.gz
