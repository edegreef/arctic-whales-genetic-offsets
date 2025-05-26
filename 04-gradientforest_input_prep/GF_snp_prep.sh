#!/bin/bash

####### Script for selecting target SNPs (identified by "LEA_lfmms.R") and then formatting file for inputting into GradientForest
## Need vcftools (Danecek et al., 20211) installed, and the "vcf2forR.sh" needs to be in working directory.
#sudo apt-get install vcftools

# narwhal file: narwhal_snps.filtered.n57.imputed.thin.vcf.gz
# beluga file: beluga_snps.filtered.n140.imputed.thin.vcf.gz
# N/A for bowhead whale

# set up snp prefix variable
snp_prefix=narwhal_snps.filtered.n57.imputed.thin
list=narwhal_K3_top0.01_scan_4env_merged_snplist.txt
output_prefix=narwhal_snps_filtered_K3_top0.01_scan_4env

# Filtering for just the snp list we determined from genome scan and lfmms). For narwhal and beluga, this should contain 47.5K SNPs
vcftools --gzvcf $prefix.vcf.gz --positions $list --recode --recode-INFO-all --stdout > $output_prefix.vcf

# Run vcf2forR.sh to convert to format for GradientForest
./vcf2forR.sh $output_prefix
