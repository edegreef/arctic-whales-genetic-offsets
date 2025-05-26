#!/bin/bash

####### Script for imputing and thinning SNPs to prepare for gradient forest analyses
## Need vcftools (Danecek et al., 2011), bcftools (Danecek et al., 2021), and beagle (Browning et al., 2018, 2021) installed.
#sudo apt-get install vcftools
#sudo apt-get install bcftools

# Beagle program can be downloaded from https://faculty.washington.edu/browning/beagle/b5_2.html
# The specific jar file used in this script is beagle.28Jun21.220.jar (https://faculty.washington.edu/browning/beagle/beagle.28Jun21.220.jar) and requires Java version 8.

# Noting here that beagle runs into error when there are scaffolds with just 1 snp. This was applicable in the narwhal and beluga datasets of the three whales -- see "remove_single_snp_scaf_notes.txt" for preparing narwhal and beluga dataset for this prior to imputation. 
# If don't want to run the extra steps from "remove_single_snp_scaf_notes.txt", can run line 13 for narwhal, and line 14 for beluga whale, as these scaffolds have been identified.
#vcftools --vcf narwhal_snps.filtered.n57.vcf --not-chr SIHG01000985.1 --not-chr SIHG01003364.1 --not-chr SIHG01003795.1 --not-chr SIHG01005513.1 --recode --recode-INFO-all --stdout > narwhal_snps.filtered.n57.sm.vcf
#vcftools --vcf beluga_snps.filtered.n140.vcf --not-chr HiC_scaffold_1354 --recode --recode-INFO-all --stdout > beluga_snps.filtered.n140.sm.vcf

# bowhead/narwhal/beluga starting file names:
# Narwhal: narwhal_snps.filtered.n57.sm.vcf 
# Bowhead: bowhead_RWmap_snps.filtered.n19.vcf
# Beluga: beluga_snps.filtered.n140.sm.vcf

# Set up snp variable
snp_input=beluga_snps.filtered.n140.sm
prefix_output=beluga_snps.filtered.n140

# Zip vcf
bgzip $snp_input.vcf
tabix -p vcf $snp_input.vcf.gz

# Run beagle for imputation! (this will phase snps as well)
java -Xmx70g -jar beagle.28Jun21.220.jar gt=$snp_input.vcf.gz iterations=20 out=$prefix_output.imputed

# Thin snps
vcftools --gzvcf $prefix_output.imputed.vcf.gz --recode --recode-INFO-all --thin 1000 --stdout > $prefix_output.imputed.thin.vcf

# Zip & index
bgzip $prefix_output.imputed.thin.vcf
tabix -p vcf $prefix_output.imputed.thin.vcf.gz

# Convert to ped/map file for LEA
bcftools view -H $prefix_output.imputed.thin.vcf.gz | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom-map.txt
vcftools --gzvcf $prefix_output.imputed.thin.vcf.gz --plink --chrom-map chrom-map.txt --out $prefix_output.imputed.thin
