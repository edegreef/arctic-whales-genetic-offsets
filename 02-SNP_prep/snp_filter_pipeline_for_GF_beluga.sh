#!/bin/bash

####### Script for filtering beluga whale SNPs to prepare for gradient forest analyses
## Need bcftools (Danecek et al., 2021) and vcftools (Danecek et al., 2011) installed
#sudo apt-get install vcftools
#sudo apt-get install bcftools


# Set up variable before running
input_prefix=146Beluga
snps_prefix=beluga_snps
min_scaf_list=scaf_min100kb_beluga #list of scaffolds to keep


## 1) Remove indels - run this line if starting from raw variants. 
#vcftools --gzvcf 146Beluga_raw.ID.vcf.gz --remove-indels --recode --recode-INFO-all --stdout > $input_prefix.vcf

# Zip and index
bgzip $input_prefix.vcf
tabix -p vcf $input_prefix.vcf.gz

## 2) Keep sites with "PASS"
vcftools --gzvcf $input_prefix.vcf.gz --remove-filtered-all --recode --recode-INFO-all --stdout > $snps_prefix.PASS.vcf

## 3) Quality filter (minQ 50)
vcftools --vcf $snps_prefix.PASS.vcf --minQ 50 --recode --recode-INFO-all --stdout > $snps_prefix.PASS.minQ50.vcf


## 4) Max missingness 25%
vcftools --vcf $snps_prefix.PASS.minQ50.vcf --max-missing 0.75 --recode --recode-INFO-all --stdout > $snps_prefix.PASS.minQ50.miss.vcf


## 5) Remove non-biallic sites
vcftools --vcf $snps_prefix.PASS.minQ50.miss.vcf --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all --stdout > $snps_prefix.PASS.minQ50.miss.biallel.vcf


## 6) Remove scaffolds < 100kb - few extra steps here
# Zip & index
bgzip $snps_prefix.PASS.minQ50.miss.biallel.vcf
tabix -p vcf $snps_prefix.PASS.minQ50.miss.biallel.vcf.gz

# Convert scaffold list to one-liner to use in bcftools
awk '{print $1}' $min_scaf_list | paste -s -d, - > scaf_list_line

# Set up list
list=`cat scaf_list_line`

# Filter vcf for these scaffolds
bcftools filter --regions $list $snps_prefix.PASS.minQ50.miss.biallel.vcf.gz > $snps_prefix.PASS.minQ50.miss.biallel.min100kb.vcf

# Can make a list of scaffolds to double check
grep -v "^#" $snps_prefix.PASS.minQ50.miss.biallel.min100kb.vcf | cut -f1 | sort | uniq > filtered_contig_list_check.txt

# Zip & index
bgzip $snps_prefix.PASS.minQ50.miss.biallel.min100kb.vcf
tabix -p vcf $snps_prefix.PASS.minQ50.miss.biallel.min100kb.vcf.gz


## 7) Filter out sex-linked SNPS to create autosomal dataset - removing X chr
vcftools --gzvcf $snps_prefix.PASS.minQ50.miss.biallel.min100kb.vcf.gz --not-chr S_20_00703_ChromosomeX_scaffold_8 --not-chr S_20_00703_ChromosomeX_scaffold_9 --recode --recode-INFO-all --stdout > $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.vcf


## 8) Remove close-kins and individuals with high missingness
vcftools --vcf $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.vcf --remove-indv EHB_Kuujjuarapik_4594 --remove-indv JB_Bear_Islands_1029 --remove-indv JB_Pointe_de_Repentigny_2340 --remove-indv SL_Le_Bic_4287 --remove-indv SL_NA_55 --remove-indv SL_Sainte_Flavie_5792 --recode --recode-INFO-all --stdout > $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n140.vcf

## 9) Minor allele frequency filter
vcftools --vcf $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n140.vcf --maf 0.05 --recode --recode-INFO-all --stdout > $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n140.maf.vcf


## 10) Additional missingness filter (stricter filter prior to imputing) for 10% max missingness
vcftools --vcf $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n140.maf.vcf --max-missing 0.9 --recode --recode-INFO-all --stdout > $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n140.maf.miss01.vcf

# Lastly, can rename file
mv $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n140.maf.miss01.vcf beluga_snps.filtered.n140.vcf
