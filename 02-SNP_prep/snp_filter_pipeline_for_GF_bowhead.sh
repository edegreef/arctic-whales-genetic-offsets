#!/bin/bash

####### Script for filtering bowhead whale SNPs to prepare for gradient forest analyses
## Need bcftools (Danecek et al., 2021) and vcftools (Danecek et al., 2011) installed
#sudo apt-get install vcftools
#sudo apt-get install bcftools


# Set up variables before running
variant_prefix=bowhead_RWmap_allvariants.n20
snps_prefix=bowhead_RWmap_snps
min_scaf_list=scaf_min100kb_bowhead_RWmap #list of scaffolds to keep


## 0) Noting here that the bowhead_RWmap_allvariants.n20.vcf.gz file excludes one Newfoundland whale (WBF_2005_0298) because this individual was found stranded and is not be ideal for using locality for environmental data in this study. This is mentioned here because it is included in the raw sequencing data prep in de Greef et al. 2024 (https://doi.org/10.1111/gcb.17528), but was removed after SNP calling and before filtering for this current study.

## 1) Remove indels
vcftools --gzvcf $variant_prefix.vcf.gz --remove-indels --recode --recode-INFO-all --stdout > $snps_prefix.vcf


## 2) Keep sites with "PASS"
vcftools --vcf $snps_prefix.vcf --remove-filtered-all --recode --recode-INFO-all --stdout > $snps_prefix.PASS.vcf


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


## 7) Filter out sex-linked SNPS to create autosomal dataset - straightforward for data mapped to chr-level reference genome because we have X and Y chr
vcftools --gzvcf $snps_prefix.PASS.minQ50.miss.biallel.min100kb.vcf.gz --not-chr CM053059.1 --not-chr CM053080.1 --recode --recode-INFO-all --stdout > $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.vcf


## 8) Remove close-kins and individuals with high missingness
vcftools --vcf $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.vcf --remove-indv 99_01 --recode --recode-INFO-all --stdout > $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n19.vcf


## 9) Minor allele frequency filter
vcftools --vcf $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n19.vcf --maf 0.05 --recode --recode-INFO-all --stdout > $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n19.maf.vcf


## 10) Additional missingness filter (stricter filter prior to imputing) for 10% max missingness
vcftools --vcf $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n19.maf.vcf --max-missing 0.9 --recode --recode-INFO-all --stdout > $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n19.maf.miss01.vcf

# Lastly, can rename file
mv $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n19.maf.miss01.vcf $snps_prefix.filtered.n19.vcf

