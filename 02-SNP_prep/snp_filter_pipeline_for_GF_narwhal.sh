#!/bin/bash

####### Script for filtering narwhal SNPs to prepare for gradient forest analyses
## Need bcftools (Danecek et al., 2021) and vcftools (Danecek et al., 2011) installed
#sudo apt-get install vcftools
#sudo apt-get install bcftools


# Set up variables before running
variant_prefix=narwhal_allvariants.n60
snps_prefix=narwhal_snps
min_scaf_list=scaf_min100kb_narwhal #list of scaffolds to keep
xchr_scaf_list=final_scaffolds_Xlinked_narwhal.txt #list of scaffolds in X chromosome
ychr_scaf_list=final_scaffolds_Ylinked_narwhal.txt #list of scaffolds in Y chromosome

## 0) Noting here that the narwhal_allvariants.n60.vcf.gz file excludes two duplicate samples (94_RAHM_IQ_162 and RA_HM_IQ_121) prior to filtering. This is mentioned here because it is included in the raw sequencing data prep in de Greef et al. 2024 (https://doi.org/10.1111/gcb.17528), but was removed after SNP calling and before filtering for this current study.

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


## 7) Filter out sex-linked SNPS to create autosomal dataset

# Convert scaffold list to one-liner to use in bcftools
awk '{print $1}' $xchr_scaf_list | paste -s -d, - > xscaf_list_line
awk '{print $1}' $ychr_scaf_list | paste -s -d, - > yscaf_list_line

# Set up list
xlist=`cat xscaf_list_line`
ylist=`cat yscaf_list_line`

# Filter vcf FOR these scaffolds first
bcftools filter --regions $xlist $snps_prefix.PASS.minQ50.miss.biallel.min100kb.vcf.gz > $snps_prefix.filtered.X.vcf
bcftools filter --regions $ylist $snps_prefix.PASS.minQ50.miss.biallel.min100kb.vcf.gz > $snps_prefix.filtered.Y.vcf

# Create list of CHROM and POS for the snps in the X, also removing the ##header stuff with sed
awk '{print $1, $2}' $snps_prefix.filtered.X.vcf > X_CHROM_POS
sed '/^#/d' X_CHROM_POS > X_CHROM_POS_list
rm X_CHROM_POS

awk '{print $1, $2}' $snps_prefix.filtered.Y.vcf > Y_CHROM_POS
sed '/^#/d' Y_CHROM_POS > Y_CHROM_POS_list
rm Y_CHROM_POS

# Merge these lists for one XY list to filter out
cat X_CHROM_POS_list Y_CHROM_POS_list > XY_CHROM_POS_list


# Do the filter to remove snps on X and Y
vcftools --gzvcf $snps_prefix.PASS.minQ50.miss.biallel.min100kb.vcf.gz --exclude-positions XY_CHROM_POS_list --recode --recode-INFO-all --out $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes
mv $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.recode.vcf $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.vcf

# Can make a list of scaffolds to double check
#grep -v "^#" $snps.X.vcf | cut -f1 | sort | uniq > X_scaffold_check.txt


## 8) Remove close-kins and individuals with high missingness
vcftools --vcf $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.vcf --remove-indv ARRE_06_1164 --remove-indv B95_51_BI --remove-indv B96_393_SB --recode --recode-INFO-all --stdout > $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.vcf


## 9) Minor allele frequency filter
vcftools --vcf $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.vcf --maf 0.05 --recode --recode-INFO-all --stdout > $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.vcf


## 10) Additional missingness filter (stricter filter prior to imputing) for 10% max missingness
vcftools --vcf $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.vcf --max-missing 0.9 --recode --recode-INFO-all --stdout > $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.vcf

# Lastly, can rename file
mv $snps_prefix.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.vcf $snps_prefix.filtered.n57.vcf
