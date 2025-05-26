#!/bin/bash

# Converting vcf to format used for Gradient Forest, the code was obtained from https://github.com/pgugger/LandscapeGenomics/blob/master/2019/Exercise4.md and then prepared as a bash script here.

# to use: ./vcf2forR.sh <snps prefix>

vcftools --vcf $1.vcf --012 --out $1
cut -f2- $1.012 | sed 's/-1/NA/g' >$1.temp
tr -d '\t' <$1.012.pos | tr '\n' '\t' | sed 's/[[:space:]]*$//' >header
paste <(echo "ID" | cat - $1.012.indv) <(echo "" | cat header - $1.temp) > $1.forR
rm header $1.temp
