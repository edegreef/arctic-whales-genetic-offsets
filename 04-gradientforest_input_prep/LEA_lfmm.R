# Script for preparing the putative adaptive SNP datasets to use in Gradient Forest models. This was completed for narwhal and beluga whale, and was not applicable to the bowhead whale dataset (I tried the bowhead whale data in this pipeline, but due to lack of population structure and small sample size we were not able to properly identify adaptive SNPs).

# This script is divided into 3 main steps. The starting lines for each step are listed below.
# Step 1) Run genome scan with snmf.pvalues (starting on line 23)
# Step 2) Run latent factor mixed models (LFMMs) with environmental predictors (starting on line 90)
# Step 3) Prepare list of SNPs to extract to use for Gradient Forest (starting on line 145)

# Starting input files are in ped/map format for snmf (completed with "impute_and_thin_snps.sh")

# R-package "LEA" and "qvalue" need to be installed through BiocManager:
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("LEA")
#BiocManager::install("qvalue")

# Load packages
library(LEA) # (Frichot & Francois, 2015)
library(vcfR) # (Knaus & Grunwald, 2017)
library(tidyverse) #(Wickham et al., 2019)
library(qvalue) # (Storey et al., 2025)

##########  
##### Step 1): Run genome scan with SNMF
##########  

# Convert ped file to geno file (automatically outputs in working directory)
# narwhal: narwhal_snps.filtered.n57.imputed.thin.ped
# beluga: beluga_snps.filtered.n140.imputed.thin.ped

ped2geno("beluga_snps.filtered.n140.imputed.thin.ped")

# Run snmf for K 1-6, if want to save time can just use a specific K based on previous snmf run and/or fewer repetitions. From previous work (de Greef et al. 2024; Montana et al, 2024; and Muller et al. unpublished), we expect that beluga has K=6 population, narwhal has K=3 subpopulations.
project=NULL
project=snmf("beluga_snps.filtered.n140.imputed.thin.geno", K=1-6, entropy=TRUE, repetitions=10, project="new")

# Can also load snmf project if already run/saved:
#project=load.snmfProject("beluga_snps.filtered.n140.imputed.thin.snmfProject")

# Genome scan for selection: population differentiation tests. #Using K=6 for beluga, and K=3 for narwhal
p <- snmf.pvalues(project, entropy=TRUE, ploidy=2, K=6)
pvalues <- p$pvalues

# Plot histogram of pvalues
par(mfrow=c(1,1))
hist(pvalues, col="orange")
plot(-log10(pvalues), pch=19, col="blue", cex=.5)

# Look at data frame and log vals
pvals <- as.data.frame(pvalues)
pvals$log <- -log10(pvals$pvalues)

# How many with -log10 minimum 3?
pvals_log10_3 <- subset(pvals, log >= 3)

# Because not many (though comparable for beluga) decided to extract the top 1%

# Adding snp info to the p values
# Load map file
# narwhal: narwhal_snps.filtered.n57.imputed.thin.map
# beluga: beluga_snps.filtered.n140.imputed.thin.map
snp_info <- read.table("beluga_snps.filtered.n140.imputed.thin.map")

# Pull out only CHROM and POS
snp_info <- select(snp_info, V1, V4)

# Extract the log p's, add snp_info, order from large to small
logp <- as.data.frame(-log10(pvalues))
logp$original_order <- 1:nrow(logp)
snps_logp <- cbind(snp_info, logp)
colnames(snps_logp) <- c("CHROM", "POS", "logp", "order")

snps_logp2 <- snps_logp[,c(4,1,2,3)]

# Putting in order of logp (higher more significant)
snps_logp_order <- snps_logp2[order(-snps_logp2$logp),]

# Set variable for number of snps
count <- nrow(snps_logp_order)

# Extract top x%, using 1% here
top_count <- floor(0.01*count)
snp_subset <- as.data.frame(snps_logp_order[1:top_count,])

# Save list of snps - base genome scan, save top snps as .txt or .csv
write_delim(snp_subset, "beluga_K6.imputed.thin1000.top0.01logp.subset.txt")
write.csv(snp_subset, "beluga_K6.imputed.thin1000.top0.01logp.subset.csv", row.names = FALSE, quote=FALSE)


##########  
##### Step 2) Run latent factor mixed models (LFMMs) with environmental predictors
########## 

# Select best run for K from the snmf models (beluga K=6, and narwhal K=3)
best <- which.min(cross.entropy(project, K=6))

# Make lfmm file of the snps
ped2lfmm("beluga_snps.filtered.n140.imputed.thin.ped")

## Run lfmm! This is set up for doing 1 environmental variable at a time, then have to re-run this section for each variable (and each species)

# Manually made .env files for each environmental variable from the env_data_by_sample_narwhal.csv and env_data_by_sample_beluga.csv file

# Lines 105-138 need to be re-run for each environmental variable. Set the respective env file on line 107 (sst, ice, chloro, salinity), and update line 119 and 138 (and title in 141 if want to save plot) with the environmental variable name for the respective output filename.

# Set up variables
lfmm_data <- "beluga_snps.filtered.n140.imputed.thin.lfmm"
lfmm_env <- "env_files/belgua_salinity_100KM.env"

# Run lfmm2 with selected K (k=6 for beluga, k=3 for narwhal)
mod <- lfmm2(input=lfmm_data, env=lfmm_env, K=6)

# Computing P-values and plotting their minus log10 values 
pv <- lfmm2.test(object=mod, 
                 input=lfmm_data, 
                 env=lfmm_env,
                 linear=TRUE, genomic.control=TRUE)

# Optionally save workspace
save.image("beluga_snps.filtered.n140.imputed.thin_salinity_100KM.RData")

# Re-load if needed if resuming from here
#load("beluga_snps.filtered.n140.imputed.thin_salinity_100KM.RData")

# Setting up to add snp info to the p values
snp_info <- read.table("beluga_snps.filtered.n140.imputed.thin.map")

# Pull out only CHROM and POS
snp_info <- select(snp_info, V2)
colnames(snp_info) <- "snp"

pvals <- data.frame("snp"=snp_info,"zscore"=pv$zscores,"pvalue"=pv$pvalues)

# Use qvalue to use FDR rave 0.05 to mark significant snps
pvals$qval <- qvalue(pvals$pvalue, fdr.level=0.05)$signif
lfmm.snps <- pvals[pvals$qval == TRUE,"snp"]

# Not many in this dataset, so saving the whole thing instead-- need to do this for each variable
write.csv(pvals, "beluga_snps.filtered.n140.imputed.thin_salinity_100KM.pvalues.csv")

# Can also make manhattan plot
plot(-log10(pv$pvalues), col="grey", cex=.4, pch=19, main="salinity_100KM_present")


##########  
##### Step 3) Prepare list of SNPs to extract to use for Gradient Forest
##########  

# After running lfmm2 for each variable & species separately, then put results together into one dataframe.
# Not many sig snps so going with the top 1% again for each
# Last part of step 3 is to make list of unique snps (could be some overlap so want to avoid, though maybe it doesn't matter if the list for filtering later has overlap?) Need to run lines 170-189 for each environmental variable (specifically changing variable name on line 170 and 189).

# Read in top0.01 snps from genome scan:
scan <- read.csv("beluga_K6.imputed.thin1000.top0.01logp.subset.csv")

# Probably not the most efficient code, but next extracting the top 1% signiticant SNPs like we did for the snmf scan.
# One environmental variable at a time, sea surface temp, ice thickness, salinity, chlorophyll (need to run Lines 170-189 for each variable)

sst <- read.csv("beluga_snps.filtered.n140.imputed.thin_sst_100KM.pvalues.csv")
icethick <- read.csv("beluga_snps.filtered.n140.imputed.thin_icethick_100KM.pvalues.csv")
salinity <- read.csv("beluga_snps.filtered.n140.imputed.thin_salinity_100KM.pvalues.csv")
chloro <- read.csv("beluga_snps.filtered.n140.imputed.thin_chloro_100KM.pvalues.csv")

# Read in the map file to add snp info to the p values
snp_info <- read.table("beluga_snps.filtered.n140.imputed.thin.map")

# Pull out only CHROM and POS
snp_info <- select(snp_info, V1, V4)

# Get the log p's, add snp_info, order from large to small, then pull top 1%. Need to change the variable to either "sst", "icethick", "salinity", or "chloro" depending which variable is being extracted from
logp <- as.data.frame(-log10(salinity$pvalue))
logp$original_order <- 1:nrow(logp) #order number can also act as snp id later since this is for all snps
snps_logp <- cbind(snp_info, logp)
colnames(snps_logp) <- c("CHROM", "POS", "logp", "order")
snps_logp2 <- snps_logp[,c(4,1,2,3)] #reordering the columns

# Putting in order of logp (higher more significant)
snps_logp_order <- snps_logp2[order(-snps_logp2$logp),]

# Set variable for number of snps
count <- nrow(snps_logp_order)

# Extract top x%, using 1% here
top_count <- floor(0.01*count)
snp_subset <- as.data.frame(snps_logp_order[1:top_count,])

topsnps <- snp_subset

# Save it as csv too so dont have to redo later
write.csv(topsnps, "beluga_K6.imputed.thin1000.lfmm.salinity_100KM.top0.01logp.subset.csv")

#####
# Continue or reload the "..subset.csv" files for all the top 1% snps.
scan <- read.csv("beluga_K6.imputed.thin1000.top0.01logp.subset.csv")
sst <- read.csv("beluga_K6.imputed.thin1000.lfmm.sst_100KM.top0.01logp.subset.csv")
icethick <- read.csv("beluga_K6.imputed.thin1000.lfmm.icethick_100KM.top0.01logp.subset.csv")
salinity <- read.csv("beluga_K6.imputed.thin1000.lfmm.salinity_100KM.top0.01logp.subset.csv")
chloro <- read.csv("beluga_K6.imputed.thin1000.lfmm.chloro_100KM.top0.01logp.subset.csv")

# See how many snps overlap, and make a nice snp list to use for extracting these snps from vcf later
scan <- select(scan, -logp)
sst <- select(sst, c(-X,-logp))
icethick <- select(icethick, c(-X,-logp))
salinity <- select(salinity, c(-X,-logp))
chloro <- select(chloro, c(-X,-logp))

# Combine snps form all subsets
merged <- rbind(scan, sst, icethick, salinity, chloro)

# Retain only unique snp IDs
merged_nodup <- merged[!duplicated(merged), ]

# Map snp CHROM POS list for vcftools filtering
merged_nodup_order <- merged_nodup[order(merged_nodup$order),]
snp_list <- select(merged_nodup_order, -order)

# Save this list to use in vcftools filtering to extract the putative adaptive snps
write_delim(snp_list, file="beluga_K6_top0.01_scan_4env_merged_snplist.txt", delim = "\t", col_names=FALSE)
