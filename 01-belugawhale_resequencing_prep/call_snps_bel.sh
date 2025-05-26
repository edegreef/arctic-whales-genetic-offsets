#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --nodes=1
#SBATCH --job-name=platypus
#SBATCH --output=%x-%j.out

# Script to call variants with Platypus program (Rimmer et al., 2014)

# Load modules
module load nixpkgs/16.09 gcc/7.3.0 platypus/0.8.1

# Run platypus to call variants from bams & reference genome 
python /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/gcc7.3/platypus/0.8.1/bin/Platypus.py callVariants \
--bamFiles=beluga_bam_list.txt \
--refFile=S_20_00703_MOBELS_0202.22xii22.FINAL_scaffolds.fsa \
--output=146Beluga_raw.vcf \
--nCPU=32 --minReads=2
