#!/bin/bash

#SBATCH --time=4-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=downsample
#SBATCH --output=%x-%j.out

# Downsampling select samples with GATK (McKenna et al. 2010)

# Load gatk
module load nixpkgs/16.09 gatk/4.1.2.0

# 5 samples keeping 1/3 of the reads, and 70 samples keeping 1/2 of the reads. 
keep1=0.3333
keep2=0.5

# Could do all in one script, or split in batches if want to do faster.
# The samples to keep 1/3
gatk --java-options "-Xmx50G" DownsampleSam -I CS_Pagnirtung_1022_1_R1.fastq.gz.deDupRG.bam -O CS_Pagnirtung_1022_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep1
gatk --java-options "-Xmx50G" DownsampleSam -I CS_Pagnirtung_1025_1_R1.fastq.gz.deDupRG.bam -O CS_Pagnirtung_1025_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep1
gatk --java-options "-Xmx50G" DownsampleSam -I CS_Pagnirtung_1549_1_R1.fastq.gz.deDupRG.bam -O CS_Pagnirtung_1549_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep1
gatk --java-options "-Xmx50G" DownsampleSam -I CS_Pagnirtung_1553_1_R1.fastq.gz.deDupRG.bam -O CS_Pagnirtung_1553_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep1
gatk --java-options "-Xmx50G" DownsampleSam -I CS_Pagnirtung_1558_1_R1.fastq.gz.deDupRG.bam -O CS_Pagnirtung_1558_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep1

# The samples to keep 1/2
# 8 samples from CS
gatk --java-options "-Xmx50G" DownsampleSam -I CS_Pangnirtung_12_1_R1.fastq.gz.deDupRG.bam -O CS_Pangnirtung_12_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I CS_Pangnirtung_14_1_R1.fastq.gz.deDupRG.bam -O CS_Pangnirtung_14_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I CS_Pangnirtung_15_1_R1.fastq.gz.deDupRG.bam -O CS_Pangnirtung_15_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I CS_Pangnirtung_43_1_R1.fastq.gz.deDupRG.bam -O CS_Pangnirtung_43_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I CS_Pangnirtung_44_1_R1.fastq.gz.deDupRG.bam -O CS_Pangnirtung_44_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I CS_Pangnirtung_45_1_R1.fastq.gz.deDupRG.bam -O CS_Pangnirtung_45_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I CS_Pangnirtung_46_1_R1.fastq.gz.deDupRG.bam -O CS_Pangnirtung_46_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I CS_Pangnirtung_60_3_R1.fastq.gz.deDupRG.bam -O CS_Pangnirtung_60_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2

# 14 samples from EHA
gatk --java-options "-Xmx50G" DownsampleSam -I EHA_Cunningham_Inlet_101_3.fastq.gz.deDupRG.bam -O EHA_Cunningham_Inlet_101_3.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I EHA_Cunningham_Inlet_103_3.fastq.gz.deDupRG.bam -O EHA_Cunningham_Inlet_103_3.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I EHA_Cunningham_Inlet_104_3.fastq.gz.deDupRG.bam -O EHA_Cunningham_Inlet_104_3.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I EHA_Cunningham_Inlet_105_3.fastq.gz.deDupRG.bam -O EHA_Cunningham_Inlet_105_3.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I EHA_Cunningham_Inlet_106_3.fastq.gz.deDupRG.bam -O EHA_Cunningham_Inlet_106_3.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I EHA_Cunningham_Inlet_107_3.fastq.gz.deDupRG.bam -O EHA_Cunningham_Inlet_107_3.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I EHA_Cunningham_Inlet_108_3.fastq.gz.deDupRG.bam -O EHA_Cunningham_Inlet_108_3.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I EHA_Cunningham_Inlet_109_3.fastq.gz.deDupRG.bam -O EHA_Cunningham_Inlet_109_3.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I EHA_Cunningham_Inlet_110_3.fastq.gz.deDupRG.bam -O EHA_Cunningham_Inlet_110_3.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I EHA_Cunningham_Inlet_111_3.fastq.gz.deDupRG.bam -O EHA_Cunningham_Inlet_111_3.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I EHA_Cunningham_Inlet_63_3.fastq.gz.deDupRG.bam -O EHA_Cunningham_Inlet_63_3.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I EHA_Cunningham_Inlet_65_3.fastq.gz.deDupRG.bam -O EHA_Cunningham_Inlet_65_3.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I EHA_Cunningham_Inlet_94_3.fastq.gz.deDupRG.bam -O EHA_Cunningham_Inlet_94_3.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I EHA_Cunningham_Inlet_95_3.fastq.gz.deDupRG.bam -O EHA_Cunningham_Inlet_95_3.fastq.gz.deDupRG.downsampled.bam -P $keep2

# 15 samples from HB
gatk --java-options "-Xmx50G" DownsampleSam -I HB_Sanikiluaq_25_1_R1.fastq.gz.deDupRG.bam -O HB_Sanikiluaq_25_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I HB_Sanikiluaq_26_1_R1.fastq.gz.deDupRG.bam -O HB_Sanikiluaq_26_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I HB_Sanikiluaq_27_1_R1.fastq.gz.deDupRG.bam -O HB_Sanikiluaq_27_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I HB_Sanikiluaq_28_1_R1.fastq.gz.deDupRG.bam -O HB_Sanikiluaq_28_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I HB_Sanikiluaq_29_1_R1.fastq.gz.deDupRG.bam -O HB_Sanikiluaq_29_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I HB_Sanikiluaq_30_1_R1.fastq.gz.deDupRG.bam -O HB_Sanikiluaq_30_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I HB_Sanikiluaq_31_1_R1.fastq.gz.deDupRG.bam -O HB_Sanikiluaq_31_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I HB_Sanikiluaq_32_1_R1.fastq.gz.deDupRG.bam -O HB_Sanikiluaq_32_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I HB_Sanikiluaq_33_1_R1.fastq.gz.deDupRG.bam -O HB_Sanikiluaq_33_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I HB_Sanikiluaq_34_1_R1.fastq.gz.deDupRG.bam -O HB_Sanikiluaq_34_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I HB_Sanikiluaq_35_1_R1.fastq.gz.deDupRG.bam -O HB_Sanikiluaq_35_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I HB_Sanikiluaq_36_1_R1.fastq.gz.deDupRG.bam -O HB_Sanikiluaq_36_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I HB_Sanikiluaq_37_1_R1.fastq.gz.deDupRG.bam -O HB_Sanikiluaq_37_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I HB_Sanikiluaq_38_1_R1.fastq.gz.deDupRG.bam -O HB_Sanikiluaq_38_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I HB_Sanikiluaq_39_1_R1.fastq.gz.deDupRG.bam -O HB_Sanikiluaq_39_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2

# 10 samples from Iqualuit
gatk --java-options "-Xmx50G" DownsampleSam -I NA_Iqaluit_58_3_R1.fastq.gz.deDupRG.bam -O NA_Iqaluit_58_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I NA_Iqaluit_64_3_R1.fastq.gz.deDupRG.bam -O NA_Iqaluit_64_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I NA_Iqaluit_77_3_R1.fastq.gz.deDupRG.bam -O NA_Iqaluit_77_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I NA_Iqaluit_78_3_R1.fastq.gz.deDupRG.bam -O NA_Iqaluit_78_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I NA_Iqaluit_79_3_R1.fastq.gz.deDupRG.bam -O NA_Iqaluit_79_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I NA_Iqaluit_80_3_R1.fastq.gz.deDupRG.bam -O NA_Iqaluit_80_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I NA_Iqaluit_81_3_R1.fastq.gz.deDupRG.bam -O NA_Iqaluit_81_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I NA_Iqaluit_82_3_R1.fastq.gz.deDupRG.bam -O NA_Iqaluit_82_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I NA_Iqaluit_83_3_R1.fastq.gz.deDupRG.bam -O NA_Iqaluit_83_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I NA_Iqaluit_84_3_R1.fastq.gz.deDupRG.bam -O NA_Iqaluit_84_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2

# 9 samples from SL
gatk --java-options "-Xmx50G" DownsampleSam -I SL_NA_113_3_R1.fastq.gz.deDupRG.bam -O SL_NA_113_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I SL_NA_114_3_R1.fastq.gz.deDupRG.bam -O SL_NA_114_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I SL_NA_49_3_R1.fastq.gz.deDupRG.bam -O SL_NA_49_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I SL_NA_50_3_R1.fastq.gz.deDupRG.bam -O SL_NA_50_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I SL_NA_51_3_R1.fastq.gz.deDupRG.bam -O SL_NA_51_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I SL_NA_53_3_R1.fastq.gz.deDupRG.bam -O SL_NA_53_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I SL_NA_54_3_R1.fastq.gz.deDupRG.bam -O SL_NA_54_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I SL_NA_56_3_R1.fastq.gz.deDupRG.bam -O SL_NA_56_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I SL_NA_57_3_R1.fastq.gz.deDupRG.bam -O SL_NA_57_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2

# 14 samples from WHB
gatk --java-options "-Xmx50G" DownsampleSam -I WHB_Arviat_10_1_R1.fastq.gz.deDupRG.bam -O WHB_Arviat_10_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I WHB_Arviat_11_1_R1.fastq.gz.deDupRG.bam -O WHB_Arviat_11_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I WHB_Arviat_17_1_R1.fastq.gz.deDupRG.bam -O WHB_Arviat_17_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I WHB_Arviat_18_1_R1.fastq.gz.deDupRG.bam -O WHB_Arviat_18_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I WHB_Arviat_19_1_R1.fastq.gz.deDupRG.bam -O WHB_Arviat_19_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I WHB_Arviat_20_1_R1.fastq.gz.deDupRG.bam -O WHB_Arviat_20_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I WHB_Arviat_21_1_R1.fastq.gz.deDupRG.bam -O WHB_Arviat_21_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I WHB_Arviat_22_1_R1.fastq.gz.deDupRG.bam -O WHB_Arviat_22_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I WHB_Arviat_23_1_R1.fastq.gz.deDupRG.bam -O WHB_Arviat_23_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I WHB_Arviat_24_1_R1.fastq.gz.deDupRG.bam -O WHB_Arviat_24_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I WHB_Arviat_40_1_R1.fastq.gz.deDupRG.bam -O WHB_Arviat_40_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I WHB_Arviat_41_1_R1.fastq.gz.deDupRG.bam -O WHB_Arviat_41_1_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I WHB_Arviat_97_3_R1.fastq.gz.deDupRG.bam -O WHB_Arviat_97_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2
gatk --java-options "-Xmx50G" DownsampleSam -I WHB_Arviat_98_3_R1.fastq.gz.deDupRG.bam -O WHB_Arviat_98_3_R1.fastq.gz.deDupRG.downsampled.bam -P $keep2


# Index downsampled bams:
ls *.downsampled.bam | sed 's/.downsampled.bam$//' | parallel --jobs 4 'samtools index {}.downsampled.bam'
