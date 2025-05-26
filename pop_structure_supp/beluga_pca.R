# Running PCA using help from the PCAdapt vignette: https://bcm-uga.github.io/pcadapt/articles/pcadapt.html
# This is for beluga whale specifically since Muller et al. is not published yet.
# PCAs for narwhal and bowhead whale are in de Greef et al. (2024); https://doi.org/10.1111/gcb.17528

# Load modules
library(pcadapt) # (Prive et al., 2020)
library(ggplot2) # (Wickham, 2016)
library(patchwork) # (Pedersen, 2024)
library(dplyr) # (Wickham et al., 2023)

# Load beluga whale data file beluga_snps.filtered.n140.imputed.thin.vcf.gz.
# Loading as a vcf is deprecated now it seems? so need to use bfile format, which can be converted on linux through "plink --allow-extra-chr --make-bed --vcf beluga_snps.filtered.n140.imputed.thin.vcf.gz --set-missing-var-ids @:#\$1,\$2 --out beluga_snps.filtered.n140.imputed.thin --double-id"
snp_data <- read.pcadapt("beluga_snps.filtered.n140.imputed.thin.bed", type = "bed")

# Load sample info
sample_info <- read.csv("sample_info_beluga.csv", header=T)
sample_info <- subset(sample_info, final_vcf=="yes") 

# Need genetic pop info too for color labels
pop_info <- read.csv("sample_locations_all_species.csv")
pop_info <- subset(pop_info, Beluga=="Y")
pop_info <- pop_info[,c("Location_abbrev", "beluga_pop")]
sample_info <- left_join(sample_info, pop_info, by=c("location_ID"="Location_abbrev"))


# Run pcadapt. K value will be how many eigenvectors to be produced
x <- pcadapt(input = snp_data, K = 8)

# Screeplot
plot(x, option = "screeplot")

# Plot quick PCA
plot(x, option = "scores", pop = sample_info$beluga_pop)

# Plot other eigenvectors
plot(x, option = "scores", i = 3, j = 4, pop = sample_info$beluga_pop)

# Look at pca scores
scores <- as.data.frame(x$scores)

# Look at loadings
loadings <- as.data.frame(x$loadings)

# Z scores
z_scores <- as.data.frame(x$zscores)

# Look at proportion variance
proportion <- as.data.frame(x$singular.values)
proportion$squared <- proportion$`x$singular.values`* proportion$`x$singular.values`
prop_var <- as.data.frame(proportion$squared)
PC1_proportion <- (round(prop_var[1,], digits=4))*100
PC2_proportion <- (round(prop_var[2,], digits=4))*100
PC3_proportion <- (round(prop_var[3,], digits=4))*100
PC4_proportion <- (round(prop_var[4,], digits=4))*100

# Save the scores as separate data file to adjust the pca plot for colors, etc
evec <- cbind(sample_info$Vcf_ID, scores)
colnames(evec)[1] <- "sample"

# Plot in ggplot, this is part of supplemental figure S1
# PC1 and PC2
pca <- ggplot(data=evec, aes(x=V1,y=V2))+
  geom_point(aes(color=sample_info$beluga_pop), alpha=0.7, size=2.5)+  
  theme_bw()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  xlab(paste("PC1 (", PC1_proportion, "%)", sep=""))+
  ylab(paste("PC2 (", PC2_proportion, "%)", sep=""))+
  labs(color= "Population")+
  scale_colour_manual(values=c("#F0E442", "#009E73", "#0072B2", "#56B4E9", "#CC79A7", "#E69F00"),             
                      breaks=c("EHA", "CS", "HB", "LWR", "JB", "SL"))

pca

# PC3 and PC4
pca2 <- ggplot(data=evec, aes(x=V3,y=V4))+
  geom_point(aes(color=sample_info$beluga_pop), alpha=0.7, size=2.5)+  
  theme_bw()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  xlab(paste("PC3 (", PC3_proportion, "%)", sep=""))+
  ylab(paste("PC4 (", PC4_proportion, "%)", sep=""))+
  labs(color= "Population")+
  scale_colour_manual(values=c("#F0E442", "#009E73", "#0072B2", "#56B4E9", "#CC79A7", "#E69F00"),             
                      breaks=c("EHA", "CS", "HB", "LWR", "JB", "SL"))
pca2

# Combine PCA plots
pca + pca2 + plot_layout(guides = "collect") 

ggsave("beluga_pca_PC1-4.png", width=7.5, height=3, dpi=400)
