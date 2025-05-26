# Gradient Forest analyses on Arctic whale datasets

# I got lots of help for gradient forest analyses from pgugger's github (https://github.com/pgugger/LandscapeGenomics/blob/master/2019/Exercise4.md)
# And help from gradient forest vignette: (https://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf)

# This script is divided into 3 main steps and is set up for one species/dataset at a time. Parts that require a different input or output file name for different species are noted. The starting lines for each main step are listed below.

# Step 1): Initial run with gradient forest on present data (starting on line 43)
# Step 2): Create plots for biological space and geographical space (starting on line 180)
# Step 3): Estimate genetic offsets using future environmental data (starting on line 306)

# Installation prep:
# Need to make sure rtools is downloaded (https://cran.rstudio.com/bin/windows/Rtools/rtools42/rtools.html)
# Packages "extendedForest" and "gradientForest" can be installed through R-project.org:
#install.packages("extendedForest", repos="http://R-Forge.R-project.org")
#install.packages("gradientForest", repos="http://R-Forge.R-project.org")

# Load libraries
library(gradientForest) # (Ellis et al., 2012)
library(dplyr) # (Wickham et al., 2022)
library(vegan) # (Oksanen et al., 2025)
library(data.table) # (Barrett et al., 2025)
library(marmap) # (Pante & Simon-Bouhet, 2013)
library(rworldmap) # (South, 2023)
library(raster) # (Hijmans, 2025)
library(tidyverse) #(Wickham et al., 2019)
library(ggspatial) # (Dunnington et al., 2023)

# Additional ones for climate data stuff
library(sdmpredictors) # (Bosch et al., 2023)
library(sp) # (Pebesma & Bivand, 2005)
library(sf) # (Pebesma & Bivand, 2023)

# Some more to help with editing map plots
library(RStoolbox) # (Muller et al., 2024)
library(rnaturalearth) # (Massicotte et al., 2023)
library(tidyterra) # (Hernangomex et al., 2025)
library(terra) # (Hijmans et al., 2025)
library(svglite) # (Wickham et al., 2025)


##########  
##### Step 1: Initial run with gradient forest on present data
##########  
## 1.1 prepare data files

# Load snp data. Can use read.table(), but fread() is faster and easier.
# File names for each species: 
# narwhal - "narwhal_snps_filtered_K3_top0.01_scan_4env.forR"
# bowhead whale - "bowhead_RWmap_snps.filtered.n19.imputed.thin.forR"
# beluga whale - "beluga_snps_filtered_K6_top0.01_scan_4env.forR"

snp_temp <- fread("beluga_snps_filtered_K6_top0.01_scan_4env.forR", header=T)

snp_temp <- as.data.frame(snp_temp)
snp <- snp_temp[,-1]
rownames(snp) <- snp_temp[,1]

# Load environmental data for each sample
# File names for each species: 
# narwhal - "env_data_by_sample_narwhal.csv"
# bowhead whale - "env_data_by_sample_bowhead.csv"
# beluga whale - "env_data_by_sample_beluga.csv"

site.clim.dat <- read.csv("env_data_by_sample_beluga.csv")
colnames(site.clim.dat)

# Pull out sample and location IDs and present climate data
pres.clim.points <- site.clim.dat[c("Lon", "Lat", "sst","icethick", "salinity", "chloro")]

## 1.2 Generate PCNM spatial variables 

# Because sites are within water, using marmap distances instead of linear lines across land

# Coordinates #unique ones. since 
coord <- pres.clim.points[,c('Lon','Lat')]
#pcnm <- pcnm(dist(coord))  # if using linear

### Adding in marmap distances

# Download ocean depth map. Lower resolution number is more fine-scale.
ocean_map <- getNOAA.bathy(lon1 = -105, lon2 = -50, lat1 = 45, lat2 = 85, resolution = 5)
colnames(coord) <- c("x", "y")
trans <- trans.mat(ocean_map, min.depth = 0)
dist <- lc.dist(trans,coord,res="dist")
dist.mat <- as.matrix(dist)
pcnm <- pcnm(dist.mat)

###
# Keep half of positive pcnms
keep <- round(length(which(pcnm$value > 0))/2) 
pcnm.keep <- scores(pcnm)[,1:keep]
pcnm.keep

# Plot the PCNM data

# Add coord and pcnm.keep together
pcnm.for.plot <- cbind(coord, pcnm.keep)
map_outline <- getMap(resolution="high")

# Create cropping boundary for eastern Canadian Arctic
crop.extent <- extent(-110,-45, 44, 80)

# Crop outline to boundary and convert to dataframe
map_outline <- crop(map_outline, y=crop.extent) %>% fortify()

# Plot PCNMs on map- have to go one at a time here (adjust the "PCNM#" in the "fill=" on line 110). This corresponds to supplemental figures S7, S8, S9
ggplot()+
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray80", colour="gray50", linewidth=0.2)+
  geom_point(data=pcnm.for.plot, aes(x=x, y=y, fill=PCNM1), shape=21, stroke=0.7,color="black",size=5)+
  scale_fill_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027")+
  scale_x_continuous(expand=c(0,0), breaks=c(-100,-80,-60)) +
  scale_y_continuous(expand=c(0,0), breaks=c(50,60,70,80)) +
  xlab("Longitude")+
  ylab("Latitude")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(text=element_text(size=15),legend.title = element_text( size=12), legend.text=element_text(size=12), panel.border=element_rect(colour="black", fill=NA))

# Save pcnm plot
ggsave("beluga_pcnm1.png", width=4.5, height=4, dpi=400)

# Create file with climate and PCNM variables (no lat/long)
env.gf <- cbind(pres.clim.points[,c("sst", "icethick", "salinity", "chloro")], pcnm.keep)

# Define max number of splits
lev <- floor(log2(0.368*nrow(env.gf)/2))
lev

# Run GF using the climate, spatial, and snp data inputs. Used different ntree and nbin values for each species due to dataset size. Higher nbins reduces run time. If set too low for a large dataset (e.g., bowhead whale) then will crash.
gf <- gradientForest(cbind(env.gf, snp), 
                     predictor.vars=colnames(env.gf),
                     response.vars=colnames(snp), 
                     ntree=300, #500 for narwhal, 100 for bowhead, 300 for beluga
                     maxLevel=lev, 
                     trace=T, 
                     corr.threshold=0.50,
                     nbin=301) #301 for narwhal, 601 for bowhead, 301 for beluga

# Optionally save RData here so don't have to re-run gradientForest()
#save.image("beluga_K6_GF_snps_scan_4env_ntree300_nbin301_EDrun.RData")


## 1.3 Initial gradient forest plots

# Bar graphs of predictor overall variable importance
# Chose colors to highlight PCNMs and known environmental variables separately (orange for PCNM, blue for variables)

# for narwhal
plot(gf, plot.type = "O", col=c("#80ace1","#80ace1","#80ace1","#80ace1","#f6a565", "#f6a565","#f6a565","#f6a565"), lwd=1.5, cex.axis=0.9)

# for bowhead
plot(gf, plot.type = "O", col=c("#80ace1","#80ace1","#f6a565","#f6a565","#f6a565", "#80ace1","#80ace1","#f6a565","#f6a565"), lwd=1.5, cex.axis=0.9)

# for beluga
plot(gf, plot.type = "O", col=c("#80ace1","#80ace1", "#80ace1","#f6a565","#f6a565","#f6a565", "#f6a565","#80ace1","#f6a565","#f6a565"), lwd=1.5, cex.axis=0.9)

# Save these as svg through the "export" button on plot, these become panel A for supplemental figures S10, S11, S12.

# Organize variables by importance, for other plots
by.importance <- names(importance(gf))
by.importance

# List the values for R2 weighted importance
as.data.frame(importance(gf, type = "Weighted"))

# Plot turnover functions, showing how allelic composition changes along spatial or environmental gradients. this is a predictor cumulative plot
plot(gf, plot.type="C", imp.vars=by.importance, 
     show.species=F, common.scale=T, cex.axis=1, cex.lab=1.2, line.ylab=1, lwd=2,col="#7393B3",
     par.args=list(mgp=c(1.5,0.5,0), mar=c(2.5,2,2,2), omi=c(0.2,0.3,0.2,0.4)))
# can't get it to plot 4x2 instead of 2x4 but can adjust in inkscape
# Save these as svg through the "export" button on plot, these become panel B for supplemental figures S10, S11, S12.


# Plot R2 measure of the fit of the random forest model 
plot(gf, plot.type="P", show.names=T, horizontal=F, cex.axis=1, cex.labels=0.7, line=2.5)

##########  
##### Step 2: Create plots for biological space and geographical space
##########  

## 2.1 Load climate data rasters
# As noted in "extract_climate_data.R", instead of re-loading the layers from online with sdmpredictors every time I run this during load_layers(), I  downloaded and saved the rasters into a local directory so can load faster from there

# Example for setting the datadir for sdmpredictors to find the downloaded rasters, otherwise can skip this line and re-load through load_layers()
#options(sdmpredictors_datadir="C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/gradient_forest_current/climate_rasters")

# Names of present variables of interest
present.vars <- c("BO22_tempmean_ss",
                  "BO22_icethickmean_ss",
                  "BO22_salinitymean_ss", 
                  "BO22_chlomean_ss")

# Load future vars now too though won't use it until later. Here, need to adjust for target RCP (RCP45, RCP60, RCP85)
future.vars <- c("BO22_RCP45_2100_tempmean_ss",
                 "BO22_RCP45_2100_icethickmean_ss",
                 "BO22_RCP45_2100_salinitymean_ss") 
future.vars.chlo <- c("BO22_RCP45_2100_chlomean_ss") #chlorophyll has different extent so importing separately

# Load rasters
present.rasters <- load_layers(present.vars)
future.rasters.nochlo <- load_layers(future.vars)
future.rasters.chlo <- load_layers(future.vars.chlo)

# Define the eastern Canadian Arctic boundary
crop.extent <- extent(-110,-45, 44, 80)

# Crop rasters to the extent
present.rasters.crop <- crop(present.rasters, crop.extent)
future.rasters.nochlo.crop <- crop(future.rasters.nochlo, crop.extent)
future.rasters.chlo.crop <- crop(future.rasters.chlo, crop.extent)

# Add the future chlorophyll with other future variables
future.rasters.crop <- stack(future.rasters.nochlo.crop, future.rasters.chlo.crop)

# Rename variables (need to match climate data points)
names(present.rasters.crop) <- c("sst", "icethick", "salinity", "chloro")
names(future.rasters.crop) <- c("sst", "icethick", "salinity", "chloro")

## 2.2 Crop to desired area, using species' range here. 
# Note: to skip cropping by species range, skip lines 224-229, and instead run lines 231-233.

# Load species range shapefile (obtained from IUCN red list database; also described in "combine_range_maps.R")
species <- read_sf('data_0.shp')

# Crop present and future rasters to species range shape file
present.rasters.range <- mask(present.rasters.crop, species)
future.rasters.range <- mask(future.rasters.crop, species)

# If not cropping, just assign the .crop variable to .range name to avoid the need to change variable name downstream:
#present.rasters.range <- present.rasters.crop
#future.rasters.range <- future.rasters.crop

# Quick plot just to check that it worked
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))
plot(present.rasters.range,col=my.colors(1000),axes=FALSE, box=FALSE)
plot(future.rasters.range,col=my.colors(1000),axes=FALSE, box=FALSE)


## 2.3 Extract present climate data from rasters and transform environmental predictors
clim.land <- raster::extract(present.rasters.range, 1:ncell(present.rasters.range), df = TRUE)
clim.land <- na.omit(clim.land)

# Use predict function to transform environmental predictors
pred <- predict(gf, clim.land[,-1])  #note the removal of the cell ID column with [,-1])

# Convert predictions to color scale for map (following https://github.com/pgugger/LandscapeGenomics/blob/master/2019/Exercise4.md)
# Use pca on predictions
pca <- prcomp(pred, center=T, scale.=F)

# Assign pcs to colors
r <- pca$x[,1]
g <- pca$x[,2]
b <- pca$x[,3]

# Scale colors
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255

# Define raster properties with an existing one
mask<-present.rasters.crop$salinity 
#mask[]<-as.numeric(mask[]>0)

# Assign color to raster
rastR <- rastG <- rastB <- mask
rastR[clim.land$ID] <- r
rastG[clim.land$ID] <- g
rastB[clim.land$ID] <- b

# Stack color raster
rgb.rast <- stack(rastR, rastG, rastB)

####

# Plot and save the gradient forest map (these are not offsets yet). These become the base of supplemental figures S13, S14, S15.
svglite("beluga_GF_Map_updated_rangeclip.svg", width=6, height=6)
plotRGB(rgb.rast, colNA="gray80")
points(pres.clim.points$Lon, pres.clim.points$Lat, pch=19, col="black")
dev.off()

# Setting up PCA legend of biological space
nvs <- dim(pca$rotation)[1]
vec <- c("sst", "icethick", "salinity", "chloro")
lv <- length(vec)
vind <- rownames(pca$rotation) %in% vec
scal <- 40
xrng <- range(pca$x[,1], pca$rotation[,1]/scal) * 1.1
yrng <- range(pca$x[,2], pca$rotation[,2]/scal) * 1.1

# Plot and save the PCA legend associated with the gradient forest map
svglite("narwhal_GF_Map_updated_rangeclip_PCA_legend.svg", width=6, height=5)
plot((pca$x[,1:2]), xlim=xrng, ylim=yrng, 
     pch=".", cex=4, col=rgb(r,g,b, max=255),
     asp=1)
arrows(rep(0, lv), rep(0, lv), pca$rotation[vec,1]/scal, pca$rotation[vec,2]/scal, length=0.0625)
jit <- 0.0015
text(pca$rotation[vec,1]/scal + jit * sign(pca$rotation[vec,1]), pca$rotation[vec,2]/scal + jit * sign(pca$rotation[vec,2]), labels=vec)
dev.off()

# Helpful description from the gradientForest vignette:
#"The colors represent genetic variation (allelic composition) as predicted based on the modeled relationships with environmental and spatial variables. Similar colors are predicted to be more similar genetically."

##########  
##### Step 3: Estimate genetic offsets using future environmental data
##########

## 3.1 prepare data files
clim.land.future <- raster::extract(future.rasters.range, 1:ncell(future.rasters.range), df = TRUE)
clim.land.future <- na.omit(clim.land.future)

# Transform environmental variables
pred.future <- predict(gf, clim.land.future[,-1])

# Estimate genetic offset!
genetic.offset.adaptive <- sqrt((pred.future[,1]-pred[,1])^2 + 
                                  (pred.future[,2]-pred[,2])^2 + 
                                  (pred.future[,3]-pred[,3])^2 + 
                                  (pred.future[,4]-pred[,4])^2)

# Define raster properties --- the variable doesn't matter here i think
rast.offset <- future.rasters.range$salinity 

# Assign genetic offset values (difference between future and present predictions) to raster 
rast.offset[clim.land.future$ID] <- genetic.offset.adaptive

# Make color scale
offset.colors = colorRampPalette(c("#4575B4", "#abd9e9","#ffffbf","#fdae61","#f46d43","#a50026"))

# Quick plot
pdf("beluga_GF_GeneticOffset_RCP85_2100_rangeclip.pdf")
plot(rast.offset, col=offset.colors(200))
points(pres.clim.points$Lon, pres.clim.points$Lat)
dev.off()

# Add land to map
# Get map outlines from rworldmap package. For "high" resolution, need to also install package rworldxtra
map_outline <- getMap(resolution="high")

# Crop outline to boundary and convert to dataframe
map_outline <- crop(map_outline, y=crop.extent) %>% fortify()

# Create SpatRaster to work with tidyterra
rast.offset.spat <- rast(rast.offset)

# Also upload coordinates by site for map site points:
sites <- read.csv("sample_locations_all_species.csv", header=T)

# Filter coords by species, and adjust lines 351-353, and 361-371 depending which species is being run
#sites <- subset(sites, Narwhal>0) #narwhal
#sites <- subset(sites, Bowhead>0) #bowhead
sites <- subset(sites, Beluga=="Y") #beluga

# Plot and save as svg - run this part separately for each species map. Species-specific lines are labeled. This is associated with main figure 2, 3, 4, and supplemental figures S16-S21.
svglite("beluga_GF_GeneticOffset_RCP45_2100_landmap_sites_step_updated2_shape_ECA.svg", width=6, height=6)
ggplot()+
  geom_spatraster(data=rast.offset.spat)+
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray80", colour="gray50", linewidth=0.2)+
  
  # beluga points
  geom_point(data=sites, aes(x=est_longitude, y=est_latitude, shape=beluga_pop),size=2,colour="black")+ 
  scale_shape_manual(values=c(15, 8, 17, 18, 12, 16))+ 
  
  # narwhal points
 # geom_point(data=sites, aes(x=est_longitude, y=est_latitude, shape=narwhal_pop),size=2,colour="black")+
 # scale_shape_manual(values=c(17, 15, 18))+
  
  # bowhead points
  #geom_point(data=sites, aes(x=est_longitude, y=est_latitude, shape=bowhead_pop),size=2,colour="black")+
  #scale_shape_manual(values=c(16))+
  
  scale_fill_stepsn(n.breaks = 20, colours = offset.colors(12),  na.value="white")+ #limits=c(0,0.11))+
  scale_x_continuous(expand=c(0,0), breaks=c(-100,-80,-60)) +
  scale_y_continuous(expand=c(0,0), breaks=c(50,60,70,80)) +
  xlab("Longitude")+
  ylab("Latitude")+
  labs(fill="Genetic offset")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.title = element_text(hjust = 0.5))+
  annotation_north_arrow(height = unit(1, "cm"),width = unit(1, "cm"),
                         location = "tr", which_north = "true",
                         pad_x = unit(0.12, "cm"), pad_y = unit(0.3, "cm"), 
                         style = ggspatial::north_arrow_fancy_orienteering())+
  theme(text=element_text(size=15))+
  theme(legend.title = element_text( size=12), legend.text=element_text(size=12))
dev.off()


# Save the rast.offset for later too
raster::writeRaster(rast.offset, "raster.offset.RCP45.2100.beluga.easternCanada.tif")

# Note about how to import raster back in with "raster"
#rast.offset.test <- raster("raster.offset.RCP45.2100.beluga.easternCanada.tif")
