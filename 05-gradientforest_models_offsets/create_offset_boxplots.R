# Creating box plots of the genetic offset results from "gradientforest_whales_pipe.R"
# This is to extract genetic offset values by site and then summarize by (sub)population.

# It is important to note that the raster (.tif) files used in this script were not restricted by species range shape files. Range distribution is an approximation and a few sampling sites were outside of the species range. See line 222 in "gradientforest_whales_pipe.R" or README.md file for more info about not cropping to species range.

# Load libraries
library(tidyverse) #(Wickham et al., 2019)
library(raster) # (Hijmans, 2025)
library(plyr) # (Wickham, 2011)
library(reshape2) # (Wickham, 2007)

# Import saved raster back in with raster()
rast.offset <- raster("raster.offset.RCP85.2100.beluga.easternCanada.tif")

# Load site locations, subset separately by species
sites <- read.csv("sample_locations_all_species.csv", header=T)
#sites <- subset(sites, Narwhal>0) #narwhal
#sites <- subset(sites, Bowhead>0) #bowhead
sites <- subset(sites, Beluga=="Y") #beluga

x <- sites$est_longitude
y <- sites$est_latitude

pts <- cbind(x, y) %>% as_tibble()
pts <- SpatialPoints(pts, proj4string = CRS("+proj=longlat +datum=WGS84"))

# Make sure projection matches
projection(pts)==projection(rast.offset)

# Create dataframe to store extracted data
offset.data <- tibble(ID = 1:nrow(pts@coords), Lon = pts$x, Lat = pts$y)

# Extract for each layer (buffer is in meters) 
# use fun=NULL to avoid summarizing
store_data = list()
for (i in 1:nlayers(rast.offset)){
  store_data[[i]] = raster::extract(rast.offset[[i]], pts, buffer=100000, fun=NULL)
}
# Now we have genetic offset values for 100km for each site.

# The rest of the script is divided in sections by species due to having different set of site locations. Narwhal (lines 43-101), bowhead whale (lines 104-159), beluga whale (lines 162-286)

########## Narwhal
# Make dataframe for each site
# 1=AB, 2=BI, 3=CR, 4=GF, 5=IG, 6=PG, 7=PB, 8=PI, 9=RB, 10=RE, 11=SB.
# Make colnames CHA, BI, NHB for subgroups

AB <- as.data.frame(store_data[[1]][[1]])
colnames(AB) <- "AB"
BI <- as.data.frame(store_data[[1]][[2]])
colnames(BI) <- "BI"
CR <- as.data.frame(store_data[[1]][[3]])
colnames(CR) <- "CR"
GF <- as.data.frame(store_data[[1]][[4]])
colnames(GF) <- "GF"
IG <- as.data.frame(store_data[[1]][[5]])
colnames(IG) <- "IG"
PG <- as.data.frame(store_data[[1]][[6]])
colnames(PG) <- "PG"
PB <- as.data.frame(store_data[[1]][[7]])
colnames(PB) <- "PB"
PI <- as.data.frame(store_data[[1]][[8]])
colnames(PI) <- "PI"
RB <- as.data.frame(store_data[[1]][[9]])
colnames(RB) <- "RB"
RE <- as.data.frame(store_data[[1]][[10]])
colnames(RE) <- "RE"
SB <- as.data.frame(store_data[[1]][[11]])
colnames(SB) <- "SB"

# Combining data frames with different number of observations. Some sites might have NAs, but because the 100km buffer might extend over land. Use rbind.fill to merge.

CHApop <- rbind.fill(GF, RE, SB)
CHApop$pop <- "CHA"

BIpop <- rbind.fill(AB, BI, CR, IG, PG, PB, PI)
BIpop$pop <- "BI"

NHBpop <- rbind.fill(RB)
NHBpop$pop <- "NHB"

combined <- rbind.fill(CHApop, BIpop, NHBpop)
meltData <- melt(combined)

# Quick boxplot
boxplot(data=meltData, value~pop)

# Nicer boxplot. This is associated with main figure 3, and supplemental figures S18, S19.
narwhal_offset <- meltData %>%
  mutate(pop = fct_reorder(pop, value, .fun="mean")) %>%
  ggplot(aes(x=pop, y=value, fill=pop)) + 
  geom_boxplot(width=0.5, lwd=0.5, color="gray30", fill="#fdf391ff",alpha=1, outlier.size = 0.5) +
  theme(legend.position="none")+
  ylab("Genetic offset")+
  theme_bw()+
  theme(legend.position="none", panel.grid.major.x=element_blank())+
  xlab("Population")

narwhal_offset
#write.csv(meltData, "narwhal_meltdata_for_boxplot_RCP85_2100.csv")
ggsave("narwhal_boxplot_RCP85_2100_bypop.svg", width=3, height=4, dpi=1000)


########## Bowhead
# 1=AB, 2=CD, 3=CR, 4=CH, 5=GH, 6=HB, 7=IQ, 8=PG, 9=PB, 10=PI, 11=RI, 12=RB, 13=SB, 14=WB, 15=DIS 

AB <- as.data.frame(store_data[[1]][[1]])
colnames(AB) <- "AB"
CD <- as.data.frame(store_data[[1]][[2]])
colnames(CD) <- "CD"
CR <- as.data.frame(store_data[[1]][[3]])
colnames(CR) <- "CR"
CH <- as.data.frame(store_data[[1]][[4]])
colnames(CH) <- "CH"
GH <- as.data.frame(store_data[[1]][[5]])
colnames(GH) <- "GH"
HB <- as.data.frame(store_data[[1]][[6]])
colnames(HB) <- "HB"
IQ <- as.data.frame(store_data[[1]][[7]])
colnames(IQ) <- "IQ"
PG <- as.data.frame(store_data[[1]][[8]])
colnames(PG) <- "PG"
PB <- as.data.frame(store_data[[1]][[9]])
colnames(PB) <- "PB"
PI <- as.data.frame(store_data[[1]][[10]])
colnames(PI) <- "PI"
RI <- as.data.frame(store_data[[1]][[11]])
colnames(RI) <- "RI"
RB <- as.data.frame(store_data[[1]][[12]])
colnames(RB) <- "RB"
SB <- as.data.frame(store_data[[1]][[13]])
colnames(SB) <- "SB"
WB <- as.data.frame(store_data[[1]][[14]])
colnames(WB) <- "WB"
DIS <- as.data.frame(store_data[[1]][[15]])
colnames(DIS) <- "DIS"


# Combine and set as "ECA" - all in same pop
combined <- rbind.fill(AB, CD, CR, CH, GH, HB, IQ, PG, PB, PI, RI, RB, SB, WB, DIS)
combined$pop <- "ECA"
meltData <- melt(combined)

# Quick boxplot
boxplot(data=meltData, value~variable)

# Nicer boxplot. This is associated with main figure 4, and supplemental figures S20, S21.
bowhead_offset <- meltData %>%
  ggplot(aes(x=pop, y=value, fill=pop)) + 
  geom_boxplot(width=0.5, lwd=0.5, color="gray30", fill="#81c3f4ff",alpha=1, outlier.size = 0.5) +
  theme(legend.position="none")+
  ylab("Genetic offset")+
  theme_bw()+
  theme(legend.position="none", panel.grid.major.x=element_blank())+
  xlab("Population")

bowhead_offset
#write.csv(meltData, "bowhead_meltdata_for_boxplot_RCP85_2100.csv")
ggsave("bowhead_boxplot_RCP85_2100_bypop.svg", width=3, height=4, dpi=1000)


########## Beluga whale
# 1=B1, 2=B2, 3=B3.... 16=B16

B1 <- as.data.frame(store_data[[1]][[1]])
B1$location <- "B1"
B1$site <- "B1"
colnames(B1)[1] <- "data"

B2 <- as.data.frame(store_data[[1]][[2]])
B2$location <- "B2"
B2$site <- "B2"
colnames(B2)[1] <- "data"

B3 <- as.data.frame(store_data[[1]][[3]])
B3$location <- "B3"
B3$site <- "B3"
colnames(B3)[1] <- "data"

B4 <- as.data.frame(store_data[[1]][[4]])
B4$location <- "B4"
B4$site <- "B4"
colnames(B4)[1] <- "data"

B5 <- as.data.frame(store_data[[1]][[5]])
B5$location <- "B5"
B5$site <- "B5"
colnames(B5)[1] <- "data"

B6 <- as.data.frame(store_data[[1]][[6]])
B6$location <- "B6"
B6$site <- "B6"
colnames(B6)[1] <- "data"

B7 <- as.data.frame(store_data[[1]][[7]])
B7$location <- "B7"
B7$site <- "B7"
colnames(B7)[1] <- "data"

B8 <- as.data.frame(store_data[[1]][[8]])
B8$location <- "B8"
B8$site <- "B8"
colnames(B8)[1] <- "data"

B9 <- as.data.frame(store_data[[1]][[9]])
B9$location <- "B9"
B9$site <- "B9"
colnames(B9)[1] <- "data"

B10 <- as.data.frame(store_data[[1]][[10]])
B10$location <- "B10"
B10$site <- "B10"
colnames(B10)[1] <- "data"

B11 <- as.data.frame(store_data[[1]][[11]])
B11$location <- "B11"
B11$site <- "B11"
colnames(B11)[1] <- "data"

B12 <- as.data.frame(store_data[[1]][[12]])
B12$location <- "B12"
B12$site <- "B12"
colnames(B12)[1] <- "data"

B13 <- as.data.frame(store_data[[1]][[13]])
B13$location <- "B13"
B13$site <- "B13"
colnames(B13)[1] <- "data"

B14 <- as.data.frame(store_data[[1]][[14]])
B14$location <- "B14"
B14$site <- "B14"
colnames(B14)[1] <- "data"

B15 <- as.data.frame(store_data[[1]][[15]])
B15$location <- "B15"
B15$site <- "B15"
colnames(B15)[1] <- "data"

B16 <- as.data.frame(store_data[[1]][[16]])
B16$location <- "B16"
B16$site <- "B16"
colnames(B16)[1] <- "data"

#####
## See the "beluga_pop" column in the "sites" dataframe to link pop label with site ID.

# Combine values so it's by population instead of sampling site
CSpop <- rbind.fill(B1)
CSpop$pop <- "CS"

EHApop <- rbind.fill(B2)
EHApop$pop <- "EHA"

HBpop <- rbind.fill(B3,B6,B7,B8,B9,B13,B15,B16)
HBpop$pop <- "HB"

JBpop <- rbind.fill(B10,B11,B12)
JBpop$pop <- "JB"

LWRpop <- rbind.fill(B4,B5)
LWRpop$pop <- "LWR"

SLpop <- rbind.fill(B14)
SLpop$pop <- "SL"

combined <- rbind.fill(EHApop, HBpop, LWRpop, CSpop, JBpop, SLpop)
meltData <- melt(combined)

# Quick boxplot
boxplot(data=meltData, value~pop)

# Nicer boxplot. This is associated with main figure 2, and supplemental figures S16, S17
beluga_offset <- meltData %>%
  mutate(pop = fct_reorder(pop, value, .fun="mean")) %>%
  ggplot(aes(x=pop, y=value, fill=pop)) + 
  geom_boxplot(width=0.5, lwd=0.5, color="gray30", fill="#a185a9ff", alpha=1, outlier.size = 0.5) +
  theme(legend.position="none")+
  ylab("Genetic offset")+
  theme_bw()+
  theme(legend.position="none", panel.grid.major.x=element_blank())+
  xlab("Population")

beluga_offset
#write.csv(meltData, "beluga_meltdata_for_boxplot_RCP85_2100.csv")
ggsave("beluga_K6_boxplot_RCP45_2100_bypop2.svg", width=3, height=4, dpi=1000)

