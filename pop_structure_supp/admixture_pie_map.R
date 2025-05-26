# Making admixture pie maps for all three whale species, using qmatrix files from Muller et al (unpublished) and de Greef et al. (2024); https://doi.org/10.1111/gcb.17528. 
# The plots made from this script correspond to supplemental figures S1, S2, S3.
# Creating admixture pie maps were created using tutorial from Tom Jenkins: https://github.com/Tom-Jenkins/admixture_pie_chart_map_tutorial 

# Load libraries
library(tidyverse) #(Wickham et al., 2019)
library(reshape2) # (Wickham, 2007)
library(raster) # (Hijmans, 2025)
library(rworldmap) # (South, 2023)
library(ggspatial) # (Dunnington et al., 2023)
library(splitstackshape) # (Mahto, 2019)
library(svglite) # (Wickham et al., 2025)


## Load qMatrix file. Run these separately by species. Species-specific parts are noted in comments.

# narwhal: qmatrix_narwhal_K3.csv
# bowhead: qmatrix_bowhead_K2.csv
# beluga: qmatrix_beluga_K6.csv

qmatrix <- read.csv("qmatrix_beluga_K6.csv")

qlong <- melt(qmatrix, id.vars=c("Ind","Site"))
head(qlong)

## Define colour palette - run line 28, 31, or 34 depending on species
# Narwhal colors (matching admixture plots in de Greef et al. 2024)
pal <- colorRampPalette(c("#F64F2E", "#1EB4C4", "#807DBA"))

# Bowhead whale colors (matching admixture plots in de Greef et al. 2024)
pal <- colorRampPalette(c("#FFA500", "#807DBA"))

# Beluga whale colors (matching admixture plots in Muller et al. unpublished)
pal <- colorRampPalette(c("#56B4E9", "#0072B2", "#009E73", "#E69F00", "#F0E442", "#CC79A7"))

cols <- pal(length(unique(qlong$variable)))

# Calculate mean admixture proportions for each site
clusters <- grep("Cluster", names(qmatrix)) # indexes of cluster columns
avg_admix <- aggregate(qmatrix[, clusters], list(qmatrix$Site), mean)

# Order alphabetically by site
avg_admix <- avg_admix[order(as.character(avg_admix$Group.1)), ]
avg_admix

# Convert dataframe from wide to long format
avg_admix <- melt(avg_admix, id.vars = "Group.1")
head(avg_admix)

# Define a function to plot pie charts using ggplot for each site
pie_charts = function(admix_df, site, cols){
  ggplot(data = subset(admix_df, Group.1 == site),
         aes(x = "", y = value, fill = variable))+
    geom_bar(width = 1, stat = "identity", colour = "black", show.legend = FALSE)+
    coord_polar(theta = "y")+
    scale_fill_manual(values = cols)+
    theme_void()
}

# Apply function to all sites using for loop. Run line 63, 66, or 69 depending on species
# narwhal
subsites <- sort(c("AB", "BI", "CR", "GF", "IG", "PB", "PG", "PI", "RB", "RE", "SB"))

# bowhead
subsites <- sort(c("AB", "CD", "CH", "CR", "DIS", "GH", "HB", "IQ", "NF", "PB", "PG", "PI", "RB", "RI", "SB", "WB"))

# beluga
subsites <- sort(c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12", "B13", "B14", "B15", "B16"))

pies <- list()
for (i in subsites){
  pies[[i]] = pie_charts(admix_df = avg_admix, site = i, cols = cols) 
}
pies

## Import csv file containing coordinates. Run these separately by species. Species-specific parts are noted in comments.

# narwhal
sites <- read.csv("sample_locations_all_species.csv")
sites <- subset(sites, Narwhal > 0)
coords <- sites[,c("Location_abbrev", "est_latitude", "est_longitude")]
colnames(coords) <- c("location_ID", "est_latitude", "est_longitude")

# bowhead 
sites <- read.csv("sample_locations_all_species.csv")
sites <- subset(sites, Bowhead > 0)
coords <- sites[,c("Location_abbrev", "est_latitude", "est_longitude")]
colnames(coords) <- c("location_ID", "est_latitude", "est_longitude")
# Add Newfoundland coordinates because this qMatrix includes one individual that was beach-cast and is included in de Greef et al. 2024. 
new <-data.frame("NF", 50.211, -55.487)
names(new)<-c("location_ID", "est_latitude", "est_longitude")
coords <- rbind(coords, new)

# beluga
sites <- read.csv("sample_locations_all_species.csv")
sites <- subset(sites, Beluga=="Y")
coords <- sites[,c("Location_abbrev", "est_latitude", "est_longitude")]
colnames(coords) <- c("location_ID", "est_latitude", "est_longitude")

# Order alphabetically by site
coords <- coords[order(coords$location_ID), ] 
coords

# Check order matches coords order
as.character(avg_admix$Group.1) == as.character(coords$location_ID)

# Set map boundary (xmin, xmax, ymin, ymax)
boundary <- extent(-110,-45, 44, 80)

# Get map outlines from rworldmap package
map.outline <- getMap(resolution = "high")

# Crop to boundary and convert to dataframe
map.outline <- crop(map.outline, y = boundary) %>% fortify()

# Plot basemap
basemap <- ggplot()+
  geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), 
               fill="gray80", colour="gray50", linewidth=0.2)+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_x_continuous(expand=c(0,0), breaks=c(-100,-80,-60)) +
  scale_y_continuous(expand=c(0,0), breaks=c(50,60,70,80)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))+
  # annotation_scale(height=unit(0.15, "cm"), location="tr", aes(width_hint=0.15), text_cex=0.8)+
  annotation_north_arrow(height = unit(1, "cm"),width = unit(1, "cm"),
                         location = "tr", which_north = "true",
                         pad_x = unit(0.12, "cm"), pad_y = unit(0.3, "cm"), #pad_y unit 0.6 if including scale bar
                         style = ggspatial::north_arrow_fancy_orienteering())+
  theme(text=element_text(size=15))+
  theme(legend.title = element_text(size=12), legend.text=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"))
basemap

# Then add pie plots to basemap

# Extract coordinates for each site
coord.list <- list()
for (i in subsites){
  coord.list[[i]] = c(subset(coords, location_ID == i)$est_longitude, subset(coords, location_ID == i)$est_latitude)
}
coord.list

# Define pie chart sizes
radius <- 3

# Convert ggplot pie charts to annotation_custom layers
pies.ac = list()
for (i in 1:length(subsites)){
  pies.ac[[i]] = annotation_custom(grob = ggplotGrob(pies[[i]]),
                                   xmin = coord.list[[i]][[1]] - radius,
                                   xmax = coord.list[[i]][[1]] + radius,
                                   ymin = coord.list[[i]][[2]] - radius,
                                   ymax = coord.list[[i]][[2]] + radius)
}

# Add layers to basemap
pie.map <- basemap + pies.ac
pie.map

# Save figure
ggsave("admix_pie_narwhal.png", width=4.5, height=4.5, dpi=800)

# For beluga whale, save as svg to pull up LWR pie map above the others from overlapping
svglite("admix_pie_beluga.svg", width=4.5, height=4.5)
pie.map
dev.off()