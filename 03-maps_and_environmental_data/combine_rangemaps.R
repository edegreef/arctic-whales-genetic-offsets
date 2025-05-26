# Script for combining range maps for the Arctic whale species (beluga whale, narwhal, bowhead whale) within the eastern Canadian Arctic.

# For installing packages:
#list.of.packages <- c("tidyverse", "raster", "rworldmap", "mappproj", "ggspatial")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

# Load packages
library(tidyverse) #(Wickham et al., 2019)
library(raster) # (Hijmans, 2025)
library(rworldmap) # (South, 2023)
library(mapproj) # (McIlroy et al., 2025)
library(ggspatial) # (Dunnington et al., 2023)

# Load range shape files - these files are not included in the data repository because do not have permission to repost from original source. The range shape files were obtained through the IUCN web browser (beluga whale: https://www.iucnredlist.org/species/6335/50352346; narwhal https://www.iucnredlist.org/species/13704/50367651; bowhead whale https://www.iucnredlist.org/species/2467/50347659). Because each shape file is named "data_0.shp", need to manually create a folder to separate species' files (or rename the .shp file after download).
narwhal_range <- shapefile('~/redlist_species_data_narwhal/data_0.shp')
bowhead_range <- shapefile('~/redlist_species_data_bowhead/data_0.shp')
beluga_range <- shapefile('~/redlist_species_data_beluga/data_0.shp')

# Prepare base map
map_outline <- getMap(resolution="high")
crop.extent <- extent(-110,-45, 44, 80)

# Crop outline to boundary and convert to dataframe
map_outline <- crop(map_outline, y=crop.extent) %>% fortify()

# Load in data for sites points
sites <- read.csv("sample_locations_all_species.csv", header=T)
nar <- subset(sites, Narwhal>0)
bow <- subset(sites, Bowhead>0)
bel <- subset(sites, Beluga=="Y")

# Plot map (this corresponds to Figure 1)
ggplot()+
  geom_polygon(data=beluga_range, aes(x=long, y=lat, group=group), fill="#440d55c4", alpha=0.5)+
  geom_polygon(data=bowhead_range, aes(x=long, y=lat, group=group), fill="#0489eaff", alpha=0.5)+
  geom_polygon(data=narwhal_range, aes(x=long, y=lat, group=group), fill="#fce824c4", alpha=0.5)+
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray80", colour="gray50", linewidth=0.2)+
  geom_point(data=nar, aes(x=est_longitude, y=est_latitude), pch=16, size=2, alpha=1,colour="#f2ff00ff")+
  geom_point(data=bow, aes(x=est_longitude, y=est_latitude), pch=16, size=2, alpha=1,colour="#0489eaff")+
  geom_point(data=bel, aes(x=est_longitude, y=est_latitude), pch=16, size=2, alpha=1,colour="#440d55ff")+
  xlab("Longitude")+
  ylab("Latitude")+
  coord_fixed(ratio=2.1,xlim=c(-110,-45),ylim=c(44,80))+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))+
  annotation_north_arrow(height = unit(1, "cm"),width = unit(1, "cm"),
                         location = "tr", which_north = "true",
                         pad_x = unit(0.12, "cm"), pad_y = unit(0.6, "cm"),
                         style = ggspatial::north_arrow_fancy_orienteering())+
  theme(text=element_text(size=15))+
  theme(legend.title = element_text( size=12), legend.text=element_text(size=12))

# Save as SVG to touch up point aesthetics and legend in inkscape
ggsave("species_range_overlap_shapefiles_crop_updated.svg", width=6, height=7, dpi=800)  
