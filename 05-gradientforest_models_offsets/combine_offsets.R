# Stacking genetic offsets from multiple species. Input files needed are the .tif files saved from "gradientforest_whales_pipe.R". This script is to scale and then add these offset layers.

# Load libraries
library(rnaturalearth) # (Massicotte et al., 2023)
library(rworldmap) # (South, 2023)
library(ggspatial) # (Dunnington et al., 2023)
library(tidyterra) # (Hernangomex et al., 2025)
library(terra) # (Hijmans et al., 2025)
library(svglite) # (Wickham et al., 2025)
library(raster) # (Hijmans, 2025)
library(tidyverse) #(Wickham et al., 2019)

# Load offset rasters
beluga_offset <- raster("raster.offset.RCP85.2100.beluga.easternCanada.RANGECLIP.tif")
narwhal_offset <- raster("raster.offset.RCP85.2100.narwhal.easternCanada.RANGECLIP.tif")
bowhead_offset <- raster("raster.offset.RCP85.2100.bowhead.easternCanada.RANGECLIP.tif")

# Set color scale
offset.colors = colorRampPalette(c("#4575B4", "#abd9e9","#ffffbf","#fdae61","#f46d43","#a50026"))

# Test that rasters loaded in OK
plot(beluga_offset, col=offset.colors(200))
plot(narwhal_offset, col=offset.colors(200))
plot(bowhead_offset, col=offset.colors(200))

# Convert to spatraster
beluga_offset.spat <- rast(beluga_offset)
narwhal_offset.spat <- rast(narwhal_offset)
bowhead_offset.spat <- rast(bowhead_offset)

# Use the scale() function for SDs
beluga_offset_scaled <- scale(beluga_offset.spat)
narwhal_offset_scaled <- scale(narwhal_offset.spat)
bowhead_offset_scaled <- scale(bowhead_offset.spat)

# Look at plots to make sure it worked
plot(beluga_offset_scaled, col=offset.colors(200))
plot(narwhal_offset_scaled, col=offset.colors(200))
plot(bowhead_offset_scaled, col=offset.colors(200))

# Add offset values together with mosaic merge in terra
temp_combined_offset_scaled <- mosaic(beluga_offset_scaled,narwhal_offset_scaled, fun="sum")
combined_offset_scaled <- mosaic(temp_combined_offset_scaled,bowhead_offset_scaled, fun="sum")

# Quick plot
plot(combined_offset_scaled, col=offset.colors(200))

# Prepare stuff for making nicer plot
map_outline <- getMap(resolution="high")
crop.extent <- extent(-110,-45, 44, 80)

# Crop outline to boundary and convert to dataframe
map_outline <- crop(map_outline, y=crop.extent) %>% fortify()

# Final map showing sum of the 3 species. This is associated with main figure 5 (clipped by range), and also associated with supplemental figure S22 (not clipped by range)
svglite("combined_all3sp_GF_GeneticOffset_RCP85_2100_landmap_stepsn_SCALED_RANGECLIP.svg", width=6, height=6)
ggplot()+
  geom_spatraster(data=combined_offset_scaled)+
  geom_polygon(data=map_outline, aes(x=long, y=lat, group=group), fill="gray80", colour="gray50", linewidth=0.2)+
  scale_fill_stepsn(n.breaks = 20, colours = offset.colors(12), na.value="white")+#, limits=grad_lims)+
  scale_x_continuous(expand=c(0,0), breaks=c(-100,-80,-60)) +
  scale_y_continuous(expand=c(0,0), breaks=c(50,60,70,80)) +
  xlab("Longitude")+
  ylab("Latitude")+
  labs(fill="Genetic offset")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))+
  annotation_north_arrow(height = unit(1, "cm"),width = unit(1, "cm"),
                         location = "tr", which_north = "true",
                         pad_x = unit(0.12, "cm"), pad_y = unit(0.3, "cm"),
                         style = ggspatial::north_arrow_fancy_orienteering())+
  theme(text=element_text(size=15))+
  theme(legend.title = element_text( size=12), legend.text=element_text(size=12))
dev.off()
