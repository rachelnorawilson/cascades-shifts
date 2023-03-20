#################################################################################
### SCRIPT PURPOSE: plot map of NOCA plots with fire history overlain
# Modified from Angert et al. 2018, American Naturalist
# Author: Amy Angert
# last update:  20 Mar 2023

##################################################################################### LOAD LIBRARIES AND PREPARE INPUTS

## Clear workspace
rm(list = ls(all.names = TRUE))

## Libraries needed for spatial stuff 
library(tidyverse)
library(raster)
library(maptools)
library(rgdal)
library(rgeos)
library(sf)
library(maps)
library(mapdata)
library(mapproj)
library(ggstar)
library(geodata)

## Projections
prj.wgs <- "+proj=longlat + type=crs"
prj.lcc <- "+proj=lcc +lon_0=-95 +lat_1=49 +lat_2=77 +type=crs"


## Plot locations
plots <- read_csv("data/0_Lat.Long.csv") %>% 
  drop_na() %>% 
  mutate(Longitude=ifelse(Longitude>0, -Longitude, Longitude)) #input file accidentally has some longitudes in E instead of W

## Plot fire history
fire.points <- read_csv("data/0_All_Plots_Wildfire_Join.csv") %>% 
  mutate(FireHistory = ifelse(CAL_YEAR>=1983, "Burned", "Unburned")) %>% 
  mutate(FireHistory = replace_na(FireHistory, "Unburned"))

## Master plot data frame
fire.plots <- left_join(plots, fire.points, by=c("2015.Plot.Name"="Name"))
burned.plots <- fire.plots %>% filter(FireHistory=="Burned")
unburned.plots <- fire.plots %>% filter(FireHistory=="Unburned")

coordinates(burned.plots) <- ~Longitude+Latitude #convert to spatial data
projection(burned.plots) <- CRS('+proj=longlat') #define projection
burned.plots <- spTransform(burned.plots, CRS=CRS(prj.wgs)) #transform projection 
burned.plots.lcc <- spTransform(burned.plots, CRS=CRS(prj.lcc)) #transform projection 

coordinates(unburned.plots) <- ~Longitude+Latitude #convert to spatial data
projection(unburned.plots) <- CRS('+proj=longlat') #define projection
unburned.plots <- spTransform(unburned.plots, CRS=CRS(prj.wgs))
unburned.plots.lcc <- spTransform(unburned.plots, CRS=CRS(prj.lcc)) #transform projection 

# Define extent of study area
ext <- extent(min(fire.plots$Longitude)-0.1, max(fire.plots$Longitude)+0.1, min(fire.plots$Latitude)-0.1, max(fire.plots$Latitude)+0.1)
bbox <- as(ext, "SpatialPolygons") #convert coordinates to a bounding box

## Park boundary
park <- readOGR("data/shapefiles/park/NOCA_Park_boundary.shp")
park <- spTransform(park, CRS=CRS(prj.wgs))
park.lcc <- spTransform(park, CRS=CRS(prj.lcc))

## Fire polygons
fires <- readOGR("data/shapefiles/fires/NOCA_Wildfire_History.shp")
fires <- spTransform(fires, CRS=CRS(prj.wgs))
fires <- st_as_sf(fires) %>% 
  filter(CAL_YEAR>=1983)
fires.sp <- as(fires, "Spatial")
fires.lcc <- spTransform(fires.sp, CRS=CRS(prj.lcc))

burns <- readOGR("data/shapefiles/fires/Prescribed_burn_history.shp")
burns <- spTransform(burns, CRS=CRS(prj.wgs))
burns <- st_as_sf(burns) %>% 
  filter(CAL_YEAR>=1983)
burns.sp <- as(burns, "Spatial")
burns.lcc <- spTransform(burns.sp, CRS=CRS(prj.lcc))

trtmts <- readOGR("data/shapefiles/fires/Fire_treatment_history.shp")
trtmts <- spTransform(trtmts, CRS=CRS(prj.wgs))
trtmts <- st_as_sf(trtmts) %>% 
  filter(TreatYear>=1983)
trtmts.sp <- as(trtmts, "Spatial")
trtmts.lcc <- spTransform(trtmts.sp, CRS=CRS(prj.lcc))

################################################################################




################################################################################
### Pretty map 

## Set up gridlines & lat/lon labels	
frame.grd <- gridlines(sta.crop)
frame.grd.lcc <- spTransform(frame.grd, CRS=CRS(prj.lcc))
gridatt <- gridat(frame.grd, side="EN")
gridat.lcc <- spTransform(gridatt, CRS=CRS(prj.lcc))

## Set up grayscale color ramp for elevation layer
cuts=seq(0, 3000, by=500) #set breaks
pal <- colorRampPalette(c("white","black"))

## Zoomed out inset (world)
dot.plot <- data.frame(mean.lat=mean(plots$Latitude), mean.long=mean(plots$Longitude))

map_world <- borders("world", colour="black", fill="grey")

ggplot(dot.plot, aes(x = mean.long, y = mean.lat)) +
  map_world +
  geom_star(size=5, pch=21, fill="black") +
  #scale_colour_manual(values=siteColors) +
  coord_map(projection="ortho", orientation=c(48,-121,0)) +
  #scale_shape_manual(values = c(4, 24, 25)) +
  labs(x = "Longitude") +
  labs(y = "Latitude") +
  ggtitle("B") +
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank()) +
  theme_classic()
#legend("topleft", legend="B", bty="n") 
#theme(panel.background=element_rect(fill="#c7eae5"))
ggsave("figures/7_map_world_inset.png", width=8, height=5)

## Zoomed in of plots
#LCC projection
pdf(file="figures/7_map_fire_plots.pdf", width=10, height=8)
plot(park.lcc, border="black") # park boundary
plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) 
plot(burns.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) #prescribed burns layer 
plot(trtmts.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) #prescribed burns layer
plot(unburned.plots.lcc, pch=1, col="black", add=T) #add plots that didn't burn between surveys
plot(burned.plots.lcc, pch=4, col="black", cex=2, add=T) #add plots that burned between surveys
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) #add gridlines
dev.off()

# pdf(file="figures/map_fire_elev.pdf", width=10, height=8)
# plot(park.lcc, border="black") # park boundary
# plot(elev.raster.lcc, breaks=cuts, col=pal(length(cuts)-1), alpha=0.5, add=T)
# plot(park.lcc, border="black", add=T) # park boundary
# plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) 
# plot(burns.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) #prescribed burns layer 
# plot(trtmts.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) #prescribed burns layer
# plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) #add gridlines
# dev.off()

# pdf(file="figures/map_elev_plots.pdf", width=10, height=8)
# plot(park.lcc, border="black") # park boundary
# plot(elev.raster.lcc, breaks=cuts, col=pal(length(cuts)-1), alpha=0.5, add=T)
# plot(park.lcc, border="black", add=T) # park boundary
# plot(unburned.plots.lcc, pch=1, col="black", add=T) #add plots that didn't burn between surveys
# plot(burned.plots.lcc, pch=4, col="black", cex=2, add=T) #add plots that burned between surveys
# plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) #add gridlines
# dev.off()

################################################################################
