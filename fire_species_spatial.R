# Created: April 19, 2023
# Modified: 

#### This script plots fire-experiencing species on a map to visualize their distributions in the study area. Are they spatially correlated? How do distributions differ on western vs eastern slope? It also examines how properties of the fires (unplanned vs prescribed, size, severity) might have affected species dynamics.


### Load libraries
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
library(terra) 
library(stars) 

### Read in data

# List of focal species
load("data/Species.List.Rda")
species.list <- shifts %>% 
  filter(Species.Code!="MOSS") %>% 
  dplyr::select(Species=Species.Code)
species.list$Species <- as.character(species.list$Species)
species.list$Focal <- "YES"

# Plot survey data
dat <- read_csv("data/1_presence_fires_unrarefied.csv") %>% #unfortunately the plot names in this file in column Plot differ between legacy and resurvey time points
  mutate(Plot.Name.2015 = ifelse(Data.Type=="Resurvey", Plot, NA),
         Plot.Name.1980 = ifelse(Data.Type=="Legacy", Plot, NA)) %>% 
  left_join(species.list, by=c("Species.Code"="Species")) %>% 
  filter(Focal=="YES")

dat.leg <- dat %>% filter(Data.Type=="Legacy")
dat.res <- dat %>% filter(Data.Type=="Resurvey")

# Plot locations 
plots <- read_csv("data/0_Lat.Long.csv") %>% 
  drop_na() %>% 
  mutate(Longitude=ifelse(Longitude>0, -Longitude, Longitude)) #input file accidentally has some longitudes in E instead of W
colnames(plots)[2] = "Plot.Name.2015"
colnames(plots)[3] = "Plot.Name.1980"
plots$Plot.Name.1980 = as.character(plots$Plot.Name.1980)

# Join plot location and fire info to resurvey data
dat.leg <- left_join(dat.leg, plots, by=c("Plot"="Plot.Name.1980")) %>% #NOTE: this join is adding 126 rows -- TROUBLESHOOT! 
  dplyr::select(-Plot, -Plot.Name.2015.x)
colnames(dat.leg)[9] = "Plot.Name.2015"
dat.res <- left_join(dat.res, plots, by=c("Plot"="Plot.Name.2015")) %>% #this join is good
  dplyr::select(-Plot, -Plot.Name.1980.x) 
colnames(dat.res)[9] = "Plot.Name.1980"

# Re-merge time points into one frame
dat.all <- bind_rows(dat.leg, dat.res) %>% # this now has plot locations and plot names translated from legacy <--> resurvey
  drop_na() # Latitude=NA for legacy 4042, 4063, resurvey Dia234, Thor221 because these are missing from plots file

dat.pres <- filter(dat.all, Pres.Abs==1)


### Set up maps

# Define projections
prj.wgs <- "+proj=longlat + type=crs"
prj.lcc <- "+proj=lcc +lon_0=-95 +lat_1=49 +lat_2=77 +type=crs"
prj.nad <- crs(b2004)

# Transform to spatial data
coordinates(dat.pres) <- ~Longitude+Latitude #convert to spatial data
projection(dat.pres) <- CRS('+proj=longlat') #define projection
dat.pres <- spTransform(dat.pres, CRS=CRS(prj.wgs)) #transform projection 
dat.pres.lcc <- spTransform(dat.pres, CRS=CRS(prj.lcc)) #transform projection 

# Define extent of study area
ext <- extent(min(dat.pres$Longitude)-0.25, max(dat.pres$Longitude)+0.5, min(dat.pres$Latitude)-0.5, max(dat.pres$Latitude)+0.25)
bbox <- as(ext, "SpatialPolygons") #convert coordinates to a bounding box

# USA polygons (used for lat/lon gridlines) 
sta = readOGR("data/shapefiles/states/gz_2010_us_040_00_500k.shp")
projection(sta) = CRS(prj.wgs)
# Crop state lines to study area
sta.crop <- crop(sta, bbox)

# Park boundary
park <- readOGR("data/shapefiles/park/NOCA_Park_boundary.shp")
park <- spTransform(park, CRS=CRS(prj.wgs))
park.lcc <- spTransform(park, CRS=CRS(prj.lcc))

# Wildfire boundary polygons
fires <- readOGR("data/shapefiles/fires/NOCA_Wildfire_History.shp")
fires <- spTransform(fires, CRS=CRS(prj.wgs))
fires <- st_as_sf(fires) %>% 
  filter(CAL_YEAR>=1983)
fires.sp <- as(fires, "Spatial")
fires.lcc <- spTransform(fires.sp, CRS=CRS(prj.lcc))

# Prescribed burn boundary polygons
burns <- readOGR("data/shapefiles/fires/Prescribed_burn_history.shp")
burns <- spTransform(burns, CRS=CRS(prj.wgs))
burns <- st_as_sf(burns) %>% 
  filter(CAL_YEAR>=1983)
burns.sp <- as(burns, "Spatial")
burns.lcc <- spTransform(burns.sp, CRS=CRS(prj.lcc))

# Fire treatment boundary polygons (these are prescribed burns with various management treatments, like thinning)
trtmts <- readOGR("data/shapefiles/fires/Fire_treatment_history.shp")
trtmts <- spTransform(trtmts, CRS=CRS(prj.wgs))
trtmts <- st_as_sf(trtmts) %>% 
  filter(TreatYear>=1983)
trtmts.sp <- as(trtmts, "Spatial")
trtmts.lcc <- spTransform(trtmts.sp, CRS=CRS(prj.lcc))

# Burn severity maps for years in which fires affected plots
b1990 <- raster("data/shapefiles/MTBS_BSmosaics/mtbs_WA_1990.tif")
b1990.shp <- rasterToPolygons(b1990)
b1990.wgs <- spTransform(b1990.shp, CRS=CRS(prj.wgs))
b1990.crop <- crop(b1990.wgs, ext)
writeOGR(b1990.crop, dsn="data/shapefiles/MTBS_BSmosaics", layer="b1990", driver="ESRI Shapefile")
rm(list=c('b1990', 'b1990.shp', 'b1990.wgs', 'b1990.crop'))

b1993 <- raster("data/shapefiles/MTBS_BSmosaics/mtbs_WA_1993.tif")
b1993.shp <- rasterToPolygons(b1993)
b1993.wgs <- spTransform(b1993.shp, CRS=CRS(prj.wgs))
b1993.crop <- crop(b1993.wgs, ext) #last step does not work
writeOGR(b1993.crop, dsn="data/shapefiles/MTBS_BSmosaics", layer="b1993", driver="ESRI Shapefile")
rm(list=c('b1993', 'b1993.shp', 'b1993.wgs', 'b1993.crop'))

b1994 <- raster("data/shapefiles/MTBS_BSmosaics/mtbs_WA_1994.tif")
b1994.shp <- rasterToPolygons(b1994)
b1994.wgs <- spTransform(b1994.shp, CRS=CRS(prj.wgs))
b1994.crop <- crop(b1994.wgs, ext) 
writeOGR(b1994.crop, dsn="data/shapefiles/MTBS_BSmosaics", layer="b1994", driver="ESRI Shapefile")
rm(list=c('b1994', 'b1994.shp', 'b1994.wgs', 'b1994.crop'))

b1997 <- raster("data/shapefiles/MTBS_BSmosaics/mtbs_WA_1997.tif")
b1997.shp <- rasterToPolygons(b1997)
b1997.wgs <- spTransform(b1997.shp, CRS=CRS(prj.wgs))
b1997.crop <- crop(b1997.wgs, ext) #this step fails
writeOGR(b1997.crop, dsn="data/shapefiles/MTBS_BSmosaics", layer="b1997", driver="ESRI Shapefile")
rm(list=c('b1997', 'b1997.shp', 'b1997.wgs', 'b1997.crop'))

b2002 <- raster("data/shapefiles/MTBS_BSmosaics/mtbs_WA_2002.tif")
b2002.shp <- rasterToPolygons(b2002)
b2002.wgs <- spTransform(b2002.shp, CRS=CRS(prj.wgs))
b2002.crop <- crop(b2002.wgs, ext) 
writeOGR(b2002.crop, dsn="data/shapefiles/MTBS_BSmosaics", layer="b2002", driver="ESRI Shapefile")
rm(list=c('b2002', 'b2002.shp', 'b2002.wgs', 'b2002.crop'))

b2004 <- raster("data/shapefiles/MTBS_BSmosaics/mtbs_WA_2004.tif")
b2004.shp <- rasterToPolygons(b2004)
b2004.wgs <- spTransform(b2004.shp, CRS=CRS(prj.wgs))
b2004.crop <- crop(b2004.wgs, ext) 
writeOGR(b2004.crop, dsn="data/shapefiles/MTBS_BSmosaics", layer="b2004", driver="ESRI Shapefile")
rm(list=c('b2004', 'b2004.shp', 'b2004.wgs', 'b2004.crop'))

b2006 <- raster("data/shapefiles/MTBS_BSmosaics/mtbs_WA_2006.tif")
b2006.shp <- rasterToPolygons(b2006)
b2006.wgs <- spTransform(b2006.shp, CRS=CRS(prj.wgs))
b2006.crop <- crop(b2006.wgs, ext) 
writeOGR(b2006.crop, dsn="data/shapefiles/MTBS_BSmosaics", layer="b2006", driver="ESRI Shapefile")
rm(list=c('b2006', 'b2006.shp', 'b2006.wgs', 'b2006.crop'))

b2009 <- raster("data/shapefiles/MTBS_BSmosaics/mtbs_WA_2009.tif")
b2009.wgs <- projectRaster(b2009, crs=prj.wgs)
b2009.shp <- rasterToPolygons(b2009)
b2009.wgs <- spTransform(b2009.shp, CRS=CRS(prj.wgs))
b2009.crop <- crop(b2009.wgs, ext) 
rm(list=c('b2009', 'b2009.shp', 'b2009.wgs'))

b2010 <- raster("data/shapefiles/MTBS_BSmosaics/mtbs_WA_2010.tif")
b2010.shp <- rasterToPolygons(b2010)
b2010.wgs <- spTransform(b2010.shp, CRS=CRS(prj.wgs))
b2010.crop <- crop(b2010.wgs, ext) 
rm(list=c('b2010', 'b2010.shp', 'b2010.wgs'))

b2014 <- raster("data/shapefiles/MTBS_BSmosaics/mtbs_WA_2014.tif")
b2014.shp <- rasterToPolygons(b2014)
b2014.wgs <- spTransform(b2014.shp, CRS=CRS(prj.wgs))
b2014.crop <- crop(b2014.wgs, ext) 
rm(list=c('b2014', 'b2014.shp', 'b2014.wgs'))

#stack
severity <- stack(b1990, b1993)

  
# Ecoregions polygons
ecoreg <- readOGR("data/shapefiles/Ecoregions_of_the_Pacific_Northwest/Ecoregions_of_the_Pacific_Northwest.shp")
ecoreg <- spTransform(ecoreg, CRS=CRS(prj.wgs))
ecoreg.crop <- crop(ecoreg, bbox)
ecoreg.crop.lcc <- spTransform(ecoreg.crop, CRS=CRS(prj.lcc))

# Set up gridlines & lat/lon labels	
frame.grd <- gridlines(ecoreg.crop)
frame.grd.lcc <- spTransform(frame.grd, CRS=CRS(prj.lcc))
gridatt <- gridat(frame.grd, side="EN")
gridat.lcc <- spTransform(gridatt, CRS=CRS(prj.lcc))


### Maps of fire-experiencing species

# ACMI
pdf(file="figures/map_ACMI_plots.pdf", width=10, height=8)
plot(park.lcc, border="black") # park boundary
plot(ecoreg.crop.lcc, col=rgb(.5,0.8,1,0.3), border="blue", add=T) # ecoregions polygons
plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) # wildfire polygons
plot(burns.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # prescribed burns polygons 
plot(trtmts.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # fire treatments polygons
plot(subset(dat.pres.lcc, Species.Code=="ACMI" & Data.Type=="Legacy"), pch=1, cex=1, col="black", add=T) # historical presences of ACMI
plot(subset(dat.pres.lcc, Species.Code=="ACMI"& Data.Type=="Resurvey"), pch=4, cex=1, col="black", add=T) # contemporary presences of ACMI
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) # gridlines
dev.off()

# ARUV
pdf(file="figures/map_ARUV_plots.pdf", width=10, height=8)
plot(park.lcc, border="black") # park boundary
plot(ecoreg.crop.lcc, col=rgb(.5,0.8,1,0.3), border="blue", add=T) # ecoregions polygons
plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) # wildfire polygons
plot(burns.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # prescribed burns polygons 
plot(trtmts.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # fire treatments polygons
plot(subset(dat.pres.lcc, Species.Code=="ARUV" & Data.Type=="Legacy"), pch=1, cex=1, col="black", add=T) # historical presences of ARUV
plot(subset(dat.pres.lcc, Species.Code=="ARUV"& Data.Type=="Resurvey"), pch=4, cex=1, col="black", add=T) # contemporary presences of ARUV
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) # gridlines
dev.off()

# VAME
pdf(file="figures/map_VAME_plots.pdf", width=10, height=8)
plot(park.lcc, border="black") # park boundary
plot(ecoreg.crop.lcc, col=rgb(.5,0.8,1,0.3), border="blue", add=T) # ecoregions polygons
plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) # wildfire polygons
plot(burns.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # prescribed burns polygons 
plot(trtmts.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # fire treatments polygons
plot(subset(dat.pres.lcc, Species.Code=="VAME" & Data.Type=="Legacy"), pch=1, cex=1, col="black", add=T) # historical presences of VAME
plot(subset(dat.pres.lcc, Species.Code=="VAME"& Data.Type=="Resurvey"), pch=4, cex=1, col="black", add=T) # contemporary presences of VAME
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) # gridlines
dev.off()

# CARU
pdf(file="figures/map_CARU_plots.pdf", width=10, height=8)
plot(park.lcc, border="black") # park boundary
plot(ecoreg.crop.lcc, col=rgb(.5,0.8,1,0.3), border="blue", add=T) # ecoregions polygons
plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) # wildfire polygons
plot(burns.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # prescribed burns polygons 
plot(trtmts.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # fire treatments polygons
plot(subset(dat.pres.lcc, Species.Code=="CARU" & Data.Type=="Legacy"), pch=1, cex=1, col="black", add=T) # historical presences of CARU
plot(subset(dat.pres.lcc, Species.Code=="CARU"& Data.Type=="Resurvey"), pch=4, cex=1, col="black", add=T) # contemporary presences of CARU
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) # gridlines
dev.off()

# CEVE
pdf(file="figures/map_CEVE_plots.pdf", width=10, height=8)
plot(park.lcc, border="black") # park boundary
plot(ecoreg.crop.lcc, col=rgb(.5,0.8,1,0.3), border="blue", add=T) # ecoregions polygons
plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) # wildfire polygons
plot(burns.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # prescribed burns polygons 
plot(trtmts.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # fire treatments polygons
plot(subset(dat.pres.lcc, Species.Code=="CEVE" & Data.Type=="Legacy"), pch=1, cex=1, col="black", add=T) # historical presences of CEVE
plot(subset(dat.pres.lcc, Species.Code=="CEVE"& Data.Type=="Resurvey"), pch=4, cex=1, col="black", add=T) # contemporary presences of CEVE
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) # gridlines
dev.off()

# EPAN (CHAN)
pdf(file="figures/map_EPAN_plots.pdf", width=10, height=8)
plot(park.lcc, border="black") # park boundary
plot(ecoreg.crop.lcc, col=rgb(.5,0.8,1,0.3), border="blue", add=T) # ecoregions polygons
plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) # wildfire polygons
plot(burns.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # prescribed burns polygons 
plot(trtmts.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # fire treatments polygons
plot(subset(dat.pres.lcc, Species.Code=="EPAN" & Data.Type=="Legacy"), pch=1, cex=1, col="black", add=T) # historical presences of EPAN
plot(subset(dat.pres.lcc, Species.Code=="EPAN"& Data.Type=="Resurvey"), pch=4, cex=1, col="black", add=T) # contemporary presences of EPAN
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) # gridlines
dev.off()

# PAMY
pdf(file="figures/map_PAMY_plots.pdf", width=10, height=8)
plot(park.lcc, border="black") # park boundary
plot(ecoreg.crop.lcc, col=rgb(.5,0.8,1,0.3), border="blue", add=T) # ecoregions polygons
plot(fires.lcc, col=rgb(1,0,0,0.3), border="red4", add=T) # wildfire polygons
plot(burns.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # prescribed burns polygons 
plot(trtmts.lcc, col=rgb(1,0.7,0,0.3), border="orange", add=T) # fire treatments polygons
plot(subset(dat.pres.lcc, Species.Code=="PAMY" & Data.Type=="Legacy"), pch=1, cex=1, col="black", add=T) # historical presences of PAMY
plot(subset(dat.pres.lcc, Species.Code=="PAMY"& Data.Type=="Resurvey"), pch=4, cex=1, col="black", add=T) # contemporary presences of PAMY
plot(frame.grd.lcc, add=TRUE, lty="dashed", col="grey", lwd=1) # gridlines
dev.off()


### Who had legacy records on the west vs east side (i.e., what was opportunity to get to burned areas?)

# Spatial vector of polygons created from ecoregions shapefile
ecoreg.vect <- vect(ecoreg) 

# Subset to ecoregions of interest
regions <- ecoreg.vect[ecoreg.vect$NAME == 'North Cascades' | 
                         ecoreg.vect$NAME == 'East Cascades', ]
# Check projection
prj.reg <- crs(regions)                     

# Transform plot presences to same projection
dat.pres.reg <- spTransform(dat.pres, CRS=CRS(prj.reg)) #transform projection 

# Transform to sf objects
pres_pt <- st_as_sf(x = dat.pres.reg)#, 
ecoregion_poly <- st_as_sf(regions)

# Check that points map onto polygon
plot(ecoregion_poly)
points(pres_pt, cex=2) # not showing up

# New object with points that are within ecoregion polygon
pres_ecoregion <- st_intersection(ecoregion_poly, pres_pt) 

# Convert back to data frame
# Extract columns of interest
pres_reg_df <- pres_ecoregion[, c('NAME', 'geometry', 'Species.Code', 'Pres.Abs', 'Fires', 'Elevation.m', 'Data.Type', 'Plot.Name.1980', 'Plot.Name.2015')] 
# Extrat lat-long
pres_reg_df$lon.lat <- substr(as.character(pres_reg_df$geometry), 3, (nchar(as.character(pres_reg_df$geometry))-1)) #remove parentheses

foo <- separate(pres_reg_df, lon.lat, into = c('Longitude', 'Latitude'), sep = ', ') #separate by lon, lat
pres_reg_df <- as.data.frame(foo) #remove spatial attributes

pres_reg_df$Longitude <- as.numeric(pres_reg_df$Longitude) #turn lon-lat to numeric
pres_reg_df$Latitude <- as.numeric(pres_reg_df$Latitude)

table(pres_reg_df$NAME, 
      pres_reg_df$Species.Code, 
      pres_reg_df$Data.Type)


### What is distribution of distances between a species' historical presences and burned plots? (i.e., what was opportunity to colonize burned areas?)

# Pick out burned plots (n=38)
dat.burn <- dat.all %>% 
  filter(Fires=="Burned") %>% 
  dplyr::select(Latitude, Longitude) %>% 
  unique() 

# Transform to spatial data
coordinates(dat.burn) <- ~Longitude+Latitude #convert to spatial data
projection(dat.burn) <- CRS('+proj=longlat') #define projection
dat.burn <- spTransform(dat.burn, CRS=CRS(prj.wgs)) #transform projection 

# All presence records found during legacy survey (n=2109)
dat.pres.leg <- dat.all %>% 
  filter(Data.Type=="Legacy",
         Pres.Abs==1)

# Transform to spatial data
coordinates(dat.pres.leg) <- ~Longitude+Latitude #convert to spatial data
projection(dat.pres.leg) <- CRS('+proj=longlat') #define projection
dat.pres.leg <- spTransform(dat.pres.leg, CRS=CRS(prj.wgs)) #transform projection 

# Calculate Great Circle distance between each burned plot and every resurvey presence
dist.2.burn <- spDists(dat.pres.leg, dat.burn, longlat=T) # returns matrix with 38 columns and 2109 rows, so each column is the distance to a single burned plot and each row is one historical presence for a species

# Join back to presence data to identify species and plots
dat.burn.distances <- bind_cols(as.data.frame(dat.pres.leg), as.data.frame(dist.2.burn))

# Make taller dataframe
dat.burn.distances.tall <- pivot_longer(dat.burn.distances, cols=starts_with("V"), values_to="BurnPlotDist")

ggplot(data=dat.burn.distances.tall, aes(x=BurnPlotDist)) +
  geom_histogram() +
  facet_wrap(~Species.Code) +
  xlab("Distance from legacy presences to burned plots (km)") 


### What kinds of burns did burned plots experience? Categorize as natural, prescribed, etc.

# Spatial vector of polygons created from different fire maps
fires.vect <- vect(fires)
prescr.vect <- vect(burns) 
trtmts.vec <- vect(trtmts)

# Check projections - all the same
prj.fires <- crs(fires.vect)                     
prj.prescr <- crs(prescr.vect)                     
prj.trtmts <- crs(trtmts.vec)                     

# Transform plot presences to same projection
dat.pres.prescr <- spTransform(dat.pres, CRS=CRS(prj.prescr)) #transform projection 

# Transform to sf objects
pres_pt <- st_as_sf(x = dat.pres.prescr) 
fires_poly <- st_as_sf(fires.vect)
prescr_poly <- st_as_sf(prescr.vect)
trtmts_poly <- st_as_sf(trtmts.vec)

# New objects with points that are within different fire polygons
sf_use_s2(FALSE)
pres_fires <- st_intersection(fires_poly, pres_pt) 
pres_prescr <- st_intersection(prescr_poly, pres_pt) 
pres_trtmts <- st_intersection(trtmts_poly, pres_pt) 

# Convert back to data frames
pres_fires_df <- pres_fires[, c('CAUSE', 'CAL_YEAR', 'Acres', 'Species.Code', 'Pres.Abs', 'Fires', 'Elevation.m', 'Data.Type', 'Plot.Name.1980', 'Plot.Name.2015')] #extract columns of interest
pres_fires_df$lon.lat <- substr(as.character(pres_fires_df$geometry), 3, (nchar(as.character(pres_fires_df$geometry))-1)) #extract lat-long, remove parentheses
foo <- separate(pres_fires_df, lon.lat, into = c('Longitude', 'Latitude'), sep = ', ') #separate by lon, lat
pres_fires_df <- as.data.frame(foo) #remove spatial attributes
pres_fires_df$Longitude <- as.numeric(pres_fires_df$Longitude) #turn lon-lat to numeric
pres_fires_df$Latitude <- as.numeric(pres_fires_df$Latitude)
colnames(pres_fires_df)[1] <- "TYPE"
colnames(pres_fires_df)[3] <- "ACRES"
unique(pres_fires_df$ACRES)
unique(pres_fires_df$Plot.Name.2015)
table(pres_fires_df$TYPE,
      pres_fires_df$Plot.Name.2015,
      pres_fires_df$CAL_YEAR)
table(pres_fires_df$Plot.Name.2015,
      pres_fires_df$TYPE)
table(pres_fires_df$TYPE,
      pres_fires_df$CAL_YEAR)

pres_prescr_df <- pres_prescr[, c('NAME', 'CAL_YEAR', 'ACRES', 'Species.Code', 'Pres.Abs', 'Fires', 'Elevation.m', 'Data.Type', 'Plot.Name.1980', 'Plot.Name.2015')] #extract columns of interest
pres_prescr_df$lon.lat <- substr(as.character(pres_prescr_df$geometry), 3, (nchar(as.character(pres_prescr_df$geometry))-1)) #extract lat-long, remove parentheses
foo <- separate(pres_prescr_df, lon.lat, into = c('Longitude', 'Latitude'), sep = ', ') #separate by lon, lat
pres_prescr_df <- as.data.frame(foo) #remove spatial attributes
pres_prescr_df$Longitude <- as.numeric(pres_prescr_df$Longitude) #turn lon-lat to numeric
pres_prescr_df$Latitude <- as.numeric(pres_prescr_df$Latitude)
colnames(pres_prescr_df)[1] <- "TYPE"
unique(pres_prescr_df$ACRES)
unique(pres_prescr_df$Plot.Name.2015)
table(pres_prescr_df$Plot.Name.2015,
      pres_prescr_df$CAL_YEAR)
table(pres_prescr_df$Plot.Name.2015,
      pres_prescr_df$TYPE)

pres_trtmts_df <- pres_trtmts[, c('Treatment', 'TreatYear', 'acres', 'Species.Code', 'Pres.Abs', 'Fires', 'Elevation.m', 'Data.Type', 'Plot.Name.1980', 'Plot.Name.2015')] #extract columns of interest
pres_trtmts_df$lon.lat <- substr(as.character(pres_trtmts_df$geometry), 3, (nchar(as.character(pres_trtmts_df$geometry))-1)) #extract lat-long, remove parentheses
foo <- separate(pres_trtmts_df, lon.lat, into = c('Longitude', 'Latitude'), sep = ', ') #separate by lon, lat
pres_trtmts_df <- as.data.frame(foo) #remove spatial attributes
pres_trtmts_df$Longitude <- as.numeric(pres_trtmts_df$Longitude) #turn lon-lat to numeric
pres_trtmts_df$Latitude <- as.numeric(pres_trtmts_df$Latitude)
colnames(pres_trtmts_df)[1] <- "TYPE"
colnames(pres_trtmts_df)[2] <- "CAL_YEAR"
colnames(pres_trtmts_df)[3] <- "ACRES"
unique(pres_trtmts_df$ACRES)
unique(pres_trtmts_df$Plot.Name.2015)
table(pres_trtmts_df$Plot.Name.2015,
      pres_trtmts_df$TYPE)
table(pres_trtmts_df$Plot.Name.2015,
      pres_trtmts_df$TYPE, 
      pres_trtmts_df$CAL_YEAR)

fire.types <- rbind(pres_fires_df, pres_prescr_df, pres_trtmts_df)

table(fire.types$Plot.Name.2015, 
      fire.types$TYPE)
table(fire.types$Plot.Name.2015[fire.types$Species.Code=="ACMI"], 
      fire.types$TYPE[fire.types$Species.Code=="ACMI"],
      fire.types$CAL_YEAR[fire.types$Species.Code=="ACMI"])
table(fire.types$Plot.Name.2015[fire.types$Species.Code=="ARUV"], 
      fire.types$TYPE[fire.types$Species.Code=="ARUV"],
      fire.types$CAL_YEAR[fire.types$Species.Code=="ARUV"])
table(fire.types$Plot.Name.2015[fire.types$Species.Code=="VAME"], 
      fire.types$TYPE[fire.types$Species.Code=="VAME"],
      fire.types$CAL_YEAR[fire.types$Species.Code=="VAME"])
table(fire.types$Plot.Name.2015[fire.types$Species.Code=="CARU"], 
      fire.types$TYPE[fire.types$Species.Code=="CARU"]
      fire.types$CAL_YEAR[fire.types$Species.Code=="CARU"])
table(fire.types$Plot.Name.2015[fire.types$Species.Code=="CEVE"], 
      fire.types$TYPE[fire.types$Species.Code=="CEVE"],
      fire.types$CAL_YEAR[fire.types$Species.Code=="CEVE"])
table(fire.types$Plot.Name.2015[fire.types$Species.Code=="EPAN"], 
      fire.types$TYPE[fire.types$Species.Code=="EPAN"],
      fire.types$CAL_YEAR[fire.types$Species.Code=="EPAN"])
table(fire.types$Plot.Name.2015[fire.types$Species.Code=="PAMY"], 
      #fire.types$TYPE[fire.types$Species.Code=="PAMY"],
      fire.types$CAL_YEAR[fire.types$Species.Code=="PAMY"])





