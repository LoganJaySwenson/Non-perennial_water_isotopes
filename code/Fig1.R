#Fig 1. Map
library(lubridate)
library(sf)
library(stars)
library(raster)
library(whitebox)
library(viridis)
library(patchwork)
library(tidyverse)

#Publication theme
source("code/Theme+Settings.R")

#Read: Sites
sites <- read.csv("data/AIMS/AIMS_STIC_details.csv")
sites <- as_tibble(sites) %>%
  mutate(id = "Sampling Sites") %>%
  select(id, long, lat) %>%
  add_row(id = "USGS Gage 06879650", long = -96.59469, lat = 39.10207) %>%
  add_row(id = "NEON Tower", long = -96.61294, lat = 39.11045) %>%
  add_row(id = "NEON Sampling Site", long = -96.60274, lat = 39.10449) %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs = st_crs("+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
sites.buffer <- st_buffer(x = sites, dist = 200) #set buffer around sites to crop DEM
#sites.buffer <- st_buffer(x = sites, dist = 200) #set buffer around streams to crop DEM

#Read: Konza DEM
dem.raster <- raster::raster("data/LTER/DEM/Konza_DEM_Buffer.tif")
dem.raster <- crop(x = dem.raster, y = sites.buffer) #crop DEM to buffer set around sites
writeRaster(dem.raster, "data/spatial/site_output/dem.tif", overwrite = T) #save cropped DEM for use in wbt
dem.m  <-  rasterToPoints(dem.raster)
dem.df <-  data.frame(dem.m)
colnames(dem.df) = c("lon", "lat", "alt")

#Read: Stream network
streams <- st_read("data/spatial/stream_network/Konza_KingsCreek_network.shp")
st_crs(streams) <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
streams <- st_crop(x = streams, y = sites.buffer) #crop streams to buffer set around sites

#Plot!
sites$id <- factor(sites$id, levels = c("Sampling Sites","USGS Gage 06879650", "NEON Sampling Site", "NEON Tower"))
p1 <- 
  ggplot()+
  geom_raster(data = dem.df, aes(x = lon, y = lat, fill = alt), alpha = 0.5)+
  scale_fill_gradientn(colors = c("#00A600", "#24B300", "#7ACC00", "#ADD900", "#E6E600", "#E8C727","#EAB64E"), 
                       limits = c(318, 440), breaks = seq(320, 440, 20))+
  labs(fill = "Elevation (m)")+
  ggnewscale::new_scale_fill()+
  geom_sf(data = streams, color = "blue")+
  geom_sf(data = sites, aes(fill = id, shape = id), size=2)+
  scale_shape_manual(values = c(21, 24, 22, 23))+
  scale_fill_manual(values = c("#636363", "#636363", "#636363", "#636363"))+
  labs(fill = "", shape = "")+
  labs(x = "", y = "")+
  ggspatial::annotation_scale(location = 'br')+
  scale_x_continuous(expand = c(0,0), labels = function(x) paste0(x, '\u00B0', "W"))+
  scale_y_continuous(expand = c(0,0), labels = function(x) paste0(x, '\u00B0', "N"))+
  theme(strip.text = element_text(face = 'bold'), 
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))

guides(color = guide_legend(order = 0),
       fill  = guide_legend(order = 1))

#States map
states <- map_data("state")

#Pour point for outlet
pp <- tibble(
  x = 39.09225847,
  y = -96.58724405) %>% 
  st_as_sf(., coords = c("y","x"), crs=4326) %>% 
  st_transform(., crs = st_crs("+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

#Pour point to represent site location
pp <- st_transform(pp, crs = 4326)

#Plot!
p2 <- 
  ggplot()+
  geom_polygon(data = states, aes(x = long, y = lat, group = group), fill = NA, color = "black")+
  geom_sf(data = pp, color = "red")+
  scale_x_continuous(name = "")+
  scale_y_continuous(name = "")+
  coord_sf()+
  theme(panel.border = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), axis.ticks.y=element_blank())

p2 + p1 + plot_layout(ncol = 2, widths = c(1,2)) + plot_annotation(tag_levels = "a")
ggsave(path = "figures/", "Fig1.png", dpi=300, width = 190, height = 110, units = "mm")



#Catchment stats!
#smooth DEM
wbt_gaussian_filter(
  input = "data/LTER/DEM/Konza_DEM_Buffer.tif",
  output = "data/spatial/site_output/dem_smoothed.tif")

#breach depressions
wbt_breach_depressions(
  dem = "data/spatial/site_output/dem_smoothed.tif",
  output = "data/spatial/site_output/dem_breached.tif",
  fill_pits = F)

#flow direction raster
wbt_d8_pointer(
  dem = "data/spatial/site_output/dem_breached.tif",
  output = "data/spatial/site_output/fdr.tif")

#flow accumulation raster
wbt_d8_flow_accumulation(
  input = "data/spatial/site_output/dem_breached.tif",
  out_type = "cells",
  output = "data/spatial/site_output/fac.tif")

#delineate stream network
wbt_extract_streams(
  flow_accum = "data/spatial/site_output/fac.tif",
  output = "data/spatial/site_output/streams.tif",
  threshold = 22500)

#convert raster streams to vector
wbt_raster_streams_to_vector( 
  streams = "data/spatial/site_output/streams.tif",
  d8_pntr = "data/spatial/site_output/fdr.tif",
  output = "data/spatial/site_output/streams.shp") #convert to polygon

#pour points for Bayesian Unmixing and Young Water Fraction Approaches!
pp <- tibble(
  id = c("Bayesian Unmixing", "Young Water Fraction"),
  x = c(39.093064, 39.10449),
  y = c(-96.588494, -96.60274)) %>%
  st_as_sf(., coords = c("y","x"), crs = 4326)%>%
  st_transform(., crs = st_crs("+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) %>%
  st_write("data/spatial/site_output/pourpoints.shp", append = F)

#snap pour points
wbt_jenson_snap_pour_points(
  pour_pts = "data/spatial/site_output/pourpoints.shp",
  streams = "data/spatial/site_output/streams.tif",
  snap_dist = 100,
  output = "data/spatial/site_output/snap.shp")

#delineate watershed
wbt_watershed(
  d8_pntr = "data/spatial/site_output/fdr.tif",
  pour_pts = "data/spatial/site_output/snap.shp",
  output = "data/spatial/site_output/watershed.tif")

#read: watersheds
watersheds <- raster("data/spatial/site_output/watershed.tif")
watersheds <- watersheds %>%
  st_as_stars() %>%
  st_as_sf(., merge = TRUE)

#read: stream network
streams <- st_read("data/spatial/stream_network/Konza_KingsCreek_network.shp")
st_crs(streams) <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#plot!
ggplot()+
  geom_sf(data = watersheds)+
  geom_sf(data = streams, color = "blue")+
  geom_sf(data = pp, color = "red")+
  scale_x_continuous(expand = c(0,0), labels = function(x) paste0(x, '\u00B0', "W"))+
  scale_y_continuous(expand = c(0,0), labels = function(x) paste0(x, '\u00B0', "N"))+
  coord_sf()+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))

#read: flow accumulation raster
flow.accumulation.raster <- raster("data/spatial/site_output/fac.tif")

#calculate: contributing area for snapped pour points!
pp.contributing.area <- raster::extract(x = flow.accumulation.raster,
                                        y = pp,
                                        method = "simple")
pp$ContributingArea_ha <- pp.contributing.area * 1/10000
pp