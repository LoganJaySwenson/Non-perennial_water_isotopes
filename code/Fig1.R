# Fig 1. Map of Konza Prairie
library(sf)
library(stars)
library(raster)
library(whitebox)
library(cowplot)
library(lubridate)
library(tidyverse)

# Publication theme
source("code/Theme+Settings.R")

crs <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Read: Sampling sites
sites <- read_csv("data/AIMS/AIMS_STIC_details.csv") %>%
  mutate(id = "Synoptic Sampling Sites") %>%
  select(id, long, lat) %>%
  add_row(id = "USGS Gage 06879650", long = -96.59469, lat = 39.10207) %>%
  add_row(id = "NEON Tower", long = -96.61294, lat = 39.11045) %>%
  add_row(id = "NEON Sampling Site", long = -96.60274, lat = 39.10449) %>% 
  add_row(id = "Groundwater Wells", long = -96.584221, lat = 39.084785) %>%
  add_row(id = "Meteorological Station", long = -96.608184, lat = 39.100785) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs = crs)

sites_buffer <- st_buffer(x = sites, dist = 200) # set buffer around sites to crop DEM

# Read: Stream network
streams <- st_read("data/spatial/stream_network/Konza_KingsCreek_network.shp") %>%
  st_transform(st_crs(crs)) %>%
  st_crop(sites_buffer) %>% # crop streams to buffer set around sites
  filter(NAME == "Kings Creek")

# Read: Watershed at outlet
watershed <- st_read("data/spatial/outlet/usgs_watershed.shp") %>%
  st_transform(crs = crs)

sites_buffer <- st_union(watershed, sites_buffer)
sites_buffer <- st_bbox(sites_buffer)

# Read: Konza DEM
dem <- raster::raster("data/LTER/DEM/Konza_DEM_Buffer.tif") %>%
  # crop DEM to buffer set around sites
  crop(sites_buffer) %>% 
  # save cropped DEM for use in wbt
  writeRaster("data/spatial/site_output/dem.tif", overwrite = T)

dem_df <- as_tibble(rasterToPoints(dem)) 
colnames(dem_df) = c("lon", "lat", "elevation")

# Plot!
sites$id <- factor(sites$id, levels = c("Synoptic Sampling Sites", "NEON Sampling Site", "NEON Tower", "USGS Gage 06879650", "Groundwater Wells", "Meteorological Station"))
p1 <- 
  ggplot()+
  geom_raster(data = dem_df, aes(x = lon, y = lat, fill = elevation), alpha = 0.7)+
  scale_fill_gradientn(colors = c("#00A600", "#24B300", "#7ACC00", "#ADD900", "#E6E600", "#E8C727","#EAB64E"), 
                       limits = c(318, 444), breaks = seq(320, 440, 20))+
  labs(fill = "Elevation (m)")+
  ggnewscale::new_scale_fill()+
  geom_sf(data = watershed, fill = NA, color = "black", linewidth=0.4)+
  geom_sf(data = streams, color = "blue")+
  geom_sf(data = sites, aes(fill = id, shape = id), size=2.5)+
  scale_shape_manual(values = c(21, 21, 22, 24, 25, 23))+
  scale_fill_manual(values = c("#A84268", "#FCB97D", "#377eb8", "#FF5E5B", "#f781bf", "#984ea3"))+
  labs(fill = NULL , shape = NULL)+
  labs(x = "", y = "")+
  ggspatial::annotation_scale(location = "bl")+
  ggspatial::annotation_north_arrow(location = "tr", height = unit(8, "mm"), width = unit(8, "mm"))+
  scale_x_continuous(expand = c(0,0), labels = function(x) paste0(x, '\u00B0', "W"))+
  scale_y_continuous(expand = c(0,0), labels = function(x) paste0(x, '\u00B0', "N"))+
  theme(strip.text = element_text(face = 'bold'), 
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8), 
        legend.box = "horizontal",
        legend.justification = "bottom")

# Read: States map
states <- st_read("data/spatial/us_state/us_state.shp") %>%
  filter(!STATE %in% c("VI", "PR", "MP", "HI", "AS", "AK", "GU", "PW", "MH", "FM")) %>%
  st_transform(crs = crs)

# Pour point for outlet
pp <- tibble(
  x = 39.09225847,
  y = -96.58724405) %>% 
  st_as_sf(coords = c("y","x"), crs=4326) %>% 
  st_transform(crs = crs)

# Plot!
p2 <- 
  ggplot()+
  geom_sf(data = states, fill = NA, color = "black")+
  geom_sf(data = pp, color = "red")+
  scale_x_continuous(name = "")+
  scale_y_continuous(name = "")+
  coord_sf()+
  theme(panel.border = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = margin(0,0,-15,-15, "pt"))

ggdraw()+
  draw_plot(p1)+
  draw_plot(p2, x = 0.66, y = 0.54, width = 0.32, height = 0.32)

ggsave(path = "figures/", "Fig1.png", dpi=300, width = 190, height = 110, units = "mm")
ggsave(path = "figures/", "Figure1.pdf", dpi=300, width = 190, height = 110, units = "mm")



# Catchment stats & generate drainage area for USGS stream gage for use in Fig. 1
# smooth DEM
wbt_gaussian_filter(
  input = "data/LTER/DEM/Konza_DEM_Buffer.tif",
  output = "data/spatial/site_output/dem_smoothed.tif")

# breach depressions
wbt_breach_depressions(
  dem = "data/spatial/site_output/dem_smoothed.tif",
  output = "data/spatial/site_output/dem_breached.tif",
  fill_pits = F)

# flow direction raster
wbt_d8_pointer(
  dem = "data/spatial/site_output/dem_breached.tif",
  output = "data/spatial/site_output/fdr.tif")

# flow accumulation raster
wbt_d8_flow_accumulation(
  input = "data/spatial/site_output/dem_breached.tif",
  out_type = "cells",
  output = "data/spatial/site_output/fac.tif")

# delineate stream network
wbt_extract_streams(
  flow_accum = "data/spatial/site_output/fac.tif",
  output = "data/spatial/site_output/streams.tif",
  threshold = 22500)

# convert raster streams to vector
wbt_raster_streams_to_vector( 
  streams = "data/spatial/site_output/streams.tif",
  d8_pntr = "data/spatial/site_output/fdr.tif",
  output = "data/spatial/site_output/streams.shp") #convert to polygon

# pour points for Bayesian Unmixing and Young Water Fraction Approaches!
pp <- tibble(
  id = c("Bayesian Unmixing", "Young Water Fraction", "USGS Gage 06879650"),
  x = c(39.093064, 39.10449, 39.10207),
  y = c(-96.588494, -96.60274, -96.59469)) %>%
  st_as_sf(coords = c("y","x"), crs = 4326)%>%
  st_transform(crs = crs) %>%
  st_write("data/spatial/site_output/pourpoints.shp", append = F)

# snap pour points
wbt_jenson_snap_pour_points(
  pour_pts = "data/spatial/site_output/pourpoints.shp",
  streams = "data/spatial/site_output/streams.tif",
  snap_dist = 100,
  output = "data/spatial/site_output/snap.shp")

# delineate watershed
wbt_watershed(
  d8_pntr = "data/spatial/site_output/fdr.tif",
  pour_pts = "data/spatial/site_output/snap.shp",
  output = "data/spatial/site_output/watershed.tif")

# read: watersheds
watersheds <- raster("data/spatial/site_output/watershed.tif")
watersheds <- watersheds %>%
  st_as_stars() %>%
  st_as_sf(merge = T)

# read: stream network
streams <- st_read("data/spatial/stream_network/Konza_KingsCreek_network.shp") %>%
  st_transform(st_crs(crs))

# Plot to see drainage areas
ggplot()+
  geom_sf(data = subset(watersheds, watershed == 1 | watershed == 3))+
  geom_sf(data = streams, color = "blue")+
  geom_sf(data = pp, color = "red")+
  scale_x_continuous(expand = c(0,0), labels = function(x) paste0(x, '\u00B0', "W"))+
  scale_y_continuous(expand = c(0,0), labels = function(x) paste0(x, '\u00B0', "N"))+
  coord_sf()+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))

# Read: Flow accumulation raster
flow_accumulation_raster <- raster::raster("data/spatial/site_output/fac.tif")

# Get drainage area for snapped pour points
pp_drainage_area <- raster::extract(flow_accumulation_raster,
                                        pp,
                                        method = "simple")

pp$DrainageArea_ha <- (pp_drainage_area * 9.36261^2)* (1/1e4) 
pp

# Save: Drainage area for USGS 06879650 
usgs_watershed <- watersheds %>%
  subset(watershed %in% c(1,3)) %>%
  st_union() %>%
  st_write("data/spatial/outlet/usgs_watershed.shp")