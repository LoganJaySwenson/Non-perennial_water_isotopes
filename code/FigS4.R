#Fig S4. Spatial variation in δ18O during the summer dry-down period with limestone units.
library(lubridate)
library(sf)
library(raster)
library(viridis)
library(tidyverse)

#Publication theme
source("code/Theme+Settings.R")

#Read: AIMS isotopes
AIMS_isotopes <- as_tibble(read.csv("data/AIMS_isotopes_coordinates.csv"))
AIMS_isotopes <- st_as_sf(x = AIMS_isotopes, coords = c("long", "lat"), crs = 4326)

#Read: Stream network
streams <- st_read("data/spatial/stream_network/Konza_stream_network.shp")
st_crs(streams) <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
streams_buffer <- st_buffer(x = streams, dist = 200)

#Top of Cottonwood Limestone is 04M03
sites <- as_tibble(read.csv("data/AIMS/AIMS_STIC_details.csv"))
ref <- round(sites %>% filter(siteID == "04M03") %>% pull(Elevation_m), digits = 2)

#Limestone/shale members
limestones <- as_tibble(read.csv("data/LTER/Konza_GW_Well_Info.csv")) %>%
  mutate(order = 1:21) %>%
  arrange(-order) %>%
  mutate(cumulative = cumsum(thickness_m),
         top = cumulative + ref - 17.0,
         bottom = top - thickness_m) %>%
  filter(type == "Ls") %>%
  print(., n=11)

#Read: Konza DEM
dem.raster <- raster::raster("data/LTER/DEM/Konza_DEM_Buffer.tif")
dem.raster <- crop(x = dem.raster, y = streams_buffer) #crop DEM to buffer set around sites
dem.m  <-  rasterToPoints(dem.raster)
dem.df <-  data.frame(dem.m)
colnames(dem.df) = c("lon", "lat", "elevation")
dem.df <- dem.df %>%
  mutate(elevation.binned = case_when(
    elevation <= 429.02 & elevation >= 423.47 ~ "Florence Ls",
    elevation <= 416.52 & elevation >= 412.57 ~ "Kinney Ls",
    elevation <= 412.57 & elevation >= 405.80 ~ "Wymore Ls",
    elevation <= 400.92 & elevation >= 398.44 ~ "Threemile Ls",
    elevation <= 394.08 & elevation >= 391.23 ~ "Funston Ls",
    elevation <= 385.12 & elevation >= 382.68 ~ "Crouse Ls",
    elevation <= 378.32 & elevation >= 377.34 ~ "Middleburg Ls",
    elevation <= 375.16 & elevation >= 373.22 ~ "Eiss Ls",
    elevation <= 369.07 & elevation >= 368.31 ~ "Morrill Ls",
    elevation <= 364.96 & elevation >= 363.21 ~ "Cottonwood Ls",
    elevation <= 355.06 & elevation >= 350.82 ~ "Neva Ls",)) %>%
  drop_na(elevation.binned)

#Plot!
months <- c("Jun" = "June", "Jul" = "July", "Aug" = "August")
AIMS_isotopes$month <- factor(AIMS_isotopes$month, levels = c( "Jun", "Jul", "Aug"))
AIMS_isotopes$d18OWaterBinned <- cut(AIMS_isotopes$d18OWater, breaks = c(seq(-6.5, -4.5, 0.25), 0.2))
ggplot(dem.df)+
  #geom_raster(aes(x = lon, y = lat, fill = elevation.binned))+
  geom_raster(aes(x = lon, y = lat), fill = "#FEF3A3")+
  guides(fill = "none")+
  ggnewscale::new_scale_fill()+
  geom_sf(data = streams, color = "blue")+
  geom_sf(data = AIMS_isotopes, aes(fill = d18OWaterBinned), color = "black", pch = 21, size = 2.5)+
  scale_fill_viridis(discrete = T)+
  facet_wrap(~month, labeller = as_labeller(months))+
  labs(x = "", y = "")+
  labs(fill = expression(delta^{18}*"O (‰)"))+
  ggspatial::annotation_scale(location = 'br')+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
        strip.text = element_text(face = 'bold'))
ggsave(path = "figures/", "FigS4.png", dpi=300, width = 190, height = 90, units = "mm")