#Fig 2. Spatial variation in δ18O during the summer dry-down period.
library(lubridate)
library(sf)
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

#Plot!
months <- c("Jun" = "June", "Jul" = "July", "Aug" = "August")
AIMS_isotopes$month <- factor(AIMS_isotopes$month, levels = c( "Jun", "Jul", "Aug"))
AIMS_isotopes$d18OWaterBinned <- cut(AIMS_isotopes$d18OWater, breaks = c(seq(-6.5, -4.5, 0.25), 0.2))
ggplot()+
  geom_sf(data = streams, color = "blue")+
  geom_sf(data = AIMS_isotopes, aes(fill = d18OWaterBinned), color = "black", pch = 21, size = 2.5)+
  scale_fill_viridis(discrete = T)+
  facet_wrap(~month, labeller = as_labeller(months))+
  labs(x = "", y = "")+
  labs(fill = expression(delta^{18}*"O (‰)"))+
  ggspatial::annotation_scale(location = 'br')+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
        strip.text = element_text(face = 'bold'))
ggsave(path = "figures/", "Fig3.png", dpi=300, width = 190, height = 90, units = "mm")