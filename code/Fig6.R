#Fig 6. Spatial- % of streamwater less than ~3 months in age. 
#streamwater age is defined as the mean of the posterior distribution of source mixtures.
library(lubridate)
library(sf)
library(viridis)
library(tidyverse)

#Publication theme
source("code/Theme+Settings.R")

#Read: AIMS isotopes
AIMS_isotopes <- as_tibble(read.csv("data/AIMS_isotopes_results.csv"))
AIMS_isotopes <- st_as_sf(x = AIMS_isotopes, coords = c("long", "lat"), crs = 4326)

#Read: Stream network
streams <- st_read("data/spatial/stream_network/Konza_stream_network.shp")
st_crs(streams) <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#Plot!
months <- c("Jun" = "June", "Jul" = "July", "Aug" = "August")
AIMS_isotopes$month <- factor(AIMS_isotopes$month, levels = c( "Jun", "Jul", "Aug"))
AIMS_isotopes$s1Binned <- cut(AIMS_isotopes$s1, breaks = c(seq(from = 36, to = 64, by = 4)))
ggplot()+
  geom_sf(data = streams, color = "blue")+
  geom_sf(data = AIMS_isotopes, aes(fill = s1Binned), color = "black", pch = 21, size = 2.5)+
  scale_fill_viridis(direction = -1, discrete = T)+
  facet_wrap(~month, labeller = as_labeller(months))+
  labs(x = "", y = "")+
  labs(fill = "Streamwater\n< 3 months\nin age (%)")+
  ggspatial::annotation_scale(location = 'br')+
  guides(fill = guide_legend(reverse=T))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
        strip.text = element_text(face = 'bold'))
ggsave(path = "figures/", "Fig6.png", dpi=300, width = 190, height = 90, units = "mm")