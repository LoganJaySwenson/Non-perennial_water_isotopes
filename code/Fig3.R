# Fig 3. Spatial variation in δ18O during summer dry-down period.
library(sf)
library(viridis)
library(lubridate)
library(tidyverse)

# Publication theme
source("code/Theme+Settings.R")

crs <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Read: AIMS isotopes
AIMS_isotopes <- read_csv("data/AIMS_isotopes_coordinates.csv") %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs = crs)

# Read: Stream network
streams <- st_read("data/spatial/stream_network/Konza_stream_network.shp") 
st_crs(streams) <- crs

# Plot!
months <- c("Jun" = "June", "Jul" = "July", "Aug" = "August")
AIMS_isotopes$month <- factor(AIMS_isotopes$month, levels = c( "Jun", "Jul", "Aug"))
AIMS_isotopes$d18OWaterBinned <- cut(AIMS_isotopes$d18OWater, breaks = c(seq(-6.5, -4.5, 0.3), 0.2))
ggplot()+
  geom_sf(data = streams, color = "blue")+
  geom_sf(data = AIMS_isotopes, aes(fill = d18OWaterBinned), color = "black", pch = 21, size = 2.5)+
  scale_fill_manual(values = c("#440154", "#38588c", "#25858e", "#2ab07f", "#86d549", "#fde725"), 
                    labels = c("-6.2 to -5.9", "-5.9 to -5.6", "-5.6 to -5.3", "-5.3 to -5.0", "-5.0 to -4.7", "-4.7 to  0.2"))+
  facet_wrap(~month, labeller = as_labeller(months))+
  labs(x = "", y = "")+
  labs(fill = expression(delta^{18}*"O (‰)"))+
  ggspatial::annotation_scale(location = 'br')+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
        strip.text = element_text(face = 'bold'))
ggsave(path = "figures/", "Fig3.png", dpi=300, width = 190, height = 90, units = "mm")
ggsave(path = "figures/", "Figure3.pdf", dpi=300, width = 190, height = 60, units = "mm")