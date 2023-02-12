#Fig 4. Variation in δ18O with distance to outlet + GW Levels!
library(lubridate)
library(sf)
library(raster)
library(viridis)
library(riverdist)
library(tidyverse)

#Publication theme
source("code/Theme+Settings.R")

#Read: AIMS isotopes
AIMS_isotopes <- as_tibble(read.csv("data/AIMS_isotopes_RF.csv"))
AIMS_isotopes <- AIMS_isotopes %>%
  select(1:18) %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4326)

#Read: Sample pts
pts <- as_tibble(read.csv("data/AIMS/AIMS_STIC_details.csv"))
pts <- pts %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4326) %>%
  st_transform(., crs = "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") #reproject vector to planar coordinates

#Read: Stream network
streams <- st_read("data/spatial/stream_network/Konza_stream_network.shp")
st_crs(streams) <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
streams <- streams %>%
  st_transform(., crs = "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") #reproject vector to planar coordinates

#Prep flownet
streams <- streams %>% 
  #remove z information (to make it a 'linear' feature)
  st_zm() %>% 
  #Add ID data
  mutate(
    uid = seq(1, nrow(.)),
    type = 'river')

#export streams to spatial file
st_write(streams, dsn = "data/spatial/distance_to_outlet/streams.shp", layer = "streams", append = F)

#Create flownet
flow_net <- line2network(path = "data/spatial/distance_to_outlet/streams.shp", layer="streams", tolerance = 1)

#save flow network 
save(flow_net, file = "data/spatial/distance_to_outlet/streams.rda")

#Prep sample pts
pts <- pts %>% 
  #Define coordinates
  mutate(
    x = st_coordinates(.)[,1],
    y = st_coordinates(.)[,2],) %>% 
  st_drop_geometry()

#Snap points to flow network
snap <- xy2segvert(x = pts$x, y = pts$y, rivers=flow_net)

#Define river distance between points!
output <- riverdistancemat(
  seg = snap$seg, 
  vert = snap$vert, 
  rivers = flow_net, 
  ID = pts$siteID)
outlet <- "SFM01"
output <- as_tibble((subset(output, colnames(output) %in% outlet))) %>%
  pivot_longer(., cols = everything(), names_to = "siteID", values_to = "distance_m") %>%
  arrange(., distance_m) %>%
  mutate_if(., is.numeric, round)

#Bind: Isotopes with distances
AIMS_isotopes <- left_join(AIMS_isotopes, output, by = "siteID")

#Read: GW levels
Mor4_6 <- as_tibble(read.csv("data/LTER/GW_Levels_4-6Mor.csv"))
Mor4_6 <- Mor4_6 %>%
  mutate(Date = mdy(Date)) %>%
  filter(Date >= as.Date("2020-10-01") & Date <= as.Date("2021-09-30")) %>%
  mutate(id = "Mor 4-6") %>%
  select(-c(Time, Date.Time))

LowerEis4_6 <- as_tibble(read.csv("data/LTER/GW_Levels_4-6Eis1.csv"))
LowerEis4_6 <- LowerEis4_6 %>%
  mutate(Date = mdy(Date)) %>%
  filter(Date >= as.Date("2020-10-01") & Date <= as.Date("2021-09-30")) %>%
  mutate(id = "Lower Eis 4-6") %>%
  select(-c(Time, Date.Time))

UpperEis4_6 <- as_tibble(read.csv("data/LTER/GW_Levels_4-6Eis2.csv"))
UpperEis4_6 <- UpperEis4_6 %>%
  mutate(Date = mdy(Date)) %>%
  filter(Date >= as.Date("2020-10-01") & Date <= as.Date("2021-09-30")) %>%
  mutate(id = "Upper Eis 4-6") %>%
  select(-c(Time, Date...Time))

#Plot!
months <- c("Jun" = "June", "Jul" = "July", "Aug" = "August")
AIMS_isotopes$month <- factor(AIMS_isotopes$month, levels = c( "Jun", "Jul", "Aug"))
AIMS_isotopes$d18OWaterBinned <- cut(AIMS_isotopes$d18OWater, breaks = c(seq(-6.5, -4.5, 0.25), 0.2))
AIMS_isotopes$flowing <- factor(AIMS_isotopes$flowing, levels = c("y", "n"))
ggplot()+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = min(UpperEis4_6$LEVEL), ymax = max(UpperEis4_6$LEVEL), fill = "#4DAF4A", alpha = 0.2)+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = min(LowerEis4_6$LEVEL), ymax = max(LowerEis4_6$LEVEL), fill = "#F781BF", alpha = 0.2)+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = min(Mor4_6$LEVEL), ymax = max(Mor4_6$LEVEL), fill = "#999999", alpha = 0.2)+
  geom_point(data = na.omit(AIMS_isotopes), aes(distance_m, Elevation_m, fill = d18OWaterBinned, shape = flowing), size=2)+
  facet_wrap(~month, labeller = as_labeller(months))+
  scale_fill_manual(values = c("#440154", "#46327e", "#365c8d", "#277f8e", "#1fa187", "#4ac16d", "#a0da39", "#fde725"))+
  scale_shape_manual(values = c(21,24), labels = c("Flowing", "Pooled"))+
  labs(x = "Distance to outlet (m)", y = "Elevation (m)")+
  labs(fill = expression(delta^{18}*"O (‰)"))+
  labs(shape = "")+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(strip.text = element_text(face = 'bold'))
ggsave(path = "figures/", "Fig4_GroundwaterLevels.png", dpi=300, width = 190, height = 90, units = "mm")