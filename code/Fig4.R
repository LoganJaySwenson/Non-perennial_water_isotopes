# Fig 4. Variation in δ18O with distance to outlet 
library(sf)
library(raster)
library(viridis)
library(riverdist)
library(lubridate)
library(tidyverse)

# Publication theme
source("code/Theme+Settings.R")

date_start <- as.Date("2020-10-01")
date_end <- as.Date("2021-09-30")

crs <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Read: AIMS isotopes
AIMS_isotopes <- read_csv("data/AIMS_isotopes_RF.csv") %>%
  select(1:18) %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4326)

# Read: Sampling locations
pts <- read_csv("data/AIMS/AIMS_STIC_details.csv") %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs = crs) # reproject vector to planar coordinates

# Read: Stream network
streams <- st_read("data/spatial/stream_network/Konza_stream_network.shp")
st_crs(streams) <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" # reproject vector to planar coordinates

# Prep flownet
streams <- streams %>% 
  # remove z information (to make it a 'linear' feature)
  st_zm() %>% 
  # add ID 
  mutate(
    uid = seq(1, nrow(.)),
    type = 'river')

# export streams to spatial file
st_write(streams, dsn = "data/spatial/distance_to_outlet/streams.shp", layer = "streams", append = F)

# create flownet
flow_net <- line2network(path = "data/spatial/distance_to_outlet/streams.shp", layer="streams", tolerance = 1)

# save flow network 
save(flow_net, file = "data/spatial/distance_to_outlet/streams.rda")

# prep sampling locations
pts <- pts %>% 
  # define coordinates
  mutate(
    x = st_coordinates(.)[,1],
    y = st_coordinates(.)[,2],) %>% 
  st_drop_geometry()

# Snap points to flow network
snap <- xy2segvert(x = pts$x, y = pts$y, rivers=flow_net)

# Define river distance between points
output <- riverdistancemat(
  seg = snap$seg, 
  vert = snap$vert, 
  rivers = flow_net, 
  ID = pts$siteID)
outlet <- "SFM01"
output <- as_tibble((subset(output, colnames(output) %in% outlet))) %>%
  pivot_longer(cols = everything(), names_to = "siteID", values_to = "distance_m") %>%
  arrange(distance_m) %>%
  mutate_if(is.numeric, round)

# Bind: Isotopes with distances to outlet
AIMS_isotopes <- left_join(AIMS_isotopes, output, by = "siteID")

# Read: GW levels
Mor4_6 <- read_csv("data/LTER/GW_Levels_4-6Mor.csv", show_col_types=F) %>%
  rename_with(.cols = everything(), tolower) %>%
  mutate(date = mdy(date), 
         id = "Mor 4-6") %>%
  filter(date >= date_start & date <= date_end) %>%
  select(date, id, level)

LowerEis4_6 <- read_csv("data/LTER/GW_Levels_4-6Eis1.csv", show_col_types=F) %>%
  rename_with(.cols = everything(), tolower) %>%
  mutate(date = mdy(date), 
         id = "Lower Eis 4-6") %>%
  filter(date >= date_start & date <= date_end) %>%
  select(date, id, level)

UpperEis4_6 <- read_csv("data/LTER/GW_Levels_4-6Eis2.csv", show_col_types=F) %>%
  rename_with(.cols = everything(), tolower) %>%
  mutate(date = mdy(date), 
         id = "Upper Eis 4-6") %>%
  filter(date >= date_start & date <= date_end) %>%
  select(date, id, level)

# Top of Cottonwood Limestone is 04M03 (field notes)
sites <- read_csv("data/AIMS/AIMS_STIC_details.csv")
ref <- round(sites %>% filter(siteID == "04M03") %>% pull(Elevation_m), digits = 2)

# Limestone/shale members
limestones <- read_csv("data/LTER/Konza_GW_Well_Info.csv") %>%
  mutate(order = 1:21) %>%
  arrange(-order) %>%
  mutate(cumulative = cumsum(thickness_m),
         top = cumulative + ref - 17.0,
         bottom = top - thickness_m) %>%
  filter(type == "Ls")

# Plot!
months <- c("Jun" = "June", "Jul" = "July", "Aug" = "August")
AIMS_isotopes$month <- factor(AIMS_isotopes$month, levels = c( "Jun", "Jul", "Aug"))
AIMS_isotopes$d18OWaterBinned <- cut(AIMS_isotopes$d18OWater, breaks = c(seq(-6.5, -4.5, 0.3), 0.2))
AIMS_isotopes$flowing <- factor(AIMS_isotopes$flowing, levels = c("y", "n"))
ggplot()+
  geom_rect(data = limestones, aes(ymax = top, ymin = bottom, xmax = Inf, xmin = -Inf), fill = "#999999", alpha = 0.2)+
  ggnewscale::new_scale_fill()+
  geom_point(data = na.omit(AIMS_isotopes), aes(distance_m, Elevation_m, fill = d18OWaterBinned, shape = flowing), size=2)+
  facet_wrap(~month, labeller = as_labeller(months))+
  scale_fill_manual(values = c("#440154", "#38588c", "#25858e", "#2ab07f", "#86d549", "#fde725"),
                    labels = c("-6.2 to -5.9", "-5.9 to -5.6", "-5.6 to -5.3", "-5.3 to -5.0", "-5.0 to -4.7", "-4.7 to  0.2"))+
  scale_shape_manual(values = c(21,24), labels = c("Flowing", "Pooled"))+
  labs(x = "Distance to outlet (m)", y = "Elevation (m)")+
  labs(fill = expression(delta^{18}*"O (‰)"))+
  labs(shape = "")+
  scale_y_continuous(limits = c(350,410), breaks = seq(350, 410, by = 10))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1))+
  theme(strip.text = element_text(face = 'bold'))
ggsave(path = "figures/", "Fig4.png", dpi=300, width = 190, height = 90, units = "mm")
ggsave(path = "figures/", "Figure4.pdf", dpi=300, width = 190, height = 90, units = "mm")