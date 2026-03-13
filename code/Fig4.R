# Fig 4. Variation in δ18O with distance to outlet 
library(sf)
library(riverdist)
library(tidyverse)

source(file.path("code", "theme.R"))

start <- as.Date("2020-10-01")
end <- as.Date("2021-09-30")

streams <- file.path("data", "spatial", "streams.gpkg") %>% 
  read_sf(layer = "headwater_streams")

crs <- st_crs(streams)

AIMS_isotopes <- file.path("data", "AIMS", "AIMS_isotopes.csv") %>% 
  read_csv(show_col_types = F) %>%
  mutate(
    month = factor(month, levels = 6:8, labels = c("June", "July", "August")),
    flowing = factor(flowing, levels = c(T,F), labels = c("Flowing", "Pooled")),
    d18O_binned = cut(d18O, breaks = c(seq(-6.5, -4.5, 0.3), 0.2))
  ) %>% 
  st_as_sf(coords = c("lon", "lat"), crs="epsg:4326", remove=F) %>% 
  st_transform(crs) 

sites <- file.path("data", "AIMS", "AIMS_sampling_sites.csv") %>% 
  read_csv(show_col_types = F) %>% 
  st_as_sf(coords = c("lon", "lat"), crs="epsg:4326", remove=F) %>% 
  st_transform(crs) 

# prep flownet
streams <- streams %>% 
  mutate(
    UID = row_number(),
    type = "RIVER"
  ) %>% 
  select(-geom, geom) %>% 
  rename(geometry = geom)

# flownet
flownet <- line2network(streams, tolerance=1)

# snap sites to flownet
coords <- st_coordinates(sites)
sites_snapped <- xy2segvert(x = coords[,1], y = coords[,2], rivers = flownet)

# river distance between points
river_distance <- riverdistancemat(
  seg = sites_snapped$seg,
  vert = sites_snapped$vert,
  rivers = flownet,
  ID = sites$siteID
)

outlet <- "SFM01"

distances <- subset(river_distance, colnames(river_distance) %in% outlet) %>%
  as_tibble() %>% 
  pivot_longer(cols = everything(), names_to = "siteID", values_to = "distance_m") %>% 
  mutate(
    distance_m = round(distance_m, digits=2)
  )

AIMS_isotopes <- left_join(AIMS_isotopes, distances, by = "siteID")

# top of Cottonwood limestone is 04M03 (field notes)
ref <- filter(sites, siteID == "04M03")[, "elevation_m"][[1]]

stratigraphy <- file.path("data", "Konza_stratigraphy.csv") %>% 
  read_csv(show_col_types = F) %>% 
  mutate(
    order = row_number()
  ) %>%
  arrange(-order) %>% 
  mutate(
    cumulative_thickness_m = cumsum(thickness_m),
    top = cumulative_thickness_m + ref - 17.0,
    bottom = top - thickness_m
  ) %>%
  select(order, member, type, thickness_m, cumulative_thickness_m, top, bottom) 

limestone_units <- filter(stratigraphy, type == "Ls") 

ggplot()+
  geom_rect(data = limestone_units, aes(ymax=top, ymin=bottom, xmax=Inf, xmin=-Inf), fill = "gray95")+
  geom_point(data = AIMS_isotopes, aes(distance_m, elevation_m, fill = d18O_binned, shape = flowing), color="373737", size=2.5)+
  scale_fill_manual(
    name = expression(delta^{18}*"O (‰)"),
    values = c("#440154", "#38588C", "#25858E", "#2AB07F", "#86D549", "#FDE725"),
    labels = c("-6.2 to -5.9", "-5.9 to -5.6", "-5.6 to -5.3", "-5.3 to -5.0", "-5.0 to -4.7", "-4.7 to  0.2"),
    guide = guide_legend(override.aes = list(shape = 21), order = 1)
  )+
  scale_shape_manual(
    name = "",
    values = c(
      "Flowing" = 21, 
      "Pooled" = 24
      )
  )+
  labs(x = "Distance to outlet (m)", y = "Elevation (m)")+
  facet_wrap(~month)+
  theme(
    #strip.text = element_text(face = "bold"),
    strip.text = element_text(size = 12)
  )
ggsave(file.path("figures", "Fig4.png"), dpi=300, width=190, height=90, units="mm")