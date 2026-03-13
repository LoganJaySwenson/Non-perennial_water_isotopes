# Fig S4. Spatial variation in δ18O during the summer dry-down period with limestone units
library(sf)
library(terra)
library(tidyverse)

source(file.path("code", "theme.R"))

dem <- file.path("data", "spatial", "dem", "10m.tif") %>% 
  rast()

crs <- crs(dem)

AIMS_isotopes <- file.path("data", "AIMS", "AIMS_isotopes.csv") %>% 
  read_csv(show_col_types = F) %>%
  mutate(
    month = factor(month, levels = 6:8, labels = c("June", "July", "August")),
    d18O_binned = cut(d18O, breaks = c(seq(-6.5, -4.5, 0.3), 0.2))
  ) %>% 
  st_as_sf(coords = c("lon", "lat"), crs="epsg:4326", remove=F) %>% 
  st_transform(crs) 

sites <- file.path("data", "AIMS", "AIMS_sampling_sites.csv") %>% 
  read_csv(show_col_types = F) %>% 
  st_as_sf(coords = c("lon", "lat"), crs="epsg:4326", remove=F) %>% 
  st_transform(crs) 

streams <- file.path("data", "spatial", "streams.gpkg") %>% 
  read_sf(layer = "headwater_streams") %>% 
  st_transform(crs)

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

dem <- crop(dem, ext(st_buffer(streams, dist=200)))

dem_based_limestone_units <- dem %>% 
  terra::as.data.frame(xy=T) %>%
  as_tibble() %>%
  select(
    "elevation_m" = "GIS203",
    "x", "y"
    ) %>%
  mutate(
    elevation_m_binned = case_when(
      elevation_m <= 429.02 & elevation_m >= 423.47 ~ "Florence Ls",
      elevation_m <= 416.52 & elevation_m >= 412.57 ~ "Kinney Ls",
      elevation_m <= 412.57 & elevation_m >= 405.80 ~ "Wymore Ls",
      elevation_m <= 400.92 & elevation_m >= 398.44 ~ "Threemile Ls",
      elevation_m <= 394.08 & elevation_m >= 391.23 ~ "Funston Ls",
      elevation_m <= 385.12 & elevation_m >= 382.68 ~ "Crouse Ls",
      elevation_m <= 378.32 & elevation_m >= 377.34 ~ "Middleburg Ls",
      elevation_m <= 375.16 & elevation_m >= 373.22 ~ "Eiss Ls",
      elevation_m <= 369.07 & elevation_m >= 368.31 ~ "Morrill Ls",
      elevation_m <= 364.96 & elevation_m >= 363.21 ~ "Cottonwood Ls",
      elevation_m <= 355.06 & elevation_m >= 350.82 ~ "Neva Ls"
      )
    ) %>%
  drop_na(elevation_m_binned)

ggplot()+
  geom_raster(data = dem_based_limestone_units, aes(x, y), fill = "#FEF3A3")+
  ggnewscale::new_scale_fill()+
  geom_sf(data= streams, color = "#0000FF")+
  geom_sf(data = AIMS_isotopes, aes(fill = d18O_binned), color = "#373737", pch=21, size=2.5)+
  scale_fill_manual(
    name = expression(delta^{18}*"O (‰)"),
    values = c("#440154", "#38588C", "#25858E", "#2AB07F", "#86D549", "#FDE725"),
    labels = c("-6.2 to -5.9", "-5.9 to -5.6", "-5.6 to -5.3", "-5.3 to -5.0", "-5.0 to -4.7", "-4.7 to  0.2")
  )+
  labs(x = "", y = "")+
  facet_wrap(~month)+
  ggspatial::annotation_scale(data = filter(AIMS_isotopes, month == "August"), location = "br")+
  coord_sf()+
  theme(
    #strip.text = element_text(face = "bold"),
    panel.border = element_blank(),
    strip.text = element_text(size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
ggsave(file.path("figures", "FigS4.png"), dpi=300, width=190, height=60, units="mm")
