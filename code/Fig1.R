# Fig 1. Map of Konza Prairie
library(sf)
library(terra)
library(whitebox)
library(tidyverse)

source(file.path("code", "theme.R"))

dem <- file.path("data", "spatial", "dem", "10m.tif") %>% 
  rast()

crs <- crs(dem)

hillshade <- file.path("data", "spatial", "dem", "wbt", "2m_hillshade.tif") %>% 
  rast()

sites <- file.path("data", "AIMS", "AIMS_sampling_sites.csv") %>% 
  read_csv(show_col_types = F) %>% 
  mutate(
    ID = "Synoptic sampling sites"
  ) %>%
  add_row(ID = "USGS gage 06879650", lon = -96.59469, lat = 39.10207) %>%
  add_row(ID = "NEON tower", lon = -96.61294, lat = 39.11045) %>%
  add_row(ID = "NEON sampling site", lon = -96.60274, lat = 39.10449) %>% 
  add_row(ID = "Groundwater wells", lon = -96.584221, lat = 39.084785) %>%
  add_row(ID = "Meteorological station", lon = -96.608184, lat = 39.100785) %>%
  mutate(
    ID = factor(ID, levels = c("Synoptic sampling sites", "NEON sampling site", "NEON tower", "USGS gage 06879650", "Groundwater wells", "Meteorological station"))
  ) %>% 
  st_as_sf(coords = c("lon", "lat"), crs="epsg:4326", remove=F) %>% 
  st_transform(crs) %>% 
  select(ID, geometry)

watershed <- file.path("data", "spatial", "watersheds.gpkg") %>% 
  read_sf(layer = "watersheds") %>%
  filter(
    siteID == "USGS Gage 06879650"
  )

streams <- file.path("data", "spatial", "streams.gpkg") %>% 
  read_sf(layer = "streams") %>%
  rowwise() %>%
  filter(st_intersects(geom, watershed, sparse = FALSE)[1] | any(st_coordinates(geom)[,2] >= 4329000),
         st_intersects(geom, watershed, sparse = FALSE)[1] | any(st_coordinates(geom)[,1] <= 709000)
  ) %>%
  ungroup()

window <- st_bbox(
  st_union(
    st_buffer(sites, 200),
    watershed
  )
)

dem <- crop(dem, ext(window))
hillshade <- crop(hillshade, ext(window))

p1 = ggplot()+
  tidyterra::geom_spatraster(data = hillshade) +
  scale_fill_gradient(low="#000000", high="#FFFFFF", guide="none")+
  ggnewscale::new_scale_fill()+
  tidyterra::geom_spatraster(data = dem, alpha=0.8)+
  scale_fill_gradient(
    name = "Elevation (m)", 
    low="#00A600", 
    high="#EAB64E",
    guide = guide_colorbar(order = 1), 
    )+
  ggnewscale::new_scale_fill()+
  geom_sf(data = watershed, fill = NA, color = "#373737", linewidth=0.4)+
  geom_sf(data= streams, color = "#0000FF")+
  geom_sf(data = sites, aes(fill = ID, shape = ID), color = "#373737", size=3)+
  scale_fill_manual(
    name = "",
    values = c("#A84268", "#FCB97D", "#377eb8", "#FF5E5B", "#f781bf", "#984ea3"),
    guide = guide_legend(order = 2), 
    )+
  scale_shape_manual(
    name = "",
    values = c(21, 21, 22, 24, 25, 23),
    guide = guide_legend(order = 2), 
    )+
  #labs(x = "Easting [m]", y = "Northing [m]")+
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01, suffix = "\u00B0W"))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01, suffix = "\u00B0N"))+
  ggspatial::annotation_scale(location = "bl")+
  ggspatial::annotation_north_arrow(location = "tr", height = unit(8, "mm"), width = unit(8, "mm"))+
  coord_sf(
    expand = F,
    #datum = crs,
    xlim = c(window["xmin"], window["xmax"]),
    ylim = c(window["ymin"], window["ymax"]),
  )+
  theme(
    #axis.title = element_text(size = 14),
    #axis.text.y = element_text(hjust = 0.5, angle = 90)
  )

states <- tigris::states(cb = T) %>%
  filter(!STUSPS %in% c("AK", "HI", "PR", "VI", "GU", "MP", "AS")) %>%
  st_transform("epsg:5070")

pour_point <- tibble(
  ID = "KNZ",
  lon = -96.58724405,
  lat = 39.09225847
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs="epsg:4326", remove=F) %>%
  st_transform("epsg:5070")

p2 = ggplot()+
  geom_sf(data = states, fill=NA, color="#373737", lwd=0.4)+
  geom_sf(data = pour_point, fill="#E7315D", color="#373737", pch=21, size=3)+
  theme(
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p1
ggsave(file.path("figures", "Fig1a.png"), dpi=300, width=190, height=190, units="mm")

p2
ggsave(file.path("figures", "Fig1b.png"), dpi=300, width=110, height=110, units="mm")
