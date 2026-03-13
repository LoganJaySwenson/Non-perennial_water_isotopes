# Konza watershed delineation
library(sf)
library(terra)
library(whitebox)
library(tidyverse)

source(file.path("code", "theme.R"))

save_dir <- file.path("data", "spatial", "dem", "wbt")
dir.create(save_dir, showWarnings=F)

dem <- file.path("data", "spatial", "dem", "10m.tif") %>% 
  rast()

crs <- crs(dem)

wbt_hillshade(
  file.path("data", "spatial", "dem", "10m.tif"),
  file.path(save_dir, "10m_hillshade.tif"),
  azimuth=315,
  altitude=30
  )

wbt_hillshade(
  file.path("data", "spatial", "dem", "2m.tif"),
  file.path(save_dir, "2m_hillshade.tif"),
  azimuth=315,
  altitude=30
)

# fill depressions using Lindsay's (2016) algorithm
wbt_breach_depressions(
  file.path("data", "spatial", "dem", "10m.tif"),
  file.path(save_dir, "10m_filled.tif")
  )

# d8 flow direction pointer
wbt_d8_pointer(
  file.path(save_dir, "10m_filled.tif"),
  file.path(save_dir, "10m_pointer.tif")
  )

# d8 flow accumulation 
wbt_d8_flow_accumulation(
  file.path(save_dir, "10m_filled.tif"),
  file.path(save_dir, "10m_accumulation.tif")
  )

pour_points <- tibble(
  siteID = c("South Fork Kings Creek", "NEON sampling site", "USGS Gage 06879650"),
  lon = c(-96.588494, -96.60274, -96.59469),
  lat = c(39.093064, 39.10449, 39.10207)
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs="epsg:4326", remove=F) %>%
  st_transform(crs)

st_write(pour_points, file.path(save_dir, "pour_points.shp"), delete_dsn=T, quiet=T)

# delineate streams
wbt_extract_streams(
  file.path(save_dir, "10m_accumulation.tif"),
  file.path(save_dir, "10m_streams.tif"),
  threshold = 1500
  )

wbt_raster_streams_to_vector(
  file.path(save_dir, "10m_streams.tif"),
  file.path(save_dir, "10m_pointer.tif"),
  file.path(save_dir, "streams.shp"),
)

# snap pour points to streams
wbt_jenson_snap_pour_points(
  file.path(save_dir, "pour_points.shp"),
  file.path(save_dir, "10m_streams.tif"),
  file.path(save_dir, "pour_points_snapped.shp"),
  snap_dist = 100
  )

# delineate watershed
wbt_watershed(
  file.path(save_dir, "10m_pointer.tif"),
  file.path(save_dir, "pour_points_snapped.shp"),
  file.path(save_dir, "10m_watersheds.tif")
)

streams <- file.path("data", "spatial", "dem", "wbt", "streams.shp") %>% 
  read_sf() %>%
  st_set_crs(crs)

watersheds <- file.path("data", "spatial", "dem", "wbt", "10m_watersheds.tif") %>% 
  rast() %>%
  as.polygons(dissolve=T) %>%
  st_as_sf() %>% 
  slice(c(1,3,2)) %>%
  mutate(
    ID = row_number(),
  ) %>%
  select(ID, geometry)

geometries <- vector("list", nrow(watersheds))
for (i in 1:nrow(watersheds)){
  if (i == 1){
    geometries[[i]] <- watersheds$geometry[[i]]
  } else {
    geometries[[i]] <- st_union(watersheds$geometry[1:i])[[1]]
  }
}
st_geometry(watersheds) <- st_sfc(geometries, crs=crs)

watersheds <- watersheds %>%
  mutate(
    siteID = c("South Fork Kings Creek", "USGS Gage 06879650", "NEON sampling site"),
    area_km2 = round(as.double(st_area(geometry)) /1e6, digits=2)
  ) %>% 
  select(-geometry, geometry)

st_write(watersheds, file.path("data", "spatial", "watersheds.gpkg"), layer="watersheds", delete_layer=T)
st_write(streams, file.path("data", "spatial", "streams.gpkg"), layer="streams", delete_layer=T)