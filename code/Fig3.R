# Fig 3. Spatial variation in δ18O during summer dry-down period
library(sf)
library(tidyverse)

source(file.path("code", "theme.R"))

crs <- "epsg:5070"

AIMS_isotopes <- file.path("data", "AIMS", "AIMS_isotopes.csv") %>% 
  read_csv(show_col_types = F) %>%
  mutate(
    month = factor(month, levels = 6:8, labels = c("June", "July", "August")),
    d18O_binned = cut(d18O, breaks = c(seq(-6.5, -4.5, 0.3), 0.2))
  ) %>% 
  st_as_sf(coords = c("lon", "lat"), crs="epsg:4326", remove=F) %>% 
  st_transform(crs) 

streams <- file.path("data", "spatial", "streams.gpkg") %>% 
  read_sf(layer = "headwater_streams") %>% 
  st_transform(crs)

ggplot()+
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
ggsave(file.path("figures", "Fig3.png"), dpi=300, width=190, height=60, units="mm")
