# Fig 8. Bayesian-estimated young water fractions in headwaters with different priors
library(sf)
library(cowplot)
library(tidyverse)

source(file.path("code", "theme.R"))

crs <- "epsg:5070"

spatial_priors <- file.path("data", "spatial_priors.csv") %>% 
  read_csv(show_col_types = F) %>%
  mutate(
    month = factor(month, levels = 6:8, labels = c("June", "July", "August")),
    prior = factor(prior, 
                   levels = c("c(0.5, 0.5)", "c(0.145, 0.855)", "c(0.038, 0.962)")),
    s1 = s1*100,
    s2 = s2*100
  ) %>% 
  st_as_sf(coords = c("lon", "lat"), crs="epsg:4326", remove=F) %>% 
  st_transform(crs) 

priors <- levels(spatial_priors$prior)

streams <- file.path("data", "spatial", "streams.gpkg") %>% 
  read_sf(layer = "headwater_streams") %>% 
  st_transform(crs)

spatial_summary <- spatial_priors %>%
  st_drop_geometry() %>%
  group_by(prior, month) %>%
  summarise(
    mean_s1 = mean(s1, na.rm=T),
    sd_s1 = sd(s1, na.rm=T),
    min_s1 = min(s1, na.rm=T),
    max_s1 = max(s1, na.rm=T),
    label = paste0(sprintf("%.1f", mean_s1), "% ± ", sprintf("%.1f", sd_s1)),
    .groups = "drop"
  ) %>%
  group_by(prior) %>%
  summarise(
    label = paste0(
      "June: ", label[month == "June"], "\n",
      "July: ", label[month == "July"], "\n",
      "August: ", label[month == "August"]
      ),
    .groups = "drop"
  )

titles <- list(
  expression("Uninformed (0.5/0.5)"),
  expression("Flow-weighted" ~ F["yw"] ~ "(0.15/0.85)"),
  expression("Time-weighted" ~ F["yw"] ~ "(0.04/0.96)")
)

pp <- list()
for (i in seq_along(priors)){
  
  spatial_prior <- filter(spatial_priors, prior == priors[i])

  vmin <- min(spatial_prior$s1, na.rm=T)
  vmax <- max(spatial_prior$s1, na.rm=T)
  label <- filter(spatial_summary, prior == priors[i])[["label"]]

  p = ggplot() +
    geom_sf(data= streams, color = "#0000FF")+
    geom_sf(data = spatial_prior, aes(fill = s1), color = "#373737", pch = 21, size = 2.5) +
    facet_wrap(~month) +
    labs(x = "", y = "", 
         title =  titles[[i]],
         subtitle = label # placement done in inkscape
    ) +
    scale_fill_gradient(
      name = expression(F["yw"]* " (%)"),
      breaks = c(vmin, vmax),
      labels = c(sprintf("%.1f", vmin), sprintf("%.1f", vmax)),
    ) +
    coord_sf()+
    theme(
      panel.border = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.subtitle = element_text(size = 8), 
      legend.key.height = unit(0.2, "cm"),
      legend.key.width = unit(0.1, "cm"),
    )
  
  if (i == 3)
    p <- p + ggspatial::annotation_scale(data = filter(spatial_prior, month == "August"), location = "br")
  
  pp[[i]] <- p
  
}

plot_grid(plotlist=pp, ncol=1, nrow=3)
#ggsave(file.path("figures", "Fig8.png"), dpi=300, width=190, height=190, units="mm")
ggsave(file.path("figures", "Fig8_annotated.png"), dpi=300, width=190, height=190, units="mm")
