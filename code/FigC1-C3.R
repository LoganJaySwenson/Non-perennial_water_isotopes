# Fig C1, C2, & C3. Affects of different assumed priors on estimated young water fractions
library(sf)
library(cowplot)
library(ggtext)
library(tidyverse)

source("code/Theme+Settings.R")

priors <- c("c(0.8,0.2)", "c(0.5,0.5)", "c(0.2,0.8)", "c(0.153,0.847)")

spatial_priors <- read_csv("data/spatial_priors.csv") %>%
  mutate(
    prior = str_replace_all(prior, " ", ""),
    prior = factor(prior, levels = c("c(0.8,0.2)", "c(0.5,0.5)", "c(0.2,0.8)", "c(0.153,0.847)")),
    month = factor(month, levels = c("Jun", "Jul", "Aug"))
  ) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)

spatial_summary <- spatial_priors %>%
  st_drop_geometry() %>%
  group_by(prior, month) %>%
  summarise(
    mean_s1 = mean(s1, na.rm = T),
    sd_s1 = sd(s1, na.rm = T),
    label = paste0(round(mean_s1, 1), "% ± ", round(sd_s1, 1)),
    prior_label = str_replace_all(first(prior), "c\\(([^,]+),([^)]+)\\)", "Prior: \\1/\\2"),
  )

crs <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
streams <- st_read("data/spatial/stream_network/Konza_stream_network.shp", quiet = T)
st_crs(streams) <- crs

pp <- list()

for (i in seq_along(priors)){
  
  spatial_prior <- filter(spatial_priors, prior == priors[i])
  spatial_stats <- filter(spatial_summary, prior == priors[i])
  
  min <- min(spatial_prior$s1, na.rm = T)
  max <- max(spatial_prior$s1, na.rm = T)
  
  label <- paste0(
    "June: ", filter(spatial_stats, month == "Jun") %>% .$label, "\n",
    "July: ", filter(spatial_stats, month == "Jul") %>% .$label, "\n",
    "August: ", filter(spatial_stats, month == "Aug") %>% .$label, "\n"
  )
  
  p =
    ggplot() +
    geom_sf(data = streams, color = "blue") +
    geom_sf(data = spatial_prior, aes(fill = s1), color = "#373737", pch = 21, size = 2.5) +
    facet_wrap(~month, labeller = labeller(month = c("Jun" = "June", "Jul" = "July", "Aug" = "August"))) +
    labs(x = "", y = "", 
         title =  spatial_stats$prior_label,
         #subtitle = label # placement done in inkscape
    ) +
    scale_fill_gradient(
      name = expression(F["yw"]* " (%)"),
      breaks = c(min, max),
      labels = scales::label_number()(c(min, max)),
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
  
  pp[[i]] <- p
  
}

plot_grid(plotlist = pp, nrow = 4, ncol = 1)
# ggsave(path = "figures/", "C1_no_labels.svg", dpi=600, width = 190, height = 190, units = "mm")
# ggsave(path = "figures/", "C1_labeled.svg", dpi=600, width = 190, height = 190, units = "mm")




temporal_priors <- read_csv("data/temporal_priors.csv") %>%
  # filter to 2021 water year
  filter(date >= as.Date("2020-10-01") & date <= as.Date("2021-09-30")) %>%
  mutate(
    prior = str_replace_all(prior, " ", ""),
    prior = factor(prior, levels = c("c(0.8,0.2)", "c(0.5,0.5)", "c(0.2,0.8)", "c(0.153,0.847)")),
  )

temporal_summary <- temporal_priors %>%
  group_by(prior) %>%
  summarise(
    mean_s1 = mean(s1, na.rm = T),
    sd_s1 = sd(s1, na.rm = T),
    label = paste0(round(mean_s1, 1), "% ± ", round(sd_s1, 1))
  )

ggplot()+
  geom_point(data = temporal_priors, aes(date, s1, fill = prior), color = "#373737", pch=21, size=3.5)+
  geom_hline(data = temporal_summary, aes(yintercept = mean_s1, color = prior), linetype = "dashed", lwd=0.8) +
  geom_text(data = temporal_summary, aes(as.Date("2020-11-15"), mean_s1+2.5, label = label), color = "#000000",
            hjust = 0, size = 10/.pt) +
  labs(x = "", y = expression(F["yw"]* " (%)"), fill = "Prior")+
  scale_fill_manual(
    values = c(
      "c(0.8,0.2)" = "#deebf7",
      "c(0.5,0.5)" = "#9ecae1",
      "c(0.2,0.8)" = "#4292c6",
      "c(0.153,0.847)" = "#084594"
    ),
    labels = c(
      "c(0.8,0.2)" = "0.8/0.2",
      "c(0.5,0.5)" = "0.5/0.5",
      "c(0.2,0.8)" = "0.2/0.8",
      "c(0.153,0.847)" = "0.153/0.847"
    )
  ) +
  scale_color_manual(
    values = c(
      "c(0.8,0.2)" = "#deebf7",
      "c(0.5,0.5)" = "#9ecae1",
      "c(0.2,0.8)" = "#4292c6",
      "c(0.153,0.847)" = "#084594"
    ),
    guide = "none" 
  ) +
  scale_x_date(breaks=scales::date_breaks("months"), labels=scales::date_format("%b %y"))+
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
  )
# ggsave(path = "figures/", "C2.png", dpi=300, width = 190, height = 110, units = "mm")
# ggsave(path = "figures/", "C2.pdf", dpi=300, width = 190, height = 110, units = "mm")




# Pairwise t-test to check if differences in means of s1 between months per prior are equal to 0.
# p-value < 0.05 indicates a significant difference.

pairwise_t <- list()

for (p in seq_along(priors)) {
  
  spatial_prior <- spatial_priors %>% filter(prior == priors[p])
  
  t <- pairwise.t.test(
    spatial_prior$s1,
    spatial_prior$month
  )
  
  pairwise_t[[p]] <- tibble(
    prior = priors[p],
    t1_pval = t$p.value["Jul", "Jun"],
    t1_significance = if_else(t1_pval < 0.05, T, F),
    t2_pval = t$p.value["Aug", "Jul"],
    t2_significance = if_else(t2_pval < 0.05, T, F),
    t3_pval = t$p.value["Aug", "Jun"],
    t3_significance = if_else(t3_pval < 0.05, T, F)
  )
}

pairwise_t <- bind_rows(pairwise_t)

pp <- list()

for (i in seq_along(priors)) {
  
  spatial_prior <- filter(spatial_priors, prior == priors[i])
  spatial_stats <- filter(spatial_summary, prior == priors[i])
  
  t <- filter(pairwise_t, prior == priors[i])
  
  format_p <- function(pval, sig) {
    pval_str <- format(pval, scientific = T, digits = 2)
    if (sig) {
      paste0("<b>p = ", pval_str, "*</b>")
    } else {
      paste0("p = ", pval_str)
    }
  }
  
  label <- paste0(
    "June to July: ", format_p(t$t1_pval, t$t1_significance), "<br>",
    "July to August: ", format_p(t$t2_pval, t$t2_significance), "<br>",
    "June to August: ", format_p(t$t3_pval, t$t3_significance)
  )
  
  p <- ggplot() +
    geom_boxplot(data = spatial_prior, aes(x = month, y = s1), fill = "#BBBBBB", width = 0.5, outlier.shape = 21, outlier.size = 2) +
    annotate("richtext", x = 0.5, y = -Inf, label = label,
             hjust = 0, vjust = 0, size = 9/.pt,
             fill = "transparent", label.color = NA,
             color = "#373737", label.padding = unit(0.5, "lines")) +
    labs(
      x = "",
      y = expression(F["yw"]* " (%)"),
      title = spatial_stats$prior_label
    ) +
    scale_x_discrete(labels = c("June", "July", "August"))
  
  pp[[i]] <- p
}

plot_grid(plotlist = pp, nrow = 2, ncol = 2)
# ggsave(path = "figures/", "C3.png", dpi=300, width = 150, height = 150, units = "mm")
# ggsave(path = "figures/", "C3.pdf", dpi=300, width = 150, height = 150, units = "mm")
