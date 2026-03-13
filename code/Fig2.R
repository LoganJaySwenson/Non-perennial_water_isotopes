# Fig 2. Konza Prairie conditions for 2021 Water Year
library(patchwork)
library(tidyverse)

source(file.path("code", "theme.R"))

min_q <- 0.001
start <- as.Date("2020-10-01")
end <- as.Date("2021-09-30")

# Precip at HQ weather station
precip <- file.path("data", "LTER", "LTER_HQ_weather_station.csv") %>%
  read_csv(show_col_types = F, col_select = c(siteID, date, precip_mm)) %>%
  filter(date >= start & date <= end) 

# Streamflow at USGS 06879650
flow <- file.path("data", "USGS", "USGS_06879650.csv") %>% 
  read_csv(show_col_types = F) %>% 
  filter(date >= start & date <= end) 

water_levels <- file.path("data", "LTER", "LTER_4-6wells.csv") %>% 
  read_csv(show_col_types = F) %>% 
  mutate(
    siteID = factor(siteID, levels = c("Upper Eis 4-6", "Lower Eis 4-6", "Mor 4-6"))
  ) %>% 
  filter(date >= start & date <= end) 

p1 = ggplot()+
  geom_col(data = precip, aes(date, precip_mm), color = "#984EA3")+
  geom_vline(xintercept = as.Date("2021-06-07"), color = "#E7315D", linetype = "longdash")+
  geom_vline(xintercept = as.Date("2021-07-13"), color = "#E7315D", linetype = "longdash")+
  geom_vline(xintercept = as.Date("2021-08-09"), color = "#E7315D", linetype = "longdash")+
  labs(x = "", y = "Precipitation (mm/d)")+
  scale_x_date(
    expand = c(0, 0),
    breaks = scales::date_breaks("months"), 
    labels = scales::date_format("%b %y")
    )+
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  )+
  theme(
    axis.title = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p2 = ggplot()+
  geom_line(data = flow, aes(date, if_else(flow_m3s < min_q, min_q, flow_m3s)), color = "#377EB8", lwd=0.7)+
  geom_vline(xintercept = as.Date("2021-06-07"), color = "#E7315D", linetype = "longdash")+
  geom_vline(xintercept = as.Date("2021-07-13"), color = "#E7315D", linetype = "longdash")+
  geom_vline(xintercept = as.Date("2021-08-09"), color = "#E7315D", linetype = "longdash")+
  labs(x = "", y = expression(Streamflow ~ (m^3*"/"*s)))+
  scale_x_date(
    expand = c(0, 0),
    breaks = scales::date_breaks("months"), 
    labels = scales::date_format("%b %y")
    )+
  scale_y_log10(
    expand = c(0,0),
    limits = c(0.001, 10),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    )+
  annotation_logticks(sides = "l")+
  theme(
    axis.title = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p3 = ggplot()+
  geom_line(data = water_levels, aes(datetime, water_level_masl, color = siteID), lwd=0.5)+
  geom_vline(xintercept = as.Date("2021-06-07"), color = "#E7315D", linetype = "longdash")+
  geom_vline(xintercept = as.Date("2021-07-13"), color = "#E7315D", linetype = "longdash")+
  geom_vline(xintercept = as.Date("2021-08-09"), color = "#E7315D", linetype = "longdash")+
  annotate("text", as.Date("2020-10-02"), 371.2, label = "Upper Eiss", color = "#4DAF4A", size = 10/.pt, hjust=0)+
  annotate("text", as.Date("2020-10-02"), 369.4, label = "Lower Eiss", color = "#F781BF", size = 10/.pt, hjust=0)+
  annotate("text", as.Date("2020-10-02"), 364.2, label = "Morrill", color = "#999999", size = 10/.pt, hjust=0)+
  labs(x = "", y = "Groundwater level (m.a.s.l)")+
  scale_color_manual(
    name = "",
    values = c(
      "Upper Eis 4-6" = "#4DAF4A",
       "Lower Eis 4-6" = "#F781BF",
       "Mor 4-6" = "#999999"
      ),
    labels = c(
      "Upper Eis 4-6" = "Upper Eiss",
      "Lower Eis 4-6" = "Lower Eiss",
      "Mor 4-6" = "Morrill"
    )
  )+
  guides(color = "none")+
  scale_x_date(
    expand = c(0, 0),
    breaks = seq(as.Date("2020-11-01"),
                 max(water_levels$datetime),
                 by = "1 month"),
    labels = scales::date_format("%b %Y")
  )+
  scale_y_continuous(
    breaks = seq(364, 372, 2)
  )+
  theme(
    axis.title = element_text(size = 11),
    legend.position = "bottom",
    axis.text.x = element_text(angle=45, hjust=1)
    
  )

pp <- p1 / p2 / p3
pp

ggsave(file.path("figures", "Fig2.png"), dpi=300, width=190, height=190, units="mm")
