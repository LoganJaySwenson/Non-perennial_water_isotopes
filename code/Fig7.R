# Fig 7. Temporal- % of streamwater less than ~3 months in age. 
# streamwater age is defined as the mean of the posterior distribution of source mixtures.
library(scales)
library(viridis)
library(lubridate)
library(tidyverse)

# Publication theme
source("code/Theme+Settings.R")

date_start <- as.Date("2020-10-01")
date_end <- as.Date("2021-09-30")

# Read: Streamflow at Kings Creek (USGS 06879650), Manhattan, KS
Flow <- dataRetrieval::readNWISdv(siteNumber = "06879650", parameterCd = "00060", startDate = "", endDate = "") %>%
  dataRetrieval::renameNWISColumns() %>%
  rename_with(.cols = everything(), tolower) %>%
  filter(date >= date_start & date <= date_end) %>%
  mutate(flow = (flow * 0.028316847)) #cfs to metric

# Identify start and end of flowing periods
dates_flow <- Flow$date[Flow$flow > 0]
dates_diffs <- c(999, diff(dates_flow))  # subtract dates from previous date; put first flow as 999
dates_starts <- dates_flow[which(dates_diffs != 1)]
dates_ends <- c(dates_flow[which(dates_diffs != 1) - 1], max(dates_flow))
df_flow_periods <- tibble(first_flow_date = dates_starts,
                          last_flow_date = dates_ends,
                          total_flow_days = as.numeric(dates_ends - dates_starts) + 1)

# NEON δ18O timeseries
NEON_isotopes <- read_csv("data/NEON_isotopes_results.csv") %>%
  mutate(date = ymd(date)) %>%
  filter(date >= date_start & date <= date_end) %>% #filter 2021 WY
  na.omit()

# Read: AIMS isotopes
AIMS_isotopes <- read_csv("data/AIMS_isotopes_results.csv") %>%
  mutate(date = ymd(date)) %>%
  mutate(date = replace(date, date == "2021-06-06", "2021-06-07"))

# Plot!
NEON_isotopes$s1Binned <- cut(NEON_isotopes$s1, breaks = c(seq(from = 36, to = 64, by = 4)))
AIMS_isotopes$s1Binned <- cut(AIMS_isotopes$s1, breaks = c(seq(from = 36, to = 64, by = 4)))
ggplot()+
  #geom_rect(aes(xmin = as.Date("2021-04-13"), xmax = as.Date("2021-07-07"), ymin = -Inf, ymax = Inf), fill = "#377EB8", alpha = 0.2, inherit.aes = F)+
  #geom_rect(aes(xmin = as.Date("2021-07-15"), xmax = as.Date("2021-07-24"), ymin = -Inf, ymax = Inf), fill = "#377EB8", alpha = 0.2, inherit.aes = F)+
  geom_smooth(data = NEON_isotopes, aes(date, s1,), alpha = 1, fill = "grey90")+
  geom_point(data = NEON_isotopes, aes(date, s1, fill = s1Binned), color = 'black', pch = 21, size = 2.5)+
  scale_fill_manual(values = c("#fde725", "#86d549", "#2ab07f", "#25858e", "#38588c", "#440154"), 
                    labels = c("36 to 40", "40 to 44", "44 to 48", "48 to 52", "52 to 56", "56 to 60"))+
  stat_summary(data = AIMS_isotopes, aes(date, s1),
               color = "black", 
               fill = "#E7315D",
               pch = 21,
               fun = mean,
               geom = "pointrange",
               fun.max = function(x) mean(x) + qt(.975, df = length(x)) * sd(x) / sqrt(length(x)),
               fun.min = function(x) mean(x) - qt(.975, df = length(x)) * sd(x) / sqrt(length(x)))+
  scale_y_continuous(limits = c(30,60))+
  labs(x = "", y = "Proportion of streamwater less than 3 months in age (%)")+
  labs(x = "", y = expression(F["yw"]* " (%)"))+
  #labs(fill = "Streamwater\n< 3 months\nin age (%)")+
  labs(fill = expression(F["yw"]* " (%)"))+
  guides(fill = guide_legend(reverse=T))+
  scale_x_date(breaks=date_breaks("months"), labels=date_format("%b %y"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(path = "figures/", "Fig7.png", dpi=300, width = 190, height = 110, units = "mm")
ggsave(path = "figures/", "Figure7.pdf", dpi=300, width = 190, height = 110, units = "mm")