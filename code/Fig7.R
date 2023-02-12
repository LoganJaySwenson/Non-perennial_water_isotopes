#Fig 7. Temporal- % of streamwater less than ~3 months in age. 
#streamwater age is defined as the mean of the posterior distribution of source mixtures.
library(lubridate)
library(tidyverse)
library(viridis)
library(scales)

#Publication theme
source("code/Theme+Settings.R")

#Q δ18O timeseries
Q <- read.csv("data/Q_isotopes_results.csv")
Q <- as_tibble(Q) %>%
  filter(date >= "2020-10-01" & date <= "2021-09-30") %>% #filter 2021 WY
  mutate(date = ymd(date)) %>%
  na.omit()

#Read: AIMS isotopes
AIMS_isotopes <- as_tibble(read.csv("data/AIMS_isotopes_results.csv"))
AIMS_isotopes <- AIMS_isotopes %>%
  mutate(date = ymd(date)) %>%
  mutate(date = replace(date, date == "2021-06-06", "2021-06-07"))

#Plot!
Q$s1Binned <- cut(Q$s1, breaks = c(seq(from = 36, to = 64, by = 4)))
AIMS_isotopes$s1Binned <- cut(AIMS_isotopes$s1, breaks = c(seq(from = 36, to = 64, by = 4)))
ggplot()+
  geom_smooth(data = Q, aes(date, s1,))+
  geom_hline(yintercept = 54.35, linetype = "dashed", color = "#636363")+
  annotate("text", x = as.Date("2020-10-10"), y = 56.8, label = "FYW: 54.96%", size = 8/.pt, color = "#636363")+
  geom_point(data = Q, aes(date, s1, fill = s1Binned), color = 'black', pch = 21, size = 2.5)+
  scale_fill_viridis(direction = -1, discrete = T)+
  stat_summary(data = AIMS_isotopes, aes(date, s1),
               color = "black", 
               fill = "#E7315D",
               pch = 21,
               fun = mean,
               geom = "pointrange",
               fun.max = function(x) mean(x) + qt(.975, df = length(x)) * sd(x) / sqrt(length(x)),
               fun.min = function(x) mean(x) - qt(.975, df = length(x)) * sd(x) / sqrt(length(x)))+
  scale_y_continuous(limits = c(0, 80))+
  labs(x = "", y = "Proportion of streamwater less than 3 months in age (%)")+
  labs(fill = "Streamwater\n< 3 months\nin age (%)")+
  guides(fill = guide_legend(reverse=T))+
  scale_x_date(breaks=date_breaks("months"), labels=date_format("%b %y"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(path = "figures/", "Fig7.png", dpi=300, width = 190, height = 110, units = "mm")