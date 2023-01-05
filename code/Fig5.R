#Fig 5. Variation in δ18O with d-excess & flow state
library(lubridate)
library(patchwork)
library(tidyverse)

#Publication theme
source("code/Theme+Settings.R")

#Read: AIMS isotopes
AIMS_isotopes <- read.csv("data/AIMS_isotopes_RF.csv")
AIMS_isotopes <- as_tibble(AIMS_isotopes) %>% select(1:18)

#Plot δ18O & d-excess!
AIMS_isotopes$month <- factor(AIMS_isotopes$month, levels = c( "Jun", "Jul", "Aug"))
p1 <-
  ggplot()+
  geom_point(data = AIMS_isotopes, aes(d18OWater, dexcess, fill = month), color= 'black', pch=21, size=2)+
  scale_fill_manual(values = c('#0082c8', '#3cb44b', '#e6194b'), labels = c("June", "July", "August"))+
  labs(x = "\U03B4\U00B9\U2078O (‰)", y = "d-excess (‰)")+
  labs(fill = "", shape = "")

#Plot δ18O & flow state!
AIMS_isotopes$flowing <- factor(AIMS_isotopes$flowing, levels = c("y", "n"))
p2 <-
  ggplot(AIMS_isotopes, aes(month, d18OWater, fill = flowing))+
  geom_boxplot(width = 0.5, outlier.shape = 21, outlier.size = 2)+
  labs(x = "", y = "\U03B4\U00B9\U2078O (‰)")+
  scale_fill_manual(values = c("#619CFF", "#F8766D"), name = "", labels = c("Flowing", "Pooled"))+
  scale_x_discrete(labels = c("June", "July", "August"))

p1 + p2 + plot_annotation(tag_levels = "a")
ggsave(path = "figures/", "Fig5.png", dpi=300, width = 190, height = 90, units = "mm")