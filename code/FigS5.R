#Fig S5. δ18O & δ2H isotope biplot
library(lubridate)
library(tidyverse)

#Publication theme
source("code/Theme+Settings.R")

#Read: P δ18O timeseries
P <- read.csv("data/NEON/NEON_isotopes.csv")
P <- as_tibble(P) %>%
  mutate(date = mdy(date)) %>%
  filter(type == "Precipitation") %>%
  filter(date >= "2018-11-01" & date <= "2021-12-31") %>% 
  mutate_if(is.numeric, round, digits = 2) %>%
  select(date, d18OWater, d2HWater)

#Read: AIMS isotopes
AIMS_isotopes <- as_tibble(read.csv("data/AIMS_isotopes_RF.csv"))

#MWL slope & intercept
MWL <- print(summary(lm(d2HWater ~ d18OWater, data = P)))
MWL.intercept <- MWL$coefficients[1]
MWL.slope <- MWL$coefficients[2]

# EL slope & residual standard error
EL <- print(summary(lm(d2HWater ~ d18OWater, data = AIMS_isotopes)))
EL.intercept <- EL$coefficients[1]
EL.slope <- EL$coefficients[2]

#Plot!
AIMS_isotopes$month <- factor(AIMS_isotopes$month, levels = c( "Jun", "Jul", "Aug"))
ggplot()+
  geom_abline(intercept = MWL.intercept, slope = MWL.slope)+
  annotate("text", x = -3, y = -9.0, label = "MWL", size = 9/.pt)+
  annotate("text", x = -3, y = -17.5, label = "EL", size = 9/.pt)+
  geom_abline(intercept = EL.intercept, slope = EL.slope)+
  geom_point(data = AIMS_isotopes, aes(x = d18OWater, y = d2HWater, fill = month), color = "black", pch = 21, size = 2)+
  scale_fill_manual(values = c('#0082c8', '#3cb44b', '#e6194b'), labels = c("June", "July", "August"))+
  labs(x = "\U03B4\U00B9\U2078O (‰)", y = "\U03B4\U00B2H (‰)")+
  labs(fill = "")
ggsave(path = "figures/", "FigS5.png", dpi=300, width = 90, height = 90, units = "mm")