#Fig S7. FYW Base model fit!
library(lubridate)
library(scales)
library(dataRetrieval)
library(tidyverse)

#Publication theme
source("code/Theme+Settings.R")

#Iteratively reweighted least squares (IRLS), Kirchner, ETH Zurich
source("code/IRLS.R")

#P δ18O timeseries
P <- as_tibble(read.csv("data/NEON/NEON_isotopes.csv"))
P <- P %>%
  filter(type == "Precipitation") %>%
  mutate(date = mdy(date)) %>%
  filter(date >= as.Date("2018-11-01") & date <= as.Date("2021-12-31")) %>%
  group_by(date) %>%
  mutate(d18OWater = mean(d18OWater), d2HWater = mean(d2HWater)) %>%
  slice(1) %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  ungroup() %>%
  mutate(t = as.numeric(date - ymd("2018-11-01"))) %>%
  select(date, t, d18OWater, d2HWater) %>%
  mutate(mean18O = mean(d18OWater)) %>%
  mutate(phase = (2*pi/365)*t) %>%
  mutate(cosct = cos(phase)) %>%
  mutate(sinct = sin(phase)) 

#P weights
Precip <- as_tibble(read.csv("data/LTER/Konza_Precip.csv"))
Precip$date <- as_date(Precip$date)

#add cumulative precipitation between sampling dates
rows_match <- match(P$date, Precip$date)
P$precip_sum_mm <- NA
for (i in 2:length(rows_match)){
  row_start_sum <- rows_match[i-1] + 1
  row_end_sum <- rows_match[i]
  P$precip_sum_mm[i] <- sum(Precip$precip_mm[row_start_sum:row_end_sum])
}
P <- na.omit(P)

#Fit base model
P_base_model <- IRLS(P$d18OWater, as.matrix(P[,c("cosct","sinct")]), ww=P$precip_sum_mm, type="Cauchy")
summary(P_base_model$fit)
aP <- as.numeric(P_base_model$fit$coefficients[2])
bP <- as.numeric(P_base_model$fit$coefficients[3])
P_model_fit <- cbind(P, predicted = fitted(P_base_model$fit))


#Q δ18O timeseries
Q <- as_tibble(read.csv("data/NEON/NEON_isotopes.csv"))
Q <- Q %>%
  filter(type == "Stream") %>%
  mutate(date = mdy(date)) %>%
  filter(date >= "2020-10-01" & date <= "2021-09-30") %>% #select 2021 water year
  group_by(date) %>%
  mutate(d18OWater = mean(d18OWater), d2HWater = mean(d2HWater)) %>%
  slice(1) %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  ungroup() %>%
  mutate(t = as.numeric(date - ymd("2018-11-01"))) %>%
  select(date, t, d18OWater, d2HWater) %>%
  mutate(mean18O = mean(d18OWater)) %>%
  mutate(phase = (2*pi/365)*t) %>%
  mutate(cosct = cos(phase)) %>%
  mutate(sinct = sin(phase)) 

#Q weights 
Flow <- as_tibble(readNWISdv(siteNumber = "06879650", parameterCd = "00060", startDate = "", endDate = ""))
Flow <- renameNWISColumns(Flow)
Flow <- Flow %>%
  mutate(Flow = (Flow * 0.028316847)) %>% #cfs to metric
  select(Date, Flow) %>%
  rename(date = Date, flow = Flow) 

#Q δ18O timeseries with weights
Qw <- left_join(Q, Flow, by = "date")

#Fit base model
Q_base_model <- IRLS(Qw$d18OWater, as.matrix(Qw[,c("cosct","sinct")]), ww=Qw$flow, type="Cauchy")
Q_resids_wt <- Q_base_model$wt
summary(Q_base_model$fit)
aQ <- as.numeric(Q_base_model$fit$coefficients[2])
bQ <- as.numeric(Q_base_model$fit$coefficients[3])
Qw_model_fit <- cbind(Qw, predicted = fitted(Q_base_model$fit))


#FYW
AP <- sqrt(aP^2 + bP^2)
AQ <- sqrt(aQ^2 + bQ^2)

FYW = (AQ/AP) * 100
FYW

P_model_fit <- P_model_fit %>%
  mutate(type = "P")

Qw_model_fit <- Qw_model_fit %>%
  mutate(type = "Q")

FYW_model <- bind_rows(P_model_fit, Qw_model_fit) %>%
  remove_rownames(.)

#Plot sinusoidal model fit!
ggplot()+
  geom_line(data = P_model_fit, aes(date, predicted), linetype = 'dashed', color = "#000000")+
  geom_line(data = Qw_model_fit, aes(date, predicted), linetype = 'dashed', color = "#377EB8")+
  geom_point(data = FYW_model, aes(date, d18OWater, fill = type), pch = 21)+
  annotate("text", x = as.Date("2019-04-01"), y = 4, label = paste("Mass ww method:", round(FYW, digits = 2),"%"), size = 10/.pt)+
  scale_fill_manual(values = c("#000000", "#377EB8"), labels = c(expression(δ^18*O[P]~(t)), expression(δ^18*O[Q]~(t))))+
  labs(fill = "")+
  labs(x = "", y = "\U03B4\U00B9\U2078O (‰)")+
  scale_x_date(breaks=date_breaks("years"), labels=date_format("%Y")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text=element_text(size=10),
        legend.position = c(0.922, 0.88),
        legend.background = element_rect(fill = "white", color = "black"))
ggsave(path = "figures/", "FigS7.png", dpi=300, width = 190, height = 110, units = "mm")