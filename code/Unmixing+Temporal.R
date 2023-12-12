# Bayesian model to estimate young water fractions at the outlet!
library(lubridate)
library(isoWater)
library(tidyverse)

set.seed(1)

# NEON δ18O timeseries
stream_isotopes <- read_csv("data/NEON/NEON_isotopes.csv") %>%
  filter(type == "Stream") %>%
  mutate(date = mdy(date)) %>%
  filter(date >= "2019-12-01" & date <= "2021-12-01") %>%
  select(date, d18OWater, d2HWater) %>%
  group_by(date) %>%
  mutate(d18OWater = mean(d18OWater), d2HWater = mean(d2HWater)) %>%
  slice(1) %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  ungroup()

# add columns to be filled with model results
stream_isotopes$s1 <- NA
stream_isotopes$s2 <- NA
stream_isotopes$s1_rhat <- NA
stream_isotopes$s2_rhat <- NA
stream_isotopes$s1_n.eff <- NA
stream_isotopes$s2_n.eff <- NA

# same EL slope & residual standard error
slope = c(6.00, 0.55)

# NEON precip δ18O timeseries
precip_isotopes <- read_csv("data/NEON/NEON_isotopes.csv") %>%
  mutate(date = mdy(date)) %>%
  filter(type == "Precipitation") %>%
  filter(date >= "2018-11-01" & date <= "2021-12-31") %>% 
  mutate_if(is.numeric, round, digits = 2) %>%
  select(date, d18OWater, d2HWater)

# Precip weights
precip <- as_tibble(read.csv("data/LTER/Konza_Precip.csv"))
precip$date <- as_date(precip$date)

# add cumulative precipitation between sampling dates
rows_match <- match(precip_isotopes$date, precip$date)
precip_isotopes$precip_sum_mm <- NA
for (i in 2:length(rows_match)){
  row_start_sum <- rows_match[i-1] + 1
  row_end_sum <- rows_match[i]
  precip_isotopes$precip_sum_mm[i] <- sum(precip$precip_mm[row_start_sum:row_end_sum])
}
precip_isotopes <- na.omit(precip_isotopes)

# Bayesian unmixing!
system.time(for (i in 1:length(stream_isotopes$d18OWater)){
  
  temp <- precip_isotopes[precip_isotopes$date <= stream_isotopes$date[i],]
  
  s1 <- temp[temp$date >= stream_isotopes$date[i] - 90,]
  
  if(nrow(s1) < 3){
    stream_isotopes$s1[i] <- NA
    stream_isotopes$s2[i] <- NA
    stream_isotopes$s1_rhat[i] <- NA
    stream_isotopes$s2_rhat[i] <- NA
    stream_isotopes$s1_n.eff[i] <- NA
    stream_isotopes$s2_n.eff[i] <- NA
    
  } else {
    s1_H <- sum(s1$d2HWater*s1$precip_sum_mm)/sum(s1$precip_sum_mm)
    s1_O <- sum(s1$d18OWater*s1$precip_sum_mm)/sum(s1$precip_sum_mm)
    s1_Hsd <- sd(s1$d2HWater)
    s1_Osd <- sd(s1$d18OWater)
    s1_HOc <- cov(s1$d2HWater, s1$d18OWater)
    
    s2 <- temp[temp$date < stream_isotopes$date[i] - 90,]
    
    s2_H <- sum(s2$d2HWater*s2$precip_sum_mm)/sum(s2$precip_sum_mm)
    s2_O <- sum(s2$d18OWater*s2$precip_sum_mm)/sum(s2$precip_sum_mm)
    s2_Hsd <- sd(s2$d2HWater)
    s2_Osd <- sd(s2$d18OWater)
    s2_HOc <- cov(s2$d2HWater, s2$d18OWater)
    
    temp <- data.frame("H" = as.numeric(c(s1_H,s2_H)),
                       "O" = as.numeric(c(s1_O,s2_O)),
                       "Hsd" = as.numeric(c(s1_Hsd,s2_Hsd)),
                       "Osd" = as.numeric(c(s1_Osd,s2_Osd)),
                       "HOc" = as.numeric(c(s1_HOc,s2_HOc)))
    
    sources = iso(temp$H, temp$O, temp$Hsd, temp$Osd, temp$HOc)
    
    obs <- iso(stream_isotopes$d2HWater[i], stream_isotopes$d18OWater[i], 0.5, 0.1, 0.026)
    
    mix = mixSource(obs, sources, slope, ngens = 2e5, shp = 1)
    
    stream_isotopes$s1[i] <- mix[["summary"]][4,1]
    stream_isotopes$s2[i] <- mix[["summary"]][5,1]
    stream_isotopes$s1_rhat[i] <- mix[["summary"]][4,8]
    stream_isotopes$s2_rhat[i] <- mix[["summary"]][5,8]
    stream_isotopes$s1_n.eff[i] <- mix[["summary"]][4,9]
    stream_isotopes$s2_n.eff[i] <- mix[["summary"]][5,9]
    
  }
})

stream_isotopes <- stream_isotopes %>%
  mutate(s1 = s1 * 100) %>%
  mutate(s2 = s2 * 100)

# Save to results to csv
write_csv(stream_isotopes, "data/NEON_isotopes_results.csv")