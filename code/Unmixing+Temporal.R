#Bayesian Unmixing at the outlet!
library(lubridate)
library(isoWater)
library(tidyverse)

set.seed(1)

#Q δ18O timeseries
Q <- read.csv("data/NEON/NEON_isotopes.csv")
Q <- as_tibble(Q) %>%
  filter(type == "Stream") %>%
  mutate(date = mdy(date)) %>%
  filter(date >= "2019-12-01" & date <= "2021-12-01") %>%
  select(date, d18OWater, d2HWater) %>%
  group_by(date) %>%
  mutate(d18OWater = mean(d18OWater), d2HWater = mean(d2HWater)) %>%
  slice(1) %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  ungroup()

#add columns to be filled with model results
Q$s1 <- NA
Q$s2 <- NA
Q$s1_rhat <- NA
Q$s2_rhat <- NA
Q$s1_n.eff <- NA
Q$s2_n.eff <- NA

#P δ18O timeseries
P <- read.csv("data/NEON/NEON_isotopes.csv")
P <- as_tibble(P) %>%
  mutate(date = mdy(date)) %>%
  filter(type == "Precipitation") %>%
  filter(date >= "2018-11-01" & date <= "2021-12-31") %>% 
  mutate_if(is.numeric, round, digits = 2) %>%
  select(date, d18OWater, d2HWater)

#P weights
Precip <- read.csv("data/LTER/Konza_Precip.csv")
Precip <- as_tibble(Precip) %>%
  mutate(date = ymd(date))

#same EL slope & residual standard error
slope = c(6.00, 0.55)

#add cumulative precipitation between sampling dates
rows_match <- match(P$date, Precip$date)
P$precip_sum_mm <- NA
for (i in 2:length(rows_match)){
  row_start_sum <- rows_match[i-1] + 1
  row_end_sum <- rows_match[i]
  P$precip_sum_mm[i] <- sum(Precip$precip_mm[row_start_sum:row_end_sum])
}
P <- na.omit(P)

#Run Bayesian unmixing!
system.time(for (i in 1:length(Q$d18OWater)){
  
  temp <- P[P$date <= Q$date[i],]
  
  s1 <- temp[temp$date >= Q$date[i] - 90,]
  
  if(nrow(s1) < 3){
    Q$s1[i] <- NA
    Q$s2[i] <- NA
    Q$s1_rhat[i] <- NA
    Q$s2_rhat[i] <- NA
    Q$s1_n.eff[i] <- NA
    Q$s2_n.eff[i] <- NA
    
  } else {
    s1_H <- sum(s1$d2HWater*s1$precip_sum_mm)/sum(s1$precip_sum_mm)
    s1_O <- sum(s1$d18OWater*s1$precip_sum_mm)/sum(s1$precip_sum_mm)
    s1_Hsd <- sd(s1$d2HWater)
    s1_Osd <- sd(s1$d18OWater)
    s1_HOc <- cov(s1$d2HWater, s1$d18OWater)
    
    s2 <- temp[temp$date < Q$date[i] - 90,]
    
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
    
    obs <- iso(Q$d2HWater[i], Q$d18OWater[i], 0.5, 0.1, 0.026)
    
    mix = mixSource(obs, sources, slope, ngens = 2e5, shp = 1)
    
    Q$s1[i] <- mix[["summary"]][4,1]
    Q$s2[i] <- mix[["summary"]][5,1]
    Q$s1_rhat[i] <- mix[["summary"]][4,8]
    Q$s2_rhat[i] <- mix[["summary"]][5,8]
    Q$s1_n.eff[i] <- mix[["summary"]][4,9]
    Q$s2_n.eff[i] <- mix[["summary"]][5,9]
    
  }
})

Q <- Q %>%
  mutate(s1 = s1 * 100) %>%
  mutate(s2 = s2 * 100)

#Save to results to csv
write.csv(Q, "data/Q_isotopes_results.csv", row.names = F)