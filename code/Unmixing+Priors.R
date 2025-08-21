# Sensitivity analysis of priors
library(isoWater)
library(tidyverse)

# NEON δ18O timeseries
stream_isotopes <- read_csv("data/NEON/NEON_isotopes.csv") %>%
  filter(type == "Stream") %>%
  mutate(date = mdy(date)) %>%
  filter(date >= "2019-12-01" & date <= "2021-12-01") %>%
  select(date, d18OWater, d2HWater) %>%
  # if multiple samples were collected on a single day, average to a single isotopic signal
  group_by(date) %>%
  summarise(d18OWater = mean(d18OWater), d2HWater = mean(d2HWater)) %>%
  ungroup() %>%
  mutate_if(is.numeric, round, digits = 2)

# NEON precip δ18O timeseries
precip_isotopes <- read_csv("data/NEON/NEON_isotopes.csv") %>%
  filter(type == "Precipitation") %>%
  mutate(date = mdy(date)) %>%
  filter(date >= "2018-11-01" & date <= "2021-12-31") %>% 
  select(date, d18OWater, d2HWater) %>%
  mutate_if(is.numeric, round, digits = 2) 

# precip weights
precip <- read_csv("data/LTER/Konza_Precip.csv")

# add cumulative precip between sampling dates
rows_match <- match(precip_isotopes$date, precip$date)
precip_isotopes$precip_sum_mm <- NA
for (i in 2:length(rows_match)){
  row_start_sum <- rows_match[i-1] + 1
  row_end_sum <- rows_match[i]
  precip_isotopes$precip_sum_mm[i] <- sum(precip$precip_mm[row_start_sum:row_end_sum])
}
# remove first precip isotope missing cumulative precip
precip_isotopes <- na.omit(precip_isotopes)

# same EL slope & residual standard error (estimated previously)
slope = c(6.00, 0.55)

priors <- list(
  c(0.2, 0.8), # preference towards older water
  c(0.8, 0.2), # preference towards younger water,
  c(0.153,  0.847) # prior informed by sinusoidal models
)

# Temporal unmixing to estimate young water fractions at catchment outlet
set.seed(1)
SA <- vector("list", length(priors))
for (p in seq_along(priors)) {
  
  unmixed <- stream_isotopes
  
  unmixed$prior <- priors[p]
  unmixed$s1 <- NA
  unmixed$s2 <- NA
  unmixed$s1_rhat <- NA
  unmixed$s2_rhat <- NA
  unmixed$s1_n.eff <- NA
  unmixed$s2_n.eff <- NA
  
  # Bayesian unmixing
  system.time({
    for (i in seq_along(stream_isotopes$d18OWater)) {
      
      # precip isotope signal before stream isotope sample date
      precip_signal <- precip_isotopes[precip_isotopes$date <= stream_isotopes$date[i],]
      
      # two precip sources contributing to stream isotope signal 
      s1 <- precip_signal[precip_signal$date >= stream_isotopes$date[i] - 90,]  
      s2 <- precip_signal[precip_signal$date < stream_isotopes$date[i] - 90,] 
      
      # if there are not at least 3 recent precip samples, then do not estimate source contributions
      if(nrow(s1) < 3){
        unmixed$s1[i] <- NA
        unmixed$s2[i] <- NA
        unmixed$s1_rhat[i] <- NA
        unmixed$s2_rhat[i] <- NA
        unmixed$s1_n.eff[i] <- NA
        unmixed$s2_n.eff[i] <- NA
        
        # estimate source contributions using two precip sources
      } else {
        
        # parameters defining precip distribution that fell less than 3 months ago
        s1_H <- sum(s1$d2HWater*s1$precip_sum_mm)/sum(s1$precip_sum_mm)
        s1_O <- sum(s1$d18OWater*s1$precip_sum_mm)/sum(s1$precip_sum_mm)
        s1_Hsd <- sd(s1$d2HWater)
        s1_Osd <- sd(s1$d18OWater)
        s1_HOc <- cov(s1$d2HWater, s1$d18OWater)
        
        # parameters defining precip distribution that fell longer than 3 months ago
        s2_H <- sum(s2$d2HWater*s2$precip_sum_mm)/sum(s2$precip_sum_mm)
        s2_O <- sum(s2$d18OWater*s2$precip_sum_mm)/sum(s2$precip_sum_mm)
        s2_Hsd <- sd(s2$d2HWater)
        s2_Osd <- sd(s2$d18OWater)
        s2_HOc <- cov(s2$d2HWater, s2$d18OWater)
        
        sources <- data.frame("H" = as.numeric(c(s1_H,s2_H)),
                              "O" = as.numeric(c(s1_O,s2_O)),
                              "Hsd" = as.numeric(c(s1_Hsd,s2_Hsd)),
                              "Osd" = as.numeric(c(s1_Osd,s2_Osd)),
                              "HOc" = as.numeric(c(s1_HOc,s2_HOc)))
        
        sources = iso(sources$H, sources$O, sources$Hsd, sources$Osd, sources$HOc)
        
        obs <- iso(stream_isotopes$d2HWater[i], stream_isotopes$d18OWater[i], 0.5, 0.1, 0.025)
        
        mix = mixSource(obs, sources, slope, prior=priors[[p]], ngens = 2e5, shp = 1)
        
        unmixed$s1[i] <- mix[["summary"]][4,1]
        unmixed$s2[i] <- mix[["summary"]][5,1]
        unmixed$s1_rhat[i] <- mix[["summary"]][4,8]
        unmixed$s2_rhat[i] <- mix[["summary"]][5,8]
        unmixed$s1_n.eff[i] <- mix[["summary"]][4,9]
        unmixed$s2_n.eff[i] <- mix[["summary"]][5,9]
      }
    }
  })
  SA[[p]] <- unmixed
}  

SA <- bind_rows(SA) %>%
  mutate(prior = as.character(prior),
         s1 = s1 * 100,
         s2 = s2 * 100)

uninformed <- read_csv("data/NEON_isotopes_results.csv") %>%
  mutate(prior = as.character("c(0.5, 0.5)"))

out <- bind_rows(list(uninformed, SA))

write_csv(out, "data/temporal_priors.csv")




# AIMS δ18O spatial
stream_isotopes <- read_csv("data/AIMS_isotopes_coordinates.csv") %>%
  select(-c(sampleID, replicate))

# Spatial unmixing to estimate young water fractions in headwaters
set.seed(1)
SA <- vector("list", length(priors))
for (p in seq_along(priors)) {
  
  unmixed <- stream_isotopes
  
  unmixed$prior <- priors[p]
  unmixed$s1 <- NA
  unmixed$s2 <- NA
  unmixed$s1_rhat <- NA
  unmixed$s2_rhat <- NA
  unmixed$s1_n.eff <- NA
  unmixed$s2_n.eff <- NA
  
  # Bayesian unmixing
  system.time({
    for (i in seq_along(stream_isotopes$d18OWater)) {
      
      # precip isotope signal before stream isotope sample date
      precip_signal <- precip_isotopes[precip_isotopes$date <= stream_isotopes$date[i],]
      
      # two precip sources contributing to stream isotope signal 
      s1 <- precip_signal[precip_signal$date >= stream_isotopes$date[i] - 90,]  
      s2 <- precip_signal[precip_signal$date < stream_isotopes$date[i] - 90,] 
      
      # if there are not at least 3 recent precip samples, then do not estimate source contributions
      if(nrow(s1) < 3){
        unmixed$s1[i] <- NA
        unmixed$s2[i] <- NA
        unmixed$s1_rhat[i] <- NA
        unmixed$s2_rhat[i] <- NA
        unmixed$s1_n.eff[i] <- NA
        unmixed$s2_n.eff[i] <- NA
        
        # estimate source contributions using two precip sources
      } else {
        
        # parameters defining precip distribution that fell less than 3 months ago
        s1_H <- sum(s1$d2HWater*s1$precip_sum_mm)/sum(s1$precip_sum_mm)
        s1_O <- sum(s1$d18OWater*s1$precip_sum_mm)/sum(s1$precip_sum_mm)
        s1_Hsd <- sd(s1$d2HWater)
        s1_Osd <- sd(s1$d18OWater)
        s1_HOc <- cov(s1$d2HWater, s1$d18OWater)
        
        # parameters defining precip distribution that fell longer than 3 months ago
        s2_H <- sum(s2$d2HWater*s2$precip_sum_mm)/sum(s2$precip_sum_mm)
        s2_O <- sum(s2$d18OWater*s2$precip_sum_mm)/sum(s2$precip_sum_mm)
        s2_Hsd <- sd(s2$d2HWater)
        s2_Osd <- sd(s2$d18OWater)
        s2_HOc <- cov(s2$d2HWater, s2$d18OWater)
        
        sources <- data.frame("H" = as.numeric(c(s1_H,s2_H)),
                              "O" = as.numeric(c(s1_O,s2_O)),
                              "Hsd" = as.numeric(c(s1_Hsd,s2_Hsd)),
                              "Osd" = as.numeric(c(s1_Osd,s2_Osd)),
                              "HOc" = as.numeric(c(s1_HOc,s2_HOc)))
        
        sources = iso(sources$H, sources$O, sources$Hsd, sources$Osd, sources$HOc)
        
        obs <- iso(stream_isotopes$d2HWater[i], stream_isotopes$d18OWater[i], 0.5, 0.1, 0.025)
        
        mix = mixSource(obs, sources, slope, prior=priors[[p]], ngens = 2e5, shp = 1)
        
        unmixed$s1[i] <- mix[["summary"]][4,1]
        unmixed$s2[i] <- mix[["summary"]][5,1]
        unmixed$s1_rhat[i] <- mix[["summary"]][4,8]
        unmixed$s2_rhat[i] <- mix[["summary"]][5,8]
        unmixed$s1_n.eff[i] <- mix[["summary"]][4,9]
        unmixed$s2_n.eff[i] <- mix[["summary"]][5,9]
      }
    }
  })
  SA[[p]] <- unmixed
}  

SA <- bind_rows(SA) %>%
  mutate(
    prior = as.character(prior),
    s1 = s1 * 100,
    s2 = s2 * 100
  )

uninformed <- read_csv("data/AIMS_isotopes_results.csv") %>%
  mutate(prior = as.character("c(0.5, 0.5)"))

out <- bind_rows(list(uninformed, SA))

write_csv(out, "data/spatial_priors.csv")