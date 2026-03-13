# Spatial and temporal unmixing to estimate young water fractions in headwaters and at catchment outlet
library(isoWater)
library(tidyverse)

source(file.path("code", "theme.R"))

set.seed(1)

min_q <- 0.001
end <- as.Date("2021-12-31") # data available at the time of manuscript publication

# NEON isotopes
NEON_isotopes <- file.path("data", "NEON", "NEON_isotopes.csv") %>% 
  read_csv(show_col_types = F) 

# NEON precip δ18O 
NEON_precip_isotopes <- NEON_isotopes%>%
  filter(type == "PPT", date <= end) %>% 
  mutate(
    t = as.numeric(date - date[2]), 
  )

# NEON stream δ18O at outlet for 2021 water year
NEON_stream_isotopes <- NEON_isotopes %>%
  filter(type == "SW", date >= as.Date("2020-10-01") & date <= as.Date("2021-09-30")) %>% 
  mutate(
    t = as.numeric(date - date[1]),
  )

# Precip at HQ weather station
precip <- file.path("data", "LTER", "LTER_HQ_weather_station.csv") %>%
  read_csv(show_col_types = F, col_select = c(siteID, date, precip_mm)) %>%
  filter(date >= min(NEON_precip_isotopes$date) & date <= end) 

# precip weights
rows_match <- match(NEON_precip_isotopes$date, precip$date)
NEON_precip_isotopes$precip_sum_mm <- NA
for (i in 2:length(rows_match)){
  row_start_sum <- rows_match[i-1] + 1
  row_end_sum <- rows_match[i]
  NEON_precip_isotopes$precip_sum_mm[i] <- sum(precip$precip_mm[row_start_sum:row_end_sum],na.rm=T)
}
NEON_precip_isotopes <- na.omit(NEON_precip_isotopes)

# same EL slope & residual standard error (estimated previously)
slope = c(6.00, 0.55)

priors <- list(
  c(0.5, 0.5),       # uniformed prior
  c(0.038, 0.962),   # informed by time-weighted young water fraction 
  c(0.145,  0.855)   # informed by flow-weighted young water fraction
)

out <- vector("list", length(priors))
for (p in seq_along(priors)) {
  
  stream_isotopes <- NEON_stream_isotopes
  precip_isotopes <- NEON_precip_isotopes
  
  stream_isotopes$prior <- priors[p]
  stream_isotopes$s1 <- NA
  stream_isotopes$s2 <- NA
  stream_isotopes$s1_rhat <- NA
  stream_isotopes$s2_rhat <- NA
  stream_isotopes$s1_neff <- NA
  stream_isotopes$s2_neff <- NA
  
  # Bayesian unmixing
  system.time({
    for (i in seq_along(stream_isotopes$d18O)) {
      
      # precip isotope signal before stream isotope sample date
      precip_signal <- precip_isotopes[precip_isotopes$date <= stream_isotopes$date[i],]
      
      # two precip sources contributing to stream isotope signal 
      s1 <- precip_signal[precip_signal$date >= stream_isotopes$date[i] - 90,]  
      s2 <- precip_signal[precip_signal$date < stream_isotopes$date[i] - 90,] 
      
      # if there are not at least 3 recent precip samples, then do not estimate source contributions
      if(nrow(s1) < 3){
        next
        
        # estimate source contributions using two precip sources
      } else {
        
        # parameters defining precip distribution that fell less than 3 months ago
        s1_H <- sum(s1$d2H*s1$precip_sum_mm)/sum(s1$precip_sum_mm)
        s1_O <- sum(s1$d18O*s1$precip_sum_mm)/sum(s1$precip_sum_mm)
        s1_Hsd <- sd(s1$d2H)
        s1_Osd <- sd(s1$d18O)
        s1_HOc <- cov(s1$d2H, s1$d18O)
        
        # parameters defining precip distribution that fell longer than 3 months ago
        s2_H <- sum(s2$d2H*s2$precip_sum_mm)/sum(s2$precip_sum_mm)
        s2_O <- sum(s2$d18O*s2$precip_sum_mm)/sum(s2$precip_sum_mm)
        s2_Hsd <- sd(s2$d2H)
        s2_Osd <- sd(s2$d18O)
        s2_HOc <- cov(s2$d2H, s2$d18O)
        
        sources <- data.frame("H" = as.numeric(c(s1_H,s2_H)),
                              "O" = as.numeric(c(s1_O,s2_O)),
                              "Hsd" = as.numeric(c(s1_Hsd,s2_Hsd)),
                              "Osd" = as.numeric(c(s1_Osd,s2_Osd)),
                              "HOc" = as.numeric(c(s1_HOc,s2_HOc)))
        
        sources = iso(sources$H, sources$O, sources$Hsd, sources$Osd, sources$HOc)
        
        obs <- iso(stream_isotopes$d2H[i], stream_isotopes$d18O[i], 0.5, 0.1, 0.025)
        
        mix = mixSource(obs, sources, slope, prior=priors[[p]], ngens = 1e6, shp = 1, ncores=3)
  
        stream_isotopes$s1[i] <- mix[["summary"]][4,1]
        stream_isotopes$s2[i] <- mix[["summary"]][5,1]
        stream_isotopes$s1_rhat[i] <- mix[["summary"]][4,8]
        stream_isotopes$s2_rhat[i] <- mix[["summary"]][5,8]
        stream_isotopes$s1_neff[i] <- mix[["summary"]][4,9]
        stream_isotopes$s2_neff[i] <- mix[["summary"]][5,9]
      }
    }
  })
  out[[p]] <- stream_isotopes
}  
out <- bind_rows(out)
out$prior <- as.character(out$prior)
write_csv(out, "data/temporal_priors.csv")



# AIMS δ18O in headwaters
AIMS_isotopes <- file.path("data", "AIMS", "AIMS_isotopes.csv") %>% 
  read_csv(show_col_types = F)

out <- vector("list", length(priors))
for (p in seq_along(priors)) {
  
  stream_isotopes <- AIMS_isotopes
  precip_isotopes <- NEON_precip_isotopes
  
  stream_isotopes$prior <- priors[p]
  stream_isotopes$s1 <- NA
  stream_isotopes$s2 <- NA
  stream_isotopes$s1_rhat <- NA
  stream_isotopes$s2_rhat <- NA
  stream_isotopes$s1_neff <- NA
  stream_isotopes$s2_neff <- NA
  
  # Bayesian unmixing
  system.time({
    for (i in seq_along(stream_isotopes$d18O)) {
      
      # precip isotope signal before stream isotope sample date
      precip_signal <- precip_isotopes[precip_isotopes$date <= stream_isotopes$date[i],]
      
      # two precip sources contributing to stream isotope signal 
      s1 <- precip_signal[precip_signal$date >= stream_isotopes$date[i] - 90,]  
      s2 <- precip_signal[precip_signal$date < stream_isotopes$date[i] - 90,] 
      
      # if there are not at least 3 recent precip samples, then do not estimate source contributions
      if(nrow(s1) < 3){
        next
        
        # estimate source contributions using two precip sources
      } else {
        
        # parameters defining precip distribution that fell less than 3 months ago
        s1_H <- sum(s1$d2H*s1$precip_sum_mm)/sum(s1$precip_sum_mm)
        s1_O <- sum(s1$d18O*s1$precip_sum_mm)/sum(s1$precip_sum_mm)
        s1_Hsd <- sd(s1$d2H)
        s1_Osd <- sd(s1$d18O)
        s1_HOc <- cov(s1$d2H, s1$d18O)
        
        # parameters defining precip distribution that fell longer than 3 months ago
        s2_H <- sum(s2$d2H*s2$precip_sum_mm)/sum(s2$precip_sum_mm)
        s2_O <- sum(s2$d18O*s2$precip_sum_mm)/sum(s2$precip_sum_mm)
        s2_Hsd <- sd(s2$d2H)
        s2_Osd <- sd(s2$d18O)
        s2_HOc <- cov(s2$d2H, s2$d18O)
        
        sources <- data.frame("H" = as.numeric(c(s1_H,s2_H)),
                              "O" = as.numeric(c(s1_O,s2_O)),
                              "Hsd" = as.numeric(c(s1_Hsd,s2_Hsd)),
                              "Osd" = as.numeric(c(s1_Osd,s2_Osd)),
                              "HOc" = as.numeric(c(s1_HOc,s2_HOc)))
        
        sources = iso(sources$H, sources$O, sources$Hsd, sources$Osd, sources$HOc)
        
        obs <- iso(stream_isotopes$d2H[i], stream_isotopes$d18O[i], 0.5, 0.1, 0.025)
        
        mix = mixSource(obs, sources, slope, prior=priors[[p]], ngens = 1e6, shp = 1, ncores=3)
        
        stream_isotopes$s1[i] <- mix[["summary"]][4,1]
        stream_isotopes$s2[i] <- mix[["summary"]][5,1]
        stream_isotopes$s1_rhat[i] <- mix[["summary"]][4,8]
        stream_isotopes$s2_rhat[i] <- mix[["summary"]][5,8]
        stream_isotopes$s1_neff[i] <- mix[["summary"]][4,9]
        stream_isotopes$s2_neff[i] <- mix[["summary"]][5,9]
      }
    }
  })
  out[[p]] <- stream_isotopes
}  

out <- bind_rows(out)
out$prior <- as.character(out$prior)
write_csv(out, "data/spatial_priors.csv")
