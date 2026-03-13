# Process AIMS isotopes
library(tidyverse)

AIMS_isotopes <- file.path("data", "AIMS", "AIMS_isotopes_DMP.csv") %>%
  read_csv(show_col_types = F) %>%
  separate(sampleID, c(NA, "siteID", "replicate", NA, "date"), sep="_", remove=F, convert=T) %>%
  mutate(
    date = ymd(date),
    month = month(date),
    watershed = case_when(
      str_detect(siteID, "SFM") ~ "SFKC",
      str_detect(siteID, "SFT") ~ "SFKC",
      str_detect(siteID, "01M") ~ "N01B",
      str_detect(siteID, "02M") ~ "N02B",
      str_detect(siteID, "20M") ~ "N20B",
      str_detect(siteID, "04M") ~ "N04D",
      str_detect(siteID, "04W") ~ "N04D",
      str_detect(siteID, "04T") ~ "N04D",
      str_detect(siteID, "N04D") ~ "N04D",
      TRUE ~ NA)
    ) %>%
  # remove SPRING & STILLING WELL samples
  filter(!str_detect(sampleID, "SPRING|STILLINGWELL")) %>%
  # remove sample missing SITEID
  filter(siteID != "N04D") %>% 
  # average d18O and d2H across replicates
  group_by(siteID, date, month, watershed) %>%
  summarise(
    d18O = round(mean(d18O), digits=2),
    d2H  = round(mean(d2H), digits=2),
    # Dansgaard (1964): d-excess = δ2H - 8 * δ18O
    dexcess = round(d2H - (8 * d18O), digits=2),
    # Landwehr & Coplen (2006): line-conditioned excess = δ2H - m * δ18O - b.
    lcexcess = round(d2H - (7.93 * d18O) - 10.28, digits=2),
    .groups = "drop"
  )

AIMS_sites <- file.path("data", "AIMS", "AIMS_sampling_sites.csv") %>%
  read_csv(show_col_types = F)

AIMS_sampling_notes <- file.path("data", "AIMS", "AIMS_sampling_notes.csv") %>%
  read_csv(show_col_types = F) %>%
  mutate(
    date = mdy(date),
    month = month(date),
    flowing = case_when(
      flowing == "y" ~ TRUE,
      flowing == "n" ~ FALSE,
      TRUE ~ NA)
    )

AIMS_isotopes <- AIMS_isotopes %>%
  left_join(AIMS_sites, by="siteID") %>%
  left_join(select(AIMS_sampling_notes, siteID, month, flowing), by=c("siteID", "month"))

write_csv(AIMS_isotopes, file.path("data", "AIMS", "AIMS_isotopes.csv"))


# NEON isotopes at KING and KONZ sites
save_dir <- file.path("data", "NEON")
dir.create(save_dir, showWarnings=F)

# NEON (National Ecological Observatory Network). Stable isotopes in surface water (DP1.20206.001), 
# RELEASE-2026. https://doi.org/10.48443/nd6c-h521. 
# Dataset accessed on January 28, 2026.
NEON_stream_isotopes <- neonUtilities::loadByProduct(dpID="DP1.20206.001", 
                                                     site="KING", 
                                                     package="basic", check.size=F) %$%
  asi_externalLabH2OIsotopes %>% 
  as_tibble() %>% 
  mutate(
    dpID = "DP1.20206.001",
    collectDate = ymd(collectDate),
    type = "SW"
    ) %>%
  select(
    siteID, 
    dpID,
    type,
    date = collectDate, 
    d18O = d18OWater, 
    d2H = d2HWater
    ) %>% 
  arrange(date) %>%
  group_by(siteID, date, dpID, type) %>%
  summarise(
    d18O = round(mean(d18O), digits=2),
    d2H  = round(mean(d2H), digits=2),
    .groups = "drop"
  )

# NEON (National Ecological Observatory Network). Stable isotopes in precipitation (DP1.00038.001), 
# RELEASE-2026. https://doi.org/10.48443/qm1b-6d55. 
# Dataset accessed on January 28, 2026.
NEON_precip_isotopes <- neonUtilities::loadByProduct(dpID="DP1.00038.001", 
                                                     site="KONZ", 
                                                     package="basic", check.size=F) %$%
  wdi_isoPerSample %>% 
  as_tibble() %>% 
  mutate(
    dpID = "DP1.00038.001",
    collectDate = date(collectDate),
    type = "PPT"
  ) %>%
  select(
    siteID, 
    dpID, 
    type,
    date = collectDate, 
    d18O = d18OWater, 
    d2H = d2HWater
    ) %>%
  arrange(date) %>%
  group_by(siteID, date, dpID, type) %>%
  summarise(
    d18O = round(mean(d18O), digits=2),
    d2H  = round(mean(d2H), digits=2),
    .groups = "drop"
  )

# NEON (National Ecological Observatory Network). Stable isotopes in groundwater (DP1.20276.001), 
# RELEASE-2026. https://doi.org/10.48443/cqyt-p581. 
# Dataset accessed on January 28, 2026.
NEON_groundwater_isotopes <- neonUtilities::loadByProduct(dpID="DP1.20276.001", 
                                                          site="KING", 
                                                          package="basic", check.size=F) %$%
  gsi_externalLabH2OIsotopes %>% 
  as_tibble() %>% 
  mutate(
    dpID = "DP1.20276.001",
    collectDate = ymd(collectDate),
    type = "GW"
  ) %>%
  select(
    siteID,
    dpID, 
    type,
    date = collectDate, 
    d18O = d18OWater, 
    d2H = d2HWater
    ) %>%
  arrange(date) %>%
  group_by(siteID, date, dpID, type) %>%
  summarise(
    d18O = round(mean(d18O), digits=2),
    d2H  = round(mean(d2H), digits=2),
    .groups = "drop"
  )

NEON_isotopes <- bind_rows(NEON_stream_isotopes, NEON_precip_isotopes, NEON_groundwater_isotopes)

write_csv(NEON_isotopes, file.path(save_dir, "NEON_isotopes.csv"))

# Konza Prairie HQ weather station 
# Nippert, J. 2026. AWE01 meteorological data from the konza prairie headquarters weather station 
# Environmental Data Initiative. https://doi.org/10.6073/pasta/645f6ebe7d75869d1b6747174aaa342b (Accessed 2026-01-28).
# Dataset accessed on January 28, 2026.

# air temperature at 2 m (⁰C);
# relative humidity at 2 m (%); 
# solar radiation down (0.3-3.0 μm, Wm-2); 
# wind speed at 3 m (ms-1); 
# wind direction at 3 m (degrees); 
# precipitation (mm) 
# soil temperature at 25 cm (⁰C).

save_dir <- file.path("data", "LTER")
dir.create(save_dir, showWarnings=F)

meteorology <- file.path(save_dir, "AWE012.csv") %>%
  read_csv(show_col_types = F) %>%
  rename(
    data_code = DATACODE,
    record_type = RECTYPE,
    siteID = WATERSHED,
    year = RECYEAR, 
    month = RECMONTH,
    day = RECDAY,
    doy = DAYOFYEAR,
    precip_mm = DPPT,
    max_temp_c = TMAX,
    min_temp_c= TMIN,
    average_temp_c = TAVE,
    humidity_pct = DHUMID,
    solar_radiation_wm2 = DSRAD,
    max_soil_temp_c = SMAX,
    min_soil_temp_c = SMIN,
    average_soil_temp_c = SAVE,
    wind_speed_ms = WAVE
  ) %>%
  mutate(
    date = make_date(year, month, day),
    across(where(is.character), ~ na_if(.x, ".")),
    across(c(precip_mm, max_temp_c, min_temp_c, average_temp_c, 
             humidity_pct, solar_radiation_wm2, wind_speed_ms, 
             max_soil_temp_c, min_soil_temp_c, average_soil_temp_c),
           as.numeric),
  ) %>%
  select(
    siteID, 
    date, 
    precip_mm, 
    max_temp_c, 
    min_temp_c, 
    average_temp_c, 
    humidity_pct, 
    solar_radiation_wm2, 
    wind_speed_ms,
    max_soil_temp_c, 
    min_soil_temp_c, 
    average_soil_temp_c
    )
write_csv(meteorology, file.path(save_dir, "LTER_HQ_weather_station.csv"))


# Streamflow at USGS 06879650
# Dataset accessed on January 28, 2026.
save_dir <- file.path("data", "USGS")
dir.create(save_dir, showWarnings=F)

flow <- dataRetrieval::readNWISdv(siteNumber = "06879650", 
                                  parameterCd="00060", 
                                  endDate = as.Date("2025-12-31")) %>%
  dataRetrieval::renameNWISColumns() %>%
  as_tibble() %>% 
  rename_with(tolower) %>% 
  mutate(
    flow_m3s = flow * 0.0283 # ft3/s to m3/s
    ) %>%
  select(gauge_id = site_no, date, flow_m3s, flow_cd)

write_csv(flow, file.path(save_dir, "USGS_06879650.csv"))
