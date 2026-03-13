# Precip conditions (totals, percentiles, etc.) for manuscript
library(tidyverse)

source(file.path("code", "theme.R"))

end <- as.Date("2021-12-31")  # data available at the time of manuscript publication

precip <- file.path("data", "LTER", "LTER_HQ_weather_station.csv") %>%
  read_csv(show_col_types = F, col_select = c(siteID, date, precip_mm)) %>%
  filter(date >= as.Date("1991-01-01") & date <= end) # 30-year precip conditions

precip_annual <- precip %>% 
  mutate(
    year = year(date)
  ) %>% 
  group_by(year) %>% 
  summarise(
    precip_mm = sum(precip_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(
    percentile = ecdf(precip_mm)(precip_mm) * 100
  )

precip_monthly <- precip %>% 
  mutate(
    year = year(date),
    month = month(date)
  ) %>% 
  group_by(year, month) %>% 
  summarise(
    precip_mm = sum(precip_mm, na.rm=T),
    .groups = "drop"
  ) %>% 
  group_by(month) %>% 
  mutate(
    percentile = ecdf(precip_mm)(precip_mm)*100
  ) %>%
  ungroup()

precip_seasonal <- precip %>%
  mutate(
    year = year(date),
    month = month(date),
    season = case_when(
      month %in% c(3,4,5) ~ "spring",
      month %in% c(6,7,8) ~ "summer",
      TRUE ~ NA
    )
  ) %>%
  drop_na(season) %>%
  group_by(year, season) %>%
  summarise(
    precip_mm = sum(precip_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  group_by(season) %>% 
  mutate(
    percentile = ecdf(precip_mm)(precip_mm) * 100
  ) %>%
  ungroup()

precip_average <- mean(precip_annual$precip_mm)
precip_synoptic <- precip_annual$precip_mm[precip_annual$year == 2021]
precip_percentiles_per_month <- precip_monthly$percentile[precip_monthly$year == 2021]
precip_percentiles_per_season <- precip_seasonal$percentile[precip_seasonal$year == 2021]

print(paste0("Annual average precip: ", round(precip_average), " mm"))  # average 835 mm/yr reported for Konza Prairie (Hayden, 1998)
print(paste0("Calendar 2021 precip: ", round(precip_synoptic), " mm"))
print(paste0("March percentile: ", round(precip_percentiles_per_month[3]), "%"))
print(paste0("Spring (March–May) percentile: ", round(precip_percentiles_per_season[1]), "%"))
print(paste0("Summer (June–August) percentile: ", round(precip_percentiles_per_season[2]), "%"))