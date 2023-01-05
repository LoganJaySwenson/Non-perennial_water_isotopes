#Data + Clean- AIMS
library(lubridate)
library(tidyverse)

#Read: AIMS isotopes 
AIMS_isotopes <- read.csv("data/AIMS/AIMS_isotopes.csv")
AIMS_isotopes <- as_tibble(AIMS_isotopes) %>%
  separate(sampleID, c(NA, "siteID", "replicate", NA, "date"), remove = F, extra = "drop", convert = T) %>%
  mutate(date = ymd(date)) %>%
  mutate(yday = yday(date)) %>%
  mutate(month = month(date, label = T)) %>%
  mutate(watershed = case_when(
    startsWith(siteID, "SFM") ~ "SFKC",
    startsWith(siteID, "SFT") ~ "SFKC",
    startsWith(siteID, "01M") ~ "N01B",
    startsWith(siteID, "02M") ~ "N02B",
    startsWith(siteID, "20M") ~ "N20B",
    startsWith(siteID, "04M") ~ "N04D",
    startsWith(siteID, "04W") ~ "N04D",
    startsWith(siteID, "04T") ~ "N04D",
    startsWith(siteID, "N04D") ~ "N04D")) %>%
  mutate(dexcess = d2HWater - (8 * d18OWater)) %>% #Dansgaard (1964): d-excess = δ2H - 8 * δ18O
  mutate(lcexcess = d2HWater - (7.9309 * d18OWater) - 10.2770) %>% #Landwehr & Coplen (2006): line-conditioned excess = δ2H - m * δ18O - b.
  select(sampleID, siteID, watershed, month, date, d18OWater, d2HWater, yday, replicate, dexcess, lcexcess) %>%
  filter(!grepl(c("SPRING|STILLINGWELL"), sampleID)) %>% # remove spring & stilling well samples
  filter(siteID != "N04D") %>% #remove sample missing a siteID
  group_by(siteID, month) %>% # average all reps!
  mutate(d18OWater = mean(d18OWater), d2HWater = mean(d2HWater)) %>%
  slice(1) %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  ungroup()

#Read: Isotope coordinates
pts <- as_tibble(read.csv("data/AIMS/AIMS_STIC_details.csv"))

#Bind: Isotopes with coordinates
AIMS_isotopes <- left_join(AIMS_isotopes, pts, by = "siteID")

#Read: Isotopes notes
notes <- as_tibble(read.csv("data/AIMS/AIMS_isotopes_notes.csv"))
notes <- notes %>%
  mutate(date = mdy(date)) %>%
  mutate(month = month(date, label = T)) %>%
  unite(col = "id", c("month", "site"), remove = F) %>%
  select(id, flowing, ysiWaterTemp, ysiSPC)

#Bind: Isotopes with notes
AIMS_isotopes <- unite(AIMS_isotopes, col = "id", c("month", "siteID"), remove = F)
AIMS_isotopes_notes <- left_join(AIMS_isotopes, notes, by = "id")
AIMS_isotopes_notes$id <- NULL

#Save to csv for spatial analysis
write.csv(AIMS_isotopes, "data/AIMS_isotopes_coordinates.csv", row.names = F)

#Save to csv for RF model
write.csv(AIMS_isotopes_notes, "data/AIMS_isotopes_RF.csv", row.names = F)