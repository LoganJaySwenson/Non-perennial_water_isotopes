# Fig 2. Konza Prairie conditions for 2021 Water Year
library(scales)
library(patchwork)
library(lubridate)
library(tidyverse)

# Publication theme
source("code/Theme+Settings.R")

date_start <- as.Date("2020-10-01")
date_end <- as.Date("2021-09-30")

# Read: Streamflow at Kings Creek (USGS 06879650), Konza Prairie, KS
Flow <- dataRetrieval::readNWISdv(siteNumber = "06879650", parameterCd = "00060", startDate = "", endDate = "") %>%
  dataRetrieval::renameNWISColumns() %>%
  rename_with(.cols = everything(), tolower) %>%
  filter(date >= date_start & date <= date_end) %>%
  mutate(flow = (flow * 0.028316847)) #cfs to metric

# Set minimum for plots
min_q <- 0.001
Flow$flow_forlog <- Flow$flow
Flow$flow_forlog[Flow$flow_forlog < min_q] <- min_q

# Plot!
p1 <-
  ggplot()+
  geom_vline(xintercept = as.Date("2021-06-07"), color = "#E41A1C", linetype = "dashed")+
  geom_vline(xintercept = as.Date("2021-07-13"), color = "#E41A1C", linetype = "dashed")+
  geom_vline(xintercept = as.Date("2021-08-09"), color = "#E41A1C", linetype = "dashed")+
  geom_line(data = Flow, aes(date, flow_forlog), color = "#377EB8")+
  labs(x = "", y = "Discharge (m\U00B3/s)")+
  labs(color = "")+
  scale_x_date(expand = c(0,0), limit=c(as.Date("2020-10-01"),as.Date("2021-09-30")),
               breaks=date_breaks("months"), labels=date_format("%b %y"))+
  scale_y_log10(
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.001, 10),
                expand = c(0,0)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Read: GW levels
Mor3_5 <- read_csv("data/LTER/GW_Levels_3-5Mor.csv", show_col_types=F) %>%
  rename_with(.cols = everything(), tolower) %>%
  mutate(date = mdy(date), 
         id = "Mor 3-5") %>%
  filter(date >= date_start & date <= date_end) %>%
  select(date, id, level)

Mor3_5_1 <- read_csv("data/LTER/GW_Levels_3-5-1Mor.csv", show_col_types=F) %>%
  rename_with(.cols = everything(), tolower) %>%
  mutate(date = mdy(date), 
         id = "Mor 3-5-1") %>%
  filter(date >= date_start & date <= date_end) %>%
  select(date, id, level)

Mor4_6 <- read_csv("data/LTER/GW_Levels_4-6Mor.csv", show_col_types=F) %>%
  rename_with(.cols = everything(), tolower) %>%
  mutate(date = mdy(date), 
         id = "Mor 4-6") %>%
  filter(date >= date_start & date <= date_end) %>%
  select(date, id, level)

LowerEis4_6 <- read_csv("data/LTER/GW_Levels_4-6Eis1.csv", show_col_types=F) %>%
  rename_with(.cols = everything(), tolower) %>%
  mutate(date = mdy(date), 
         id = "Lower Eis 4-6") %>%
  filter(date >= date_start & date <= date_end) %>%
  select(date, id, level)

UpperEis4_6 <- read_csv("data/LTER/GW_Levels_4-6Eis2.csv", show_col_types=F) %>%
  rename_with(.cols = everything(), tolower) %>%
  mutate(date = mdy(date), 
         id = "Upper Eis 4-6") %>%
  filter(date >= date_start & date <= date_end) %>%
  select(date, id, level)

# GW levels
GW_Levels <- bind_rows(Mor4_6, LowerEis4_6, UpperEis4_6)

# Plot!
GW_Levels$id <- factor(GW_Levels$id, levels = c("Upper Eis 4-6", "Lower Eis 4-6", "Mor 4-6"))
p2 <-
  ggplot()+
  geom_vline(xintercept = as.Date("2021-06-07"), color = "#E41A1C", linetype = "dashed")+
  geom_vline(xintercept = as.Date("2021-07-13"), color = "#E41A1C", linetype = "dashed")+
  geom_vline(xintercept = as.Date("2021-08-09"), color = "#E41A1C", linetype = "dashed")+
  geom_line(data = GW_Levels, aes(date, level, color = id))+
  labs(x = "", y = "Groundwater level (m.a.s.l)")+
  scale_color_manual(values = c("#4DAF4A", "#F781BF", "#999999"), name = "")+
  annotate("text", x = as.Date("2020-10-20"), y = 371.2, label = "Upper Eiss", size = 10/.pt, color = "#4DAF4A")+
  annotate("text", x = as.Date("2020-10-20"), y = 369.4, label = "Lower Eiss", size = 10/.pt, color = "#F781BF")+
  annotate("text", x = as.Date("2020-10-13"), y = 364.2, label = "Morrill", size = 10/.pt, color = "#999999")+
  guides(color = "none")+
  scale_x_date(expand = c(0,0), limit=c(as.Date("2020-10-01"),as.Date("2021-09-30")),
               breaks=date_breaks("months"), labels=date_format("%b %y"))+
  scale_y_continuous(breaks = c(364, 366, 368, 370, 372))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Read: Precip at Konza HQ
Precip <- read_csv("data/LTER/Konza_Precip.csv") %>%
  filter(date >= date_start & date <= date_end) %>%
  mutate(date = ymd(date))

# Plot!
p3 <-
  ggplot()+
  geom_col(data = Precip, aes(date, precip_mm), color = "#984EA3")+
  geom_vline(xintercept = as.Date("2021-06-07"), color = "#E41A1C", linetype = 'dashed')+
  geom_vline(xintercept = as.Date("2021-07-13"), color = "#E41A1C", linetype = 'dashed')+
  geom_vline(xintercept = as.Date("2021-08-09"), color = "#E41A1C", linetype = 'dashed')+
  labs(x = "", y = "Precipitation (mm/d)")+
  labs(color = "")+
  scale_x_date(expand = c(0,0), limit=c(as.Date("2020-10-01"),as.Date("2021-09-30")),
               breaks=date_breaks("months"), labels=date_format("%b %y"))+
  scale_y_continuous(expand = expansion(mult = c(0,0.05)))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p3 + p1 + p2 + plot_layout(nrow = 3)
ggsave(path = "figures/", "Fig2.png", dpi=300, width = 190, height = 190, units = "mm")
ggsave(path = "figures/", "Fig2.svg", dpi=300, width = 190, height = 190, units = "mm")
ggsave(path = "figures/", "Figure2.pdf", dpi=300, width = 190, height = 190, units = "mm")