#Fig 2. Site info
library(lubridate)
library(scales)
library(dataRetrieval)
library(patchwork)
library(tidyverse)

#Publication theme
source("code/Theme+Settings.R")

#Read: Discharge at Kings Creek (USGS gage 06879650), Manhattan, KS
Flow <- readNWISdv(siteNumber = "06879650", parameterCd = "00060", startDate = "", endDate = "")
Flow <- renameNWISColumns(Flow)
Flow <- Flow %>%
  filter(Date >= "2020-10-01" & Date <="2021-09-30") %>%
  mutate(Flow = (Flow * 0.028316847)) #cfs to metric

#Plot!
p1 <-
  ggplot()+
  geom_vline(xintercept = as.Date("2021-06-07"), color = "#E41A1C", linetype = "dashed")+
  geom_vline(xintercept = as.Date("2021-07-13"), color = "#E41A1C", linetype = "dashed")+
  geom_vline(xintercept = as.Date("2021-08-09"), color = "#E41A1C", linetype = "dashed")+
  geom_line(data = Flow, aes(Date, Flow), color = "#377EB8")+
  labs(x = "", y = "Discharge (m\U00B3/d)")+
  labs(color = "")+
  scale_x_date(expand = c(0,0), limit=c(as.Date("2020-10-01"),as.Date("2021-09-30")),
               breaks=date_breaks("months"), labels=date_format("%b %y"))+
  scale_y_continuous(expand = expansion(mult = c(0,0.05)))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Read: GW levels
Mor3_5 <- read.csv("data/LTER/GW_Levels_3-5Mor.csv")
Mor3_5$Date <- mdy(Mor3_5$Date)
Mor3_5 <- Mor3_5 %>%
  filter(Date >= "2020-10-01" & Date <="2021-09-30") %>%
  mutate(id = "Mor 3-5") %>%
  select(-c(Time, Date.Time))

Mor3_5_1 <- read.csv("data/LTER/GW_Levels_3-5-1Mor.csv")
Mor3_5_1$Date <- mdy(Mor3_5_1$Date)
Mor3_5_1 <- Mor3_5_1 %>%
  filter(Date >= "2020-10-01" & Date <="2021-09-30") %>%
  mutate(id = "Mor 3-5-1") %>%
  select(-c(Time, Date.Time))

Mor4_6 <- read.csv("data/LTER/GW_Levels_4-6Mor.csv")
Mor4_6$Date <- mdy(Mor4_6$Date)
Mor4_6 <- Mor4_6 %>%
  filter(Date >= "2020-10-01" & Date <="2021-09-30") %>%
  mutate(id = "Mor 4-6") %>%
  select(-c(Time, Date.Time))

LowerEis4_6 <- read.csv("data/LTER/GW_Levels_4-6Eis1.csv")
LowerEis4_6$Date <- mdy(LowerEis4_6$Date)
LowerEis4_6 <- LowerEis4_6 %>%
  filter(Date >= "2020-10-01" & Date <="2021-09-30") %>%
  mutate(id = "Lower Eis 4-6") %>%
  select(-c(Time, Date.Time))

UpperEis4_6 <- read.csv("data/LTER/GW_Levels_4-6Eis2.csv")
UpperEis4_6$Date <- mdy(UpperEis4_6$Date)
UpperEis4_6 <- UpperEis4_6 %>%
  filter(Date >= "2020-10-01" & Date <="2021-09-30") %>%
  mutate(id = "Upper Eis 4-6") %>%
  select(-c(Time, Date...Time))

#Bind rows: GW levels
GW_Levels <- bind_rows(Mor4_6, LowerEis4_6, UpperEis4_6)

#Plot!
GW_Levels$id <- factor(GW_Levels$id, levels = c("Upper Eis 4-6", "Lower Eis 4-6", "Mor 4-6"))
p2 <-
  ggplot()+
  geom_vline(xintercept = as.Date("2021-06-07"), color = "#E41A1C", linetype = "dashed")+
  geom_vline(xintercept = as.Date("2021-07-13"), color = "#E41A1C", linetype = "dashed")+
  geom_vline(xintercept = as.Date("2021-08-09"), color = "#E41A1C", linetype = "dashed")+
  geom_line(data = GW_Levels, aes(Date, LEVEL, color = id))+
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


#Read: Precipitation at Konza HQ
Precip <- read.csv("data/LTER/Konza_Precip.csv")
Precip <- Precip %>%
  filter(date >= "2020-10-01" & date <="2021-09-30") %>%
  mutate(date = ymd(date))

#Plot!
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