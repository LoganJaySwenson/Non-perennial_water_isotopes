#Fig S8. Uncertainty in FYW estimates!
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
base.model <- IRLS(P$d18OWater, as.matrix(P[,c("cosct","sinct")]), ww=P$precip_sum_mm, type="Cauchy")
resids <- base.model$fit$residuals
resids.sd <- sd(resids)

#Fit models
set.seed(1)
n_iter <- 10000
for (i in 1:n_iter){
  model <- IRLS(P$d18OWater + rnorm(n = length(P$d18OWater), mean = 0, sd = (max(P$d18OWater) - min(P$d18OWater)) * 0.05), as.matrix(P[,c("cosct","sinct")]), ww=P$precip_sum_mm, type="Cauchy")
  
  df <- data.frame(aP = model$fit$coefficients[2], bP = model$fit$coefficients[3])
  row.names(df) <- 1
  
  if (i == 1){
    df_fits_all <- df
  } else {
    df_fits_all <- bind_rows(df_fits_all, df)
  }
}
row.names(df_fits_all) <- 1:n_iter
P_model_fit <- df_fits_all

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
#Qw <- filter(Qw, date >= "2021-04-01" & date <= "2021-09-01") #remove all values with zero flow!
base.model <- IRLS(Qw$d18OWater, as.matrix(Qw[,c("cosct","sinct")]), ww=Qw$flow, type="Cauchy") 
resids <- base.model$fit$residuals
resids.wt <- base.model$wt
resids.sd <- sd(resids * resids.wt)

#Fit models
set.seed(1)
n_iter <- 10000

df_fits_all <- data.frame(aQ = rep(NA, n_iter),
                          bQ = rep(NA, n_iter))
for (i in 1:n_iter){
model_test =  try(model <- IRLS(Qw$d18OWater + rnorm(n = length(Qw$d18OWater), mean = 0, sd = (max(Q$d18OWater) - min(Q$d18OWater)) * 0.05), 
                                 as.matrix(Qw[,c("cosct","sinct")]), ww=Qw$flow, type="Cauchy"), silent = T) #wt = resids.sd
            if (class(model_test) == "try-error"){
              df_fits_all$aQ[i] = NA
              df_fits_all$bQ[i] = NA
            } else{
              df_fits_all$aQ[i] = model$fit$coefficients[2]
              df_fits_all$bQ[i] = model$fit$coefficients[3]
            }

             #if (i == 1){
               #df_fits_all <- df
             #} else {
               #df_fits_all <- bind_rows(df_fits_all, df)}
}
row.names(df_fits_all) <- 1:n_iter
Qw_model_fit <- df_fits_all

FYW_model <- bind_cols(P_model_fit, Qw_model_fit) %>%
  mutate(AP = sqrt(aP^2 + bP^2)) %>%
  mutate(AQ = sqrt(aQ^2 + bQ^2)) %>%
  mutate(FYW = (AQ/AP) *100) %>%
  filter(FYW >= 0 & FYW <= 100) #remove all impossible values!

#95% confidence intervals & median 
n <- length(FYW_model)
x <- mean(FYW_model$FYW)
s <- sd(FYW_model$FYW)

lowerinterval <- print(x - qt(0.975,df=n-1)*s/sqrt(n))
median <- print(median(FYW_model$FYW))
upperinterval <- print(x + qt(0.975,df=n-1)*s/sqrt(n))

#Plot distribution!
ggplot()+
  geom_histogram(data = FYW_model, aes(x = FYW, y = ..density..), binwidth = 1, colour = "black", fill = "white")+
  geom_density(data = FYW_model, aes(FYW), alpha = .5)+ 
  annotate("text", x = 25, y = 0.028, label = paste("Median:", round(median, digits = 2),"%"), size = 10/.pt)+
  annotate("text", x = 25, y = 0.030, label = paste("95% confidence interval:", round(lowerinterval, digits = 2),"% to", round(upperinterval, digits = 2), "%"), size = 10/.pt)+
  labs(x = "Young water fraction (%)", y = "Count")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = expansion(mult = c(0,0.05)))+
  theme(plot.margin=unit(c(5,5,5,5), "mm"))
ggsave(path = "figures/", "FigS8.png", dpi=300, width = 190, height = 90, units = "mm")