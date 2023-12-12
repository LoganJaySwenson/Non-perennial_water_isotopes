# Fig S7 & S8. Young water fractions with sinusoidal models
library(scales)
library(broom)
library(viridis)
library(lubridate)
library(patchwork)
library(tidyverse)

# Publication theme
source("code/Theme+Settings.R")

# Read: Streamflow at Kings Creek (USGS gage 06879650), Konza Prairie, KS
Flow <- dataRetrieval::readNWISdv(siteNumber = "06879650", parameterCd = "00060", startDate = "", endDate = "") %>%
  dataRetrieval::renameNWISColumns() %>%
  rename_with(.cols = everything(), tolower) %>%
  mutate(flow = (flow * 0.028316847)) %>% #cfs to m3/s
  select(site_no, date, flow) 

# Set minimum for plots
min_q <- 0.001
Flow <- Flow %>%
  mutate(flow_forlog = ifelse(flow < min_q, min_q, flow))

# Streamflow at Kings Creek (USGS 06879650) between Oct 2015 and Jan 2022. 
ggplot()+
  geom_line(data = Flow, aes(date, flow_forlog), color = "#377EB8")+
  labs(x = "", y = "Discharge (m\U00B3/s)")+
  labs(color = "")+
  scale_x_date(expand = c(0,0), limit=c(as.Date("2015-10-01"),as.Date("2022-01-01")),
               breaks=date_breaks("6 months"), labels=date_format("%b %y"))+
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    limits = c(0.001, 10),
    expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/KonzaPrairie_Streamflow_2015-2022.png", dpi=300, width = 190, height = 90, units = "mm")


# Stream δ18O at outlet
stream_isotopes <- as_tibble(read.csv("data/NEON/NEON_isotopes.csv")) %>%
  filter(type == "Stream") %>%
  mutate(date = mdy(date)) %>%
  group_by(date) %>%
  mutate(d18OWater = mean(d18OWater), d2HWater = mean(d2HWater)) %>%
  slice(1) %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  ungroup() %>%
  mutate(dexcess = d2HWater - (8 * d18OWater)) %>% # Dansgaard (1964): d-excess = δ2H - 8 * δ18O
  mutate(lcexcess = d2HWater - (7.9309 * d18OWater) - 10.2770) # Landwehr & Coplen (2006): line-conditioned excess = δ2H - m * δ18O - b

# Generate flow duration curve metrics
FDC <- Flow %>%
  arrange(desc(flow_forlog)) %>%
  mutate(flow_rank = seq(1:length(flow_forlog)),
         flow_prob = flow_rank/(length(flow_forlog)+1)*100,
         flow_percentile = 100 - flow_prob)

dates <- stream_isotopes$date
probs <- FDC %>%
  filter(date %in% dates) %>%
  select(date, flow, flow_forlog, flow_rank, flow_prob, flow_percentile)

stream_isotopes <- left_join(stream_isotopes, probs, by = "date")




# NEON precip δ18O 
precip_isotopes <- as_tibble(read.csv("data/NEON/NEON_isotopes.csv")) %>%
  filter(type == "Precipitation") %>%
  mutate(date = mdy(date)) %>%
  group_by(date) %>%
  mutate(d18OWater = mean(d18OWater), d2HWater = mean(d2HWater)) %>%
  slice(1) %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  ungroup() %>%
  mutate(t = as.numeric(date - ymd("2018-12-04"))) %>%
  select(date, t, d18OWater, d2HWater)

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

isotope_model <- function(t, A, phi, ks) {
  f <- (1/365.25)
  A * sin(2 * pi * f * t - phi) + ks
}

dates <- precip_isotopes$date
t <- precip_isotopes$t
y <- precip_isotopes$d18OWater
weights <- precip_isotopes$precip_sum_mm

date_range <- seq(min(dates), max(dates), by = "days")
time_range <- seq(min(t), max(t), by = 1)

p0 <- c(A = 2.5, phi = 0.0, ks = 0.0) # initial parameter guesses

fit <- nls(y ~ isotope_model(t, A, phi, ks), data = tibble(t, y), start = p0)
precip_coeffs <- tidy(fit)
precip_stats <- glance(fit)

AP_fit <- precip_coeffs %>% filter(term == "A") %>% pull(estimate)
AP_fit_se <- precip_coeffs %>% filter(term == "A") %>% pull(std.error)
phi_fit <- precip_coeffs %>% filter(term == "phi") %>% pull(estimate)
phi_fit_se <- precip_coeffs %>% filter(term == "phi") %>% pull(std.error)
ks_fit <- precip_coeffs %>% filter(term == "ks") %>% pull(estimate)
ks_fit_se <- precip_coeffs %>% filter(term == "ks") %>% pull(std.error)

y_hat <- predict(fit, tibble(t = time_range))

fit_weighted <- nls(y ~ isotope_model(t, A, phi, ks), data = tibble(t, y), start = p0, weights = weights)
precip_coeffs_weighted <- tidy(fit_weighted)
precip_stats_weighted <- glance(fit_weighted)

AP_fit_weighted <- precip_coeffs_weighted %>% filter(term == "A") %>% pull(estimate)
AP_fit_se_weighted <- precip_coeffs_weighted %>% filter(term == "A") %>% pull(std.error)
phi_fit_weighted <- precip_coeffs_weighted %>% filter(term == "phi") %>% pull(estimate)
phi_fit_se_weighted <- precip_coeffs_weighted %>% filter(term == "phi") %>% pull(std.error)
ks_fit_weighted <- precip_coeffs_weighted %>% filter(term == "ks") %>% pull(estimate)
ks_fit_se_weighted <- precip_coeffs_weighted %>% filter(term == "ks") %>% pull(std.error)

y_hat_weighted <- predict(fit_weighted, tibble(t = time_range))

parameters <- paste("AP: ", round(AP_fit, 2), "phi: ", round(phi_fit, 2), "ks: ", round(ks_fit, 2))
parameters_weighted <- paste("AP: ", round(AP_fit_weighted, 2), "phi: ", round(phi_fit_weighted, 2), "ks: ", round(ks_fit_weighted, 2))

ggplot()+
  geom_point(data = precip_isotopes, aes(date, d18OWater), color="#984EA3")+
  geom_line(data = tibble(date = date_range, y = y_hat), aes(x = date, y = y_hat))+
  geom_line(data = tibble(date = date_range, y = y_hat_weighted), aes(x = date, y = y_hat_weighted), color="#e6194b")+
  annotate("text", x = as.Date("2019-04-01"), y = 0, label = parameters, size = 8/.pt)+
  annotate("text", x = as.Date("2019-04-01"), y = -1, label = parameters_weighted, size = 8/.pt, color = "#e6194b")+
  labs(x = "", y = "Precipitation δ¹⁸O (‰)")+
  scale_x_date(breaks=date_breaks("years"), labels=date_format("%Y"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Perform Gaussian error propagation
set.seed(1)
n_iter <- 1000

precip_y_hats <- tibble()

for (i in 1:n_iter){
  
  AP_sampled <- rnorm(1, mean = AP_fit, sd = AP_fit_se)
  phi_sampled <- rnorm(1, mean = phi_fit, sd = phi_fit_se)
  ks_sampled <- rnorm(1, mean = ks_fit, sd = ks_fit_se)
  
  precip_y_hat <- isotope_model(t = time_range, A = AP_sampled, phi = phi_sampled, ks = ks_sampled)
  
  precip_y_hats <- rbind(precip_y_hats, tibble(date = date_range, y_hat = precip_y_hat, AP = AP_sampled, sim = i))
  
}

precip_y_hats_weighted <- tibble()

for (i in 1:n_iter){
  
  AP_sampled_weighted <- rnorm(1, mean = AP_fit_weighted, sd = AP_fit_se_weighted)
  phi_sampled_weighted <- rnorm(1, mean = phi_fit_weighted, sd = phi_fit_se_weighted)
  ks_sampled_weighted <- rnorm(1, mean = ks_fit_weighted, sd = ks_fit_se_weighted)
  
  precip_y_hat_weighted <- isotope_model(t = time_range, A = AP_sampled_weighted, phi = phi_sampled_weighted, ks = ks_sampled_weighted)
  
  precip_y_hats_weighted <- rbind(precip_y_hats_weighted, tibble(date = date_range, y_hat = precip_y_hat_weighted, AP = AP_sampled_weighted, sim = i))
  
}

ggplot()+
  geom_line(data = precip_y_hats, aes(date, y_hat, group = sim), color = "#999999", alpha = 0.01)+
  geom_line(data = precip_y_hats_weighted, aes(date, y_hat, group = sim), color = "#e6194b", alpha = 0.01)+
  geom_point(data = precip_isotopes, aes(date, d18OWater), color="#984EA3")+
  geom_line(data = tibble(date = date_range, y = y_hat), aes(x = date, y = y_hat))+
  geom_line(data = tibble(date = date_range, y = y_hat_weighted), aes(x = date, y = y_hat_weighted), color="#e6194b")+
  annotate("text", x = as.Date("2019-04-01"), y = 0, label = parameters, size = 8/.pt)+
  annotate("text", x = as.Date("2019-04-01"), y = -1, label = parameters_weighted, size = 8/.pt, color = "#e6194b")+
  labs(x = "", y = "Precipitation δ¹⁸O (‰)")+
  scale_x_date(breaks=date_breaks("years"), labels=date_format("%Y"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/FigS7.png", dpi=300, width = 190, height = 90, units = "mm")




# NEON stream δ18O at outlet
stream_isotopes <- stream_isotopes %>%
  mutate(t = as.numeric(date - ymd("2015-10-06")))

dates <- stream_isotopes$date
t <- stream_isotopes$t
y <- stream_isotopes$d18OWater
weights <- stream_isotopes$flow_forlog

date_range <- seq(min(dates), max(dates), by = "days")
time_range <- seq(min(t), max(t), by = 1)

p0 <- c(A = 1.0, phi = 0.0, ks = 0.0) # initial parameter guesses

fit <- nls(y ~ isotope_model(t, A, phi, ks), data = tibble(t, y), start = p0)
stream_coeffs <- tidy(fit)
stream_stats <- glance(fit)

AS_fit <- stream_coeffs %>% filter(term == "A") %>% pull(estimate)
AS_fit_se <- stream_coeffs %>% filter(term == "A") %>% pull(std.error)
phi_fit <- stream_coeffs %>% filter(term == "phi") %>% pull(estimate)
phi_fit_se <- stream_coeffs %>% filter(term == "phi") %>% pull(std.error)
ks_fit  <- stream_coeffs %>% filter(term == "ks") %>% pull(estimate)
ks_fit_se <- stream_coeffs %>% filter(term == "ks") %>% pull(std.error)

y_hat <- predict(fit, tibble(t = time_range))

fit_weighted <- nls(y ~ isotope_model(t, A, phi, ks), data = tibble(t, y), start = p0, weights = weights)
stream_coeffs_weighted <- tidy(fit_weighted)
stream_stats_weighted <- glance(fit_weighted)

AS_fit_weighted <- stream_coeffs_weighted %>% filter(term == "A") %>% pull(estimate)
AS_fit_se_weighted <- stream_coeffs_weighted %>% filter(term == "A") %>% pull(std.error)
phi_fit_weighted <- stream_coeffs_weighted %>% filter(term == "phi") %>% pull(estimate)
phi_fit_se_weighted <- stream_coeffs_weighted %>% filter(term == "phi") %>% pull(std.error)
ks_fit_weighted <- stream_coeffs_weighted %>% filter(term == "ks") %>% pull(estimate)
ks_fit_se_weighted <- stream_coeffs_weighted %>% filter(term == "ks") %>% pull(std.error)

y_hat_weighted <- predict(fit_weighted, tibble(t = time_range))

parameters <- paste("AS: ", round(AS_fit, 2), "phi: ", round(phi_fit, 2), "ks: ", round(ks_fit, 2))
parameters_weighted <- paste("AS: ", round(AS_fit_weighted, 2), "phi: ", round(phi_fit_weighted, 2), "ks: ", round(ks_fit_weighted, 2))

ggplot()+
  geom_point(data = stream_isotopes, aes(date, d18OWater), color="#984EA3")+
  geom_line(data = tibble(date = date_range, y = y_hat), aes(x = date, y = y_hat))+
  geom_line(data = tibble(date = date_range, y = y_hat_weighted), aes(x = date, y = y_hat_weighted), color="#e6194b")+
  annotate("text", x = as.Date("2016-06-01"), y = -3.4, label = parameters, size = 8/.pt)+
  annotate("text", x = as.Date("2016-06-01"), y = -3.58, label = parameters_weighted, size = 8/.pt, color = "#e6194b")+
  labs(x = "", y = "Stream δ¹⁸O (‰)")+
  scale_x_date(breaks=date_breaks("years"), labels=date_format("%Y"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Perform Gaussian error propagation
set.seed(1)
n_iter <- 1000

stream_y_hats <- tibble()

for (i in 1:n_iter){
  
  AS_sampled <- rnorm(1, mean = AS_fit, sd = AS_fit_se)
  phi_sampled <- rnorm(1, mean = phi_fit, sd = phi_fit_se)
  ks_sampled <- rnorm(1, mean = ks_fit, sd = ks_fit_se)
  
  stream_y_hat <- isotope_model(t = time_range, A = AS_sampled, phi = phi_sampled, ks = ks_sampled)
  
  stream_y_hats <- rbind(stream_y_hats, tibble(date = date_range, y_hat = stream_y_hat, AS = AS_sampled, sim = i))
  
}

stream_y_hats_weighted <- tibble()

for (i in 1:n_iter){
  
  AS_sampled_weighted <- rnorm(1, mean = AS_fit_weighted, sd = AS_fit_se_weighted)
  phi_sampled_weighted <- rnorm(1, mean = phi_fit_weighted, sd = phi_fit_se_weighted)
  ks_sampled_weighted <- rnorm(1, mean = ks_fit_weighted, sd = ks_fit_se_weighted)
  
  stream_y_hat_weighted <- isotope_model(t = time_range, A = AS_sampled_weighted, phi = phi_sampled_weighted, ks = ks_sampled_weighted)
  
  stream_y_hats_weighted <- rbind(stream_y_hats_weighted, tibble(date = date_range, y_hat = stream_y_hat_weighted, AS = AS_sampled_weighted, sim = i))
  
}

ggplot()+
  geom_line(data = stream_y_hats, aes(date, y_hat, group = sim), color = "#999999", alpha = 0.01)+
  geom_line(data = stream_y_hats_weighted, aes(date, y_hat, group = sim), color = "#e6194b", alpha = 0.01)+
  geom_point(data = stream_isotopes, aes(date, d18OWater), color="#984EA3")+
  geom_line(data = tibble(date = date_range, y = y_hat), aes(x = date, y = y_hat))+
  geom_line(data = tibble(date = date_range, y = y_hat_weighted), aes(x = date, y = y_hat_weighted), color="#e6194b")+
  annotate("text", x = as.Date("2016-06-01"), y = -3.4, label = parameters, size = 8/.pt)+
  annotate("text", x = as.Date("2016-06-01"), y = -3.58, label = parameters_weighted, size = 8/.pt, color = "#e6194b")+
  labs(x = "", y = "Stream δ¹⁸O (‰)")+
  scale_x_date(breaks=date_breaks("years"), labels=date_format("%Y"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/FigS8.png", dpi=300, width = 190, height = 90, units = "mm")


# Time-weighted young water fraction
AP <- precip_y_hats %>%
  group_by(sim) %>%
  distinct(sim, AP) %>%
  pull(AP)

AS <- stream_y_hats %>%
  group_by(sim) %>%
  distinct(sim, AS) %>%
  pull(AS)

FYW_time_weighted <- abs((AS/AP))


# Flow-weighted young water fraction
AP_weighted <- precip_y_hats_weighted %>%
  group_by(sim) %>%
  distinct(sim, AP) %>%
  pull(AP)

AS_weighted <- stream_y_hats_weighted %>%
  group_by(sim) %>%
  distinct(sim, AS) %>%
  pull(AS)

FYW_flow_weighted <- abs((AS_weighted/AP_weighted))

# Get coefficients & young water estimates
FYW <- tibble(id = c("time_weighted", "flow_weighted"),
              FYW = c(mean(FYW_time_weighted), mean(FYW_flow_weighted)),
              FYW_sd = c(sd(FYW_time_weighted), sd(FYW_flow_weighted)),
              AP = c(AP_fit, AP_fit_weighted),
              AP_se = c(AP_fit_se, AP_fit_se_weighted),
              phi = c(precip_coeffs %>% filter(term == "phi") %>% pull(estimate), precip_coeffs_weighted %>% filter(term == "phi") %>% pull(estimate)),
              phi_se = c(precip_coeffs %>% filter(term == "phi") %>% pull(std.error), precip_coeffs_weighted %>% filter(term == "phi") %>% pull(std.error)),
              ks = c(precip_coeffs %>% filter(term == "ks") %>% pull(estimate), precip_coeffs_weighted %>% filter(term == "ks") %>% pull(estimate)),
              ks_se = c(precip_coeffs %>% filter(term == "ks") %>% pull(std.error), precip_coeffs_weighted %>% filter(term == "ks") %>% pull(std.error)),
              n = c(length(stream_isotopes$d18OWater), length(stream_isotopes$d18OWater)))