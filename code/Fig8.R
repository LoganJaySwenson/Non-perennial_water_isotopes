# Fig. 8. Sensitivity of young water fractions to flow regime
# comparison of Bayesian and sinusoidal models
library(scales)
library(broom)
library(viridis)
library(lubridate)
library(tidyverse)

# Publication theme
source("code/Theme+Settings.R")

area <- 1306 # hectares

# Read: Kings Creek (USGS 06879650), Konza Prairie, KS
Flow <- dataRetrieval::readNWISdv(siteNumber = "06879650", parameterCd = "00060", startDate = "", endDate = "") %>%
  dataRetrieval::renameNWISColumns() %>%
  rename_with(.cols = everything(), tolower) %>%
  mutate(flow = flow * 0.028316847, # cfs to m3/s
         norm_flow = ((flow*1e+9)/(area*1e+10)) * 86400) %>% # m3/s to mm/d
  select(site_no, date, flow, norm_flow) # flow is m3/s norm_flow is mm/d

# Set minimum for plots
Flow <- mutate(Flow, flow_forlog = flow + 0.001)

# NEON δ18O at outlet
stream_isotopes <- read_csv("data/NEON/NEON_isotopes.csv") %>%
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
  select(date, flow, flow_forlog, flow_rank, flow_prob, flow_percentile, norm_flow)

stream_isotopes <- left_join(stream_isotopes, probs, by = "date")

date_start <- min(stream_isotopes$date)
date_end <- max(stream_isotopes$date)
stream_isotopes <- mutate(stream_isotopes, t = as.numeric(date - date_start)) # set a numeric time column for annual freq.

date_start <- as.Date("2020-10-01")
date_end <- as.Date("2021-09-30")
stream_isotopes_2021 <- filter(stream_isotopes, date >= date_start & date <= date_end) 

date_start <- min(stream_isotopes$date)
date_end <- max(stream_isotopes$date)
Flow_sub <- filter(Flow, date >= date_start & date <= date_end)




# 2. Normalized flow (m/d)

# Define variables for Eqs. 
t <- stream_isotopes$t
y <- stream_isotopes$d18OWater
Q <- stream_isotopes$norm_flow

# Eq. 8 in Gallart et al (2021)
# https://doi.org/10.5194/hess-24-1101-2020
isotope_model <- function(Q, t, F0, Sd, phi, ks) {
  AP <- -2.84
  f <- (1/365)
  AP * (1 - (1 - F0) * exp(-Q * Sd)) * sin(2 * pi * f * t - phi) + ks
}

p0 <- c(F0 = 0.1, Sd = 0.1, phi = 0.0, ks = 0.0) # initial parameter guesses

fit <- nls(y ~ isotope_model(Q, t, F0, Sd, phi, ks), data = tibble(Q, t, y), start = p0)
stream_coeffs <- tidy(fit)
stream_stats <- glance(fit)

head(stream_coeffs)

F0_fit <- stream_coeffs %>% filter(term == "F0") %>% pull(estimate)
Sd_fit <- stream_coeffs %>% filter(term == "Sd") %>% pull(estimate)
phi_fit <- stream_coeffs %>% filter(term == "phi") %>% pull(estimate)
ks_fit  <- stream_coeffs %>% filter(term == "ks") %>% pull(estimate)

# Eq. 6 in Gallart et al (2020) (F0 and Sd are both estimated from Eq. 8 above)
FYW_model <- function(Q, F0, Sd){
  1 - ((1 - F0) * exp(-Q * Sd))
}

F0 <- abs(F0_fit)
F0_se <- stream_coeffs %>% filter(term == "F0") %>% pull(std.error)
Sd <- abs(Sd_fit)
Sd_se <- stream_coeffs %>% filter(term == "F0") %>% pull(std.error)

Q_min <- min(Flow_sub$norm_flow)
Q_max <- max(Flow_sub$norm_flow)
Q_range <- seq(Q_min, Q_max, length.out = 1e5)

FYW_curve_norm <- tibble(flow = Q_range, 
                         FYW = FYW_model(Q = Q_range, F0 = F0, Sd = Sd) * 100)

# 2021 Water Year

# Define variables for Eqs. 
t <- stream_isotopes_2021$t
y <- stream_isotopes_2021$d18OWater
Q <- stream_isotopes_2021$norm_flow

p0 <- c(F0 = 0.1, Sd = 0.1, phi = 0.0, ks = 0.0) # initial parameter guesses

fit <- nls(y ~ isotope_model(Q, t, F0, Sd, phi, ks), data = tibble(Q, t, y), start = p0)
stream_coeffs_2021 <- tidy(fit)
stream_stats_2021 <- glance(fit)

head(stream_coeffs_2021)

F0_fit <- stream_coeffs_2021 %>% filter(term == "F0") %>% pull(estimate)
Sd_fit <- stream_coeffs_2021 %>% filter(term == "Sd") %>% pull(estimate)
phi_fit <- stream_coeffs_2021 %>% filter(term == "phi") %>% pull(estimate)
ks_fit  <- stream_coeffs_2021 %>% filter(term == "ks") %>% pull(estimate)

F0_2021 <- abs(F0_fit)
F0_2021_se <- stream_coeffs_2021 %>% filter(term == "F0") %>% pull(std.error)
Sd_2021 <- abs(Sd_fit)
Sd_2021_se <- stream_coeffs_2021 %>% filter(term == "Sd") %>% pull(std.error)

Q_min <- min(Flow_sub$norm_flow)
Q_max <- max(Flow_sub$norm_flow)
Q_range <- seq(Q_min, Q_max, length.out = 1e5)

FYW_curve_norm_2021 <- tibble(flow = Q_range, 
                              FYW = FYW_model(Q = Q_range, F0 = F0_2021, Sd = Sd_2021) * 100)

# Read: NEON δ18O timeseries (Bayesian estimates)
date_start <- as.Date("2020-10-01")
date_end <- as.Date("2021-09-30")

Bayesian_estimates <- read_csv("data/NEON_isotopes_results.csv") %>%
  mutate(date = ymd(date)) %>%
  filter(date >= date_start & date <= date_end) %>%
  na.omit()

dates <- Bayesian_estimates$date
probs <- FDC %>%
  filter(date %in% dates) %>%
  select(date, flow, flow_rank, flow_prob, flow_percentile, norm_flow) %>%
  mutate(flow_forlog = flow + 0.001)
Bayesian_estimates <- left_join(Bayesian_estimates, probs, by = "date")

# Eq. 6 in Gallart et al (2020) fit to Bayesian estimates (F0 and Sd are both estimated)

p0 <- c(F0 = F0, Sd = Sd) # initial parameter guesses

y <- Bayesian_estimates$s1 * 0.01
Q <- Bayesian_estimates$norm_flow

fit <- nls(y ~ FYW_model(Q, F0, Sd), data = tibble(Q), start = p0)
coeffs <- tidy(fit)
stats <- glance(fit)

head(coeffs)

F0_fit <- coeffs %>% filter(term == "F0") %>% pull(estimate)
F0_fit_se <- coeffs %>% filter(term == "F0") %>% pull(std.error)
Sd_fit <- coeffs %>% filter(term == "Sd") %>% pull(estimate)
Sd_fit_se <- coeffs %>% filter(term == "Sd") %>% pull(std.error)

Q_min <- min(Flow_sub$norm_flow)
Q_max <- max(Flow_sub$norm_flow)
Q_range <- seq(Q_min, Q_max, length.out = 1e5)

Bayesian_FYW_curve_norm <- tibble(flow = Q_range, 
                                  FYW = FYW_model(Q = Q_range, F0 = F0_fit, Sd = Sd_fit) * 100)


F0_points <- tibble(F0 = c(F0*100, F0_2021*100, F0_fit*100),
                    flow = c(Q_min,Q_min,Q_min))

parameters <- paste("F0: ", round(F0 * 100, digits=1),"%", "    Sd: ", round(Sd, 2))
parameters_2021 <- paste("F0: ", round(F0_2021 * 100, digits=1),"%", "    Sd: ", round(Sd_2021, 2))
Bayesian_parameters <- paste("F0: ", round(F0_fit * 100, digits=1),"%", "    Sd: ", round(Sd_fit, 2))


Q_min <- min(Flow_sub$norm_flow)
Q_max <- max(Flow_sub$norm_flow)
Q_range <- seq(Q_min, Q_max, length.out = 1e3)


# Perform Gaussian error propagation
set.seed(1)
n_iter <- 100

FYW_curve_norm_uncertainty <- tibble()

for (i in 1:n_iter){
  
  F0_sampled <- rnorm(1, mean = F0, sd = F0_se)
  Sd_sampled <- rnorm(1, mean = Sd, sd = Sd_se)
  
  temp <- tibble(flow = Q_range, 
                 FYW = FYW_model(Q = Q_range, F0 = F0_sampled, Sd = Sd_sampled) * 100)
  
  FYW_curve_norm_uncertainty <- rbind(FYW_curve_norm_uncertainty, tibble(temp, sim = i))
  
}

# Perform Gaussian error propagation
set.seed(1)
n_iter <- 100

FYW_curve_norm_2021_uncertainty <- tibble()

for (i in 1:n_iter){
  
  F0_sampled <- rnorm(1, mean = F0_2021, sd = F0_2021_se)
  Sd_sampled <- rnorm(1, mean = Sd_2021, sd = Sd_2021_se)
  
  temp <- tibble(flow = Q_range, 
                 FYW = FYW_model(Q = Q_range, F0 = F0_sampled, Sd = Sd_sampled) * 100)
  
  FYW_curve_norm_2021_uncertainty <- rbind(FYW_curve_norm_2021_uncertainty, tibble(temp, sim = i))
  
}

# Perform Gaussian error propagation
set.seed(1)
n_iter <- 100

Bayesian_FYW_curve_norm_uncertainty <- tibble()

for (i in 1:n_iter){
  
  F0_sampled <- rnorm(1, mean = F0_fit, sd = F0_fit_se)
  Sd_sampled <- rnorm(1, mean = Sd_fit, sd = Sd_fit_se)
  
  temp <- tibble(flow = Q_range, 
                 FYW = FYW_model(Q = Q_range, F0 = F0_sampled, Sd = Sd_sampled) * 100)
  
  Bayesian_FYW_curve_norm_uncertainty <- rbind(Bayesian_FYW_curve_norm_uncertainty, tibble(temp, sim = i))
  
}

FYW_curve_norm_uncertainty <- FYW_curve_norm_uncertainty %>%
  group_by(sim) %>%
  filter(all(FYW >= 0))

FYW_curve_norm_2021_uncertainty <- FYW_curve_norm_2021_uncertainty %>%
  group_by(sim) %>%
  filter(all(FYW >= 0))

ggplot()+
  geom_line(data = FYW_curve_norm_uncertainty, aes(flow, FYW, group = sim), color = "#999999", alpha = 0.01)+
  geom_line(data = FYW_curve_norm, aes(flow, FYW), linetype = "dashed", color = "#000000")+
  geom_line(data = FYW_curve_norm_2021_uncertainty, aes(flow, FYW, group = sim), color = "#999999", alpha = 0.01)+
  geom_line(data = FYW_curve_norm_2021, aes(flow, FYW), linetype = "dashed", color = "#000000")+
  geom_line(data = Bayesian_FYW_curve_norm, aes(flow, FYW), linetype = "dashed", color = "#e6194b")+
  geom_line(data = Bayesian_FYW_curve_norm_uncertainty, aes(flow, FYW, group = sim), color = "#e6194b", alpha = 0.01)+
  geom_point(data = Bayesian_estimates, aes(norm_flow, s1), fill = "#999999", color = "#000000", pch = 21, size = 2.5)+
  geom_point(data = F0_points, aes(flow, F0), fill = c("#000000", "#000000", "#e6194b"), color = "#000000", pch = 21, size = 2.5)+
  scale_x_continuous(limits = c(0, 80), breaks = c(seq(0,80,by=20)))+
  #annotate("text", x = 60, y = 90, label = Bayesian_parameters, size = 8/.pt, color = "#e6194b")+
  #annotate("text", x = 60, y = 85, label = parameters_2021, size = 8/.pt)+
  #annotate("text", x = 60, y = 80, label = parameters, size = 8/.pt)+
  labs(x = "Discharge (mm/d)", y = expression(F["yw"]* " (%)"))
ggsave("figures/Fig8.png", dpi=300, width = 110, height = 90, units = "mm")
ggsave("figures/Figure8.pdf", dpi=300, width = 110, height = 90, units = "mm")