# Fig S7 & S8. Young water fractions with sinusoidal models 
#              Gaussian error propagation accounting for parameter correlations between A, φ, and k.
library(tidyverse)

source(file.path("code", "theme.R"))

# Kirchner, J. W. (2016a). Aggregation in environmental systems Part 1: 
# Seasonal tracer cycles quantify young water fractions, but not mean transit times, in spatially heterogeneous catchments. 
# Hydrology and Earth System Sciences, 20(1), 279–297. https://doi.org/10.5194/hess-20-279-2016

# Kirchner, J. W. (2016b). Aggregation in environmental systems Part 2: 
# Catchment mean transit times and young water fractions under hydrologic nonstationarity. 
# Hydrology and Earth System Sciences, 20(1), 299–328. https://doi.org/10.5194/hess-20-299-2016

set.seed(1)

# eq. 4 in Kirchner (2016)
isotope_model <- function(t, A, phi, ks) {
  f <- (1/365.25)
  A * sin(2 * pi * f * t - phi) + ks
}

# eq. 8 in Kirchner (2016)
Fyw_model <- function(AS, AP){
  abs(AS/AP)
}

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

# NEON stream δ18O at outlet
NEON_stream_isotopes <- NEON_isotopes %>%
  filter(type == "SW", date <= end) %>% 
  mutate(
    t = as.numeric(date - date[1]),
  )

# Precip at HQ weather station
precip <- file.path("data", "LTER", "LTER_HQ_weather_station.csv") %>%
  read_csv(show_col_types = F, col_select = c(siteID, date, precip_mm)) %>%
  filter(date >= min(NEON_precip_isotopes$date) & date <= end) 

# Streamflow at USGS 06879650
flow <- file.path("data", "USGS", "USGS_06879650.csv") %>% 
  read_csv(show_col_types = F) %>%
  filter(date >= min(NEON_stream_isotopes$date) & date <= end) 

# precip weights
rows_match <- match(NEON_precip_isotopes$date, precip$date)
NEON_precip_isotopes$precip_sum_mm <- NA
for (i in 2:length(rows_match)){
  row_start_sum <- rows_match[i-1] + 1
  row_end_sum <- rows_match[i]
  NEON_precip_isotopes$precip_sum_mm[i] <- sum(precip$precip_mm[row_start_sum:row_end_sum],na.rm=T)
}
NEON_precip_isotopes <- na.omit(NEON_precip_isotopes)

# streamflow weights
NEON_stream_isotopes <- NEON_stream_isotopes %>%
  left_join(select(flow, date, flow_m3s), by="date") %>%
  mutate(
    flow_m3s = if_else(flow_m3s >= min_q, flow_m3s, min_q),
  )
NEON_stream_isotopes <- na.omit(NEON_stream_isotopes)



# Precipitation δ18O parameters
dates <- NEON_precip_isotopes$date
t <- NEON_precip_isotopes$t
y <- NEON_precip_isotopes$d18O
weights <- NEON_precip_isotopes$precip_sum_mm

date_range <- seq(min(dates), max(dates), by = "days")
time_range <- seq(min(t), max(t), by = 1)

p0 <- c(A = 2.5, phi = 0.0, ks = 0.0) # initial parameter guesses

# fit isotope model
fit <- nls(y ~ isotope_model(t, A, phi, ks), data = list(t=t, y=y), start = p0)
fit_weighted <- nls(y ~ isotope_model(t, A, phi, ks), data = list(t=t, y=y), start = p0, weights = weights)

# fitted parameters
coefs <- summary(fit)$coefficients
AP     <- coefs["A", "Estimate"]
AP_se  <- coefs["A", "Std. Error"]
phi    <- coefs["phi", "Estimate"]
phi_se <- coefs["phi", "Std. Error"]
ks     <- coefs["ks", "Estimate"]
ks_se  <- coefs["ks", "Std. Error"]

coefs_weighted  <- summary(fit_weighted)$coefficients
AP_weighted     <- coefs_weighted["A", "Estimate"]
AP_se_weighted  <- coefs_weighted["A", "Std. Error"]
phi_weighted    <- coefs_weighted["phi", "Estimate"]
phi_se_weighted <- coefs_weighted["phi", "Std. Error"]
ks_weighted     <- coefs_weighted["ks", "Estimate"]
ks_se_weighted  <- coefs_weighted["ks", "Std. Error"]

# parameter covariance matrix
vc <- vcov(fit)
vc_weighted <- vcov(fit_weighted)

# Perform Gaussian error propagation
n = 1e3

out <- vector("list", n)
for (i in 1:n){
  params <- MASS::mvrnorm(1, coef(fit), vc)
  
  out[[i]] <- tibble(
    n = i,
    date = date_range,
    t = time_range,
    y = isotope_model(t = time_range, A = params["A"], phi = params["phi"], ks = params["ks"]),
    AP = params["A"],
  )
  
}
precip_fit <- tibble(date = date_range, y = predict(fit, tibble(t = time_range)))
precip_fits <- bind_rows(out)

out <- vector("list", n)
for (i in 1:n){
  params <- MASS::mvrnorm(1, coef(fit_weighted), vc_weighted)
  
  out[[i]] <- tibble(
    n = i,
    date = date_range,
    t = time_range,
    y = isotope_model(t = time_range, A = params["A"], phi = params["phi"], ks = params["ks"]),
    AP = params["A"],
  )
  
}
precip_fit_weighted <- tibble(date = date_range, y = predict(fit_weighted, tibble(t = time_range)))
precip_fits_weighted <- bind_rows(out)

parameters <- paste("AP: ", round(AP, 2), "phi: ", round(phi, 2), "ks: ", round(ks, 2))
parameters_weighted <- paste("AP: ", round(AP_weighted, 2), "phi: ", round(phi_weighted, 2), "ks: ", round(ks_weighted, 2))

p1 = ggplot()+
  geom_line(data = precip_fits, aes(date, y, group = n), color = "#999999", alpha = 0.01)+
  geom_line(data = precip_fits_weighted, aes(date, y, group = n), color = "#E6194B", alpha = 0.01)+
  geom_line(data = precip_fit, aes(date, y), color = "#373737", lwd=0.6)+
  geom_line(data = precip_fit_weighted, aes(date, y), color = "#E6194B", lwd=0.6)+
  geom_point(data = NEON_precip_isotopes, aes(date, d18O), fill = "#984EA3", color = "#373737", pch=21, size=2.5)+
  annotate("text", x=as.Date("2018-12-04"), y=0.5, label=parameters, color="#373737", size=10/.pt, hjust=0)+
  annotate("text", x=as.Date("2018-12-04"), y=-0.5, label=parameters_weighted, color="#E6194B", size=10/.pt, hjust=0)+
  labs(x = "", y = expression(delta^{18}*"O (‰)"))+
  scale_x_date(
    breaks = scales::date_breaks("years"),
    labels = scales::date_format("%Y")
  )+
  theme(
    axis.title = element_text(size=12)
  )
ggsave(file.path("figures", "FigS7.png"), plot=p1, dpi=300, width=190, height=90, units="mm")



# Stream δ18O parameters
dates <- NEON_stream_isotopes$date
t <- NEON_stream_isotopes$t
y <- NEON_stream_isotopes$d18O
weights <- NEON_stream_isotopes$flow_m3s

date_range <- seq(min(dates), max(dates), by = "days")
time_range <- seq(min(t), max(t), by = 1)

p0 <- c(A = 1.0, phi = 0.0, ks = 0.0) # initial parameter guesses

# fit isotope model
fit <- nls(y ~ isotope_model(t, A, phi, ks), data = list(t=t, y=y), start = p0)
fit_weighted <- nls(y ~ isotope_model(t, A, phi, ks), data = list(t=t, y=y), start = p0, weights = weights)

# fitted parameters
coefs <- summary(fit)$coefficients
AS     <- coefs["A", "Estimate"]
AS_se  <- coefs["A", "Std. Error"]
phi    <- coefs["phi", "Estimate"]
phi_se <- coefs["phi", "Std. Error"]
ks     <- coefs["ks", "Estimate"]
ks_se  <- coefs["ks", "Std. Error"]

coefs_weighted  <- summary(fit_weighted)$coefficients
AS_weighted     <- coefs_weighted["A", "Estimate"]
AS_se_weighted  <- coefs_weighted["A", "Std. Error"]
phi_weighted    <- coefs_weighted["phi", "Estimate"]
phi_se_weighted <- coefs_weighted["phi", "Std. Error"]
ks_weighted     <- coefs_weighted["ks", "Estimate"]
ks_se_weighted  <- coefs_weighted["ks", "Std. Error"]

# parameter covariance matrix
vc <- vcov(fit)
vc_weighted <- vcov(fit_weighted)

out <- vector("list", n)
for (i in 1:n){
  
  params <- MASS::mvrnorm(1, coef(fit), vc)
  
  out[[i]] <- tibble(
    n = i,
    date = date_range,
    t = time_range,
    y = isotope_model(t = time_range, A = params["A"], phi = params["phi"], ks = params["ks"]),
    AS = params["A"],
  )
  
}
stream_fit <- tibble(date = date_range, y = predict(fit, tibble(t = time_range)))
stream_fits <- bind_rows(out)

out <- vector("list", n)
for (i in 1:n){
  
  params <- MASS::mvrnorm(1, coef(fit_weighted), vc_weighted)
  
  out[[i]] <- tibble(
    n = i,
    date = date_range,
    t = time_range,
    y = isotope_model(t = time_range, A = params["A"], phi = params["phi"], ks = params["ks"]),
    AS = params["A"],
  )
  
}
stream_fit_weighted <- tibble(date = date_range, y = predict(fit_weighted, tibble(t = time_range)))
stream_fits_weighted <- bind_rows(out)

parameters <- paste("AS: ", round(AS, 2), "phi: ", round(phi, 2), "ks: ", round(ks, 2))
parameters_weighted <- paste("AS: ", round(AS_weighted, 2), "phi: ", round(phi_weighted, 2), "ks: ", round(ks_weighted, 2))

p2 = ggplot()+
  geom_line(data = stream_fits, aes(date, y, group = n), color = "#999999", alpha = 0.01)+
  geom_line(data = stream_fits_weighted, aes(date, y, group = n), color = "#E6194B", alpha = 0.01)+
  geom_line(data = stream_fit, aes(date, y), color = "#373737", lwd=0.6)+
  geom_line(data = stream_fit_weighted, aes(date, y), color = "#E6194B", lwd=0.6)+
  geom_point(data = NEON_stream_isotopes, aes(date, d18O), fill = "#984EA3", color = "#373737", pch=21, size=2.5)+
  annotate("text", x=as.Date("2015-10-06"), y=-3.4, label=parameters, color="#373737", size=10/.pt, hjust=0)+
  annotate("text", x=as.Date("2015-10-06"), y=-3.58, label=parameters_weighted, color="#E6194B", size=10/.pt, hjust=0)+
  labs(x = "", y = expression(delta^{18}*"O (‰)"))+
  scale_x_date(
    breaks = scales::date_breaks("years"),
    labels = scales::date_format("%Y")
  )+
  theme(
    axis.title = element_text(size=12)
  )
ggsave(file.path("figures", "FigS8.png"), plot=p2, dpi=300, width=190, height=90, units="mm")

print(paste0("Time-weighted Fyw: ", round(Fyw_model(AS, AP)*100, 2), "%"))
print(paste0("Flow-weighted Fyw: ", round(Fyw_model(AS_weighted, AP_weighted)*100, 2), "%"))

precip_fits_parameters <- precip_fits %>%
  group_by(n) %>% 
  summarise(
    n = n[1],
    AP = AP[1]
  )

precip_fits_weighted_parameters <- precip_fits_weighted %>%
  group_by(n) %>% 
  summarise(
    n = n[1],
    AP_weighted = AP[1]
  )

stream_fits_parameters <- stream_fits %>%
  group_by(n) %>% 
  summarise(
    n = n[1],
    AS = AS[1]
  )

stream_fits_weighted_parameters <- stream_fits_weighted %>%
  group_by(n) %>% 
  summarise(
    n = n[1],
    AS_weighted = AS[1]
  )

fits <- list(precip_fits_parameters, precip_fits_weighted_parameters,
             stream_fits_parameters, stream_fits_weighted_parameters) %>%
  reduce(left_join, by="n") %>%
  mutate(
    Fyw = Fyw_model(AS, AP),
    Fyw_weighted = Fyw_model(AS_weighted, AP_weighted)
  )

print(paste0("Time-weighted Fyw: ", round(mean(fits$Fyw)*100, 1), "% ± ", round(sd(fits$Fyw)*100, 1)))
print(paste0("Flow-weighted Fyw: ", round(mean(fits$Fyw_weighted)*100, 1), "% ± ", round(sd(fits$Fyw_weighted)*100, 1)))
