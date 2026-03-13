# Fig 6. Sensitivity of young water fractions to flow conditions
#        Gaussian error propagation accounting for parameter correlations between F0, Sd, φ, and k.
library(tidyverse)

source(file.path("code", "theme.R"))

# Gallart, F., von Freyberg, J., Valiente, M., Kirchner, J. W., Llorens, P., & Latron, J. (2020). 
# Technical note: An improved discharge sensitivity metric for young water fractions. 
# Hydrology and Earth System Sciences, 24(3), 1101–1107. https://doi.org/10.5194/hess-24-1101-2020

set.seed(1)

# eq. 8 in Gallart et al. (2020)
isotope_model <- function(Q, t, F0, Sd, phi, ks) {
  AP <- -2.84
  f <- (1/365)
  AP * (1 - (1 - F0) * exp(-Q * Sd)) * sin(2 * pi * f * t - phi) + ks
}

# eq. 6 in Gallart et al. (2020) - F0 and Sd are estimated from eq. 8
Fyw_model <- function(Q, F0, Sd){
  F0 <- abs(F0)
  Sd <- abs(Sd)
  1 - ((1 - F0) * exp(-Q * Sd))
}

min_q <- 0.001
end <- as.Date("2021-12-31") # data available at the time of manuscript publication

area <- 1306 # drainage area at NEON sampling site

# NEON isotopes
NEON_isotopes <- file.path("data", "NEON", "NEON_isotopes.csv") %>% 
  read_csv(show_col_types = F)

# NEON stream δ18O at outlet
NEON_stream_isotopes <- NEON_isotopes %>%
  filter(type == "SW", date <= end) %>% 
  mutate(
    t = as.numeric(date - date[1]),
  )

# Streamflow at USGS 06879650
flow <- file.path("data", "USGS", "USGS_06879650.csv") %>% 
  read_csv(show_col_types = F) %>%
  filter(date >= min(NEON_stream_isotopes$date) & date <= end) %>%
  mutate(
    flow_mmd = ((flow_m3s*1e9) / (area*1e10))*86400 # m3/s to mm/d
  )

# streamflow weights
NEON_stream_isotopes <- NEON_stream_isotopes %>%
  left_join(select(flow, date, flow_m3s, flow_mmd), by="date")

NEON_stream_isotopes <- na.omit(NEON_stream_isotopes)



# Young water fraction sensitivity (2015-2022)
t <- NEON_stream_isotopes$t
y <- NEON_stream_isotopes$d18O
Q <- NEON_stream_isotopes$flow_mmd
Q_range <- seq(min(Q), max(Q), length.out = 1e4)

p0 <- c(F0 = 0.1, Sd = 0.1, phi = 0.0, ks = 0.0) # initial parameter guesses

# fit isotope model
fit <- nls(y ~ isotope_model(Q, t, F0, Sd, phi, ks), data = list(Q=Q, t=t, y=y), start = p0)

# fitted parameters
coefs  <- summary(fit)$coefficients
F0     <- coefs["F0", "Estimate"]
F0_se  <- coefs["F0", "Std. Error"]
Sd     <- coefs["Sd", "Estimate"]
Sd_se  <- coefs["Sd", "Std. Error"]

# parameter covariance matrix
vc <- vcov(fit)

# Perform Gaussian error propagation
n = 1e3

out <- vector("list", n)
for (i in 1:n){

  params <- MASS::mvrnorm(1, coef(fit), vc)
  
  out[[i]] <- tibble(
    n = i,
    flow_mmd = Q_range,
    F0 = params["F0"],
    Sd = params["Sd"],
    Fyw = Fyw_model(Q = Q_range, F0 = params["F0"], Sd = params["Sd"])
  )
  
}
sens_fit <- tibble(flow_mmd = Q_range, Fyw = Fyw_model(Q = Q_range, F0 = F0, Sd = Sd))
sens_fits <- bind_rows(out)



# Young water fraction sensitivity (2021 water year) 
NEON_stream_isotopes_2021 <- filter(NEON_stream_isotopes, date >= as.Date("2020-10-01") & date <= as.Date("2021-09-30"))

t <- NEON_stream_isotopes_2021$t
y <- NEON_stream_isotopes_2021$d18O
Q <- NEON_stream_isotopes_2021$flow_mmd

p0 <- c(F0 = 0.1, Sd = 0.1, phi = 0.0, ks = 0.0) # initial parameter guesses

# fit isotope model
fit_2021 <- nls(y ~ isotope_model(Q, t, F0, Sd, phi, ks), data = list(Q=Q, t=t, y=y), start = p0)

# fitted parameters
coefs_2021  <- summary(fit_2021)$coefficients
F0_2021     <- coefs_2021["F0", "Estimate"]
F0_se_2021  <- coefs_2021["F0", "Std. Error"]
Sd_2021     <- coefs_2021["Sd", "Estimate"]
Sd_se_2021  <- coefs_2021["Sd", "Std. Error"]

# parameter covariance matrix
vc <- vcov(fit_2021)

out <- vector("list", n)
for (i in 1:n){
  
  params <- MASS::mvrnorm(1, coef(fit_2021), vc)
  
  out[[i]] <- tibble(
    n = i,
    flow_mmd = Q_range,
    F0 = params["F0"],
    Sd = params["Sd"],
    Fyw = Fyw_model(Q = Q_range, F0 = params["F0"], Sd = params["Sd"])
  )
  
}
sens_fit_2021 <- tibble(flow_mmd = Q_range, Fyw = Fyw_model(Q = Q_range, F0 = F0_2021, Sd = Sd_2021))
sens_fits_2021 <- bind_rows(out)

F0_points <- tibble(
  period = c("2015-2022", "2021"),
  flow_mmd = c(0,0),
  F0 = c(abs(F0), abs(F0_2021)),
)

parameters <- paste0("(1) F0: ", round(abs(F0)*100, 1), "% Sd: ", round(abs(Sd), 2), " d/mm")
print(parameters)
parameters_2021 <- paste0("(2) F0: ", round(abs(F0_2021)*100, 1), "% Sd: ", round(abs(Sd_2021), 2), " d/mm")
print(parameters_2021)

p1 = ggplot()+
  #geom_line(data = sens_fits, aes(flow_mmd, Fyw*100, group = n), color = "#999999", alpha = 0.04)+
  #geom_line(data = sens_fits_2021, aes(flow_mmd, Fyw*100, group = n), color = "#E6194B", alpha = 0.04)+
  geom_line(data = sens_fit, aes(flow_mmd, Fyw*100), color= "#373737", linetype="dashed", lwd=0.6)+
  geom_line(data = sens_fit_2021, aes(flow_mmd, Fyw*100), color = "#E6194B", linetype="dashed", lwd=0.6)+
  geom_point(data = F0_points, aes(flow_mmd, F0*100, fill = period), color = "#373737", pch=21, size=2.5)+
  # annotated in inkscape
  annotate("text", x=60, y=80, label = parameters,  color = "#373737", size=10/.pt)+
  annotate("text", x=60, y=74, label = parameters_2021, color = "#E6194B", size=10/.pt)+
  scale_fill_manual(
    name = "",
    values = c(
      "2015-2022" = "#999999",
      "2021" = "#E6194B"
    ),
    guide = "none"
  )+
  labs(x = "Streamflow (mm/d)", y = expression(F["yw"]* " (%)"))+
  theme(
    axis.title.y = element_text(size=12)
  )

p1

#ggsave(file.path("figures", "Fig6.png"), plot=p1, dpi=300, width=110, height=90, units="mm")
ggsave(file.path("figures", "Fig6_annotated.png"), plot=p1, dpi=300, width=110, height=90, units="mm")

print(paste0("F0 (2015-2022): ", round(abs(F0)*100, 1), "% ± ", round(abs(F0_se)*100, 1)))
print(paste0("Sd (2015-2022): ", round(abs(Sd), 2), "d/mm ± ", round(abs(Sd_se), 2)))

print(paste0("F0 (2021): ", round(abs(F0_2021)*100, 1), "% ± ", round(abs(F0_se_2021)*100, 1)))
print(paste0("Sd (2021): ", round(abs(Sd_2021), 2), "d/mm ± ", round(abs(Sd_se_2021), 2)))
