# Fig S5. δ18O & δ2H isotope biplot
library(sf)
library(tidyverse)

source(file.path("code", "theme.R"))

end <- as.Date("2021-12-31") # data available at the time of manuscript publication

AIMS_isotopes <- file.path("data", "AIMS", "AIMS_isotopes.csv") %>% 
  read_csv(show_col_types = F) %>% 
  mutate(
    month = factor(month, levels = 6:8, labels = c("June", "July", "August"))
  )

NEON_precip_isotopes <- file.path("data", "NEON", "NEON_isotopes.csv") %>% 
  read_csv(show_col_types = F) %>%
  filter(type == "PPT") %>% 
  filter(date <= end) 

# meteoric water line (MWL) slope & intercept
MWL <- lm(d2H ~ d18O, data = NEON_precip_isotopes)
MWL_intercept <- MWL$coefficients[1]
MWL_slope <- MWL$coefficients[2]

# evaporation Line (EL) slope & intercept
EL <- lm(d2H ~ d18O, data = AIMS_isotopes)
EL_intercept <- EL$coefficients[1]
EL_slope <- EL$coefficients[2]

ggplot()+
  geom_abline(intercept = MWL_intercept, slope = MWL_slope)+
  geom_abline(intercept = EL_intercept, slope = EL_slope)+
  annotate("text", x = -3.2, y = -9.0, label = "LMWL", size = 9/.pt)+
  annotate("text", x = -3, y = -17.5, label = "EL", size = 9/.pt)+
  geom_point(data = AIMS_isotopes, aes(d18O, d2H, fill = month), color="373737", pch=21, size=2.5)+
  scale_fill_manual(
    name = "",
    values = c(
      "June" = "#0082C8",
      "July" = "#3CB44B",
      "August" = "#E6194B"
    )
  )+
  labs(x = expression(delta^{18}*"O (‰)"), y = expression(delta^{2}*"H (‰)"))
ggsave(file.path("figures", "FigS5.png"), dpi=300, width=90, height=90, units="mm")
