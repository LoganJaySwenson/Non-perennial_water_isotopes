# Fig S6. Random forest model to predict δ18O
library(party)
library(patchwork)
library(tidyverse)

source(file.path("code", "theme.R"))

set.seed(1)

AIMS_sampling_notes <- file.path("data", "AIMS", "AIMS_sampling_notes.csv") %>%
  read_csv(show_col_types = F) %>%
  mutate(
    date = mdy(date),
    month = month(date),
    ysiWaterTemp = as.double(na_if(na_if(str_remove(ysiWaterTemp, "C"), "BLANK"), "0"))
  )

AIMS_isotopes <- file.path("data", "AIMS", "AIMS_isotopes.csv") %>% 
  read_csv(show_col_types = F) %>% 
  left_join(select(AIMS_sampling_notes, siteID, month, ysiWaterTemp), by=c("siteID", "month")) %>% 
  mutate(
    yday = yday(date),
    month = factor(month, levels = 6:8, labels = c("June", "July", "August")),
    flowing = factor(flowing, levels = c(T,F), labels = c("Flowing", "Pooled")),
    watershed = factor(watershed, level = c("SFKC", "N01B", "N02B", "N04D","N20B")),
  ) %>%
  # remove samples with missing water temperature measurements (n = 3)
  drop_na()

AIMS_isotopes <- AIMS_isotopes %>%
  group_by(month) %>%
  mutate(
    set = if_else(
      row_number() %in% sample.int(n(), size = ceiling(0.8 * n())),
      "train",
      "validation"
    )
  ) %>%
  ungroup()

train <- filter(AIMS_isotopes, set == "train")
validation <- filter(AIMS_isotopes, set == "validation")

predictors <- c("Burn frequency" = "watershed",
                "Day of year" = "yday",
                "Elevation" = "elevation_m",
                "Contributing area" = "contributing_area_ha",
                "Slope" = "slope_pct",
                "Topographic wetness index" = "twi",
                "Flow state" = "flowing",
                "Water temperature" = "ysiWaterTemp"
                )

model <- cforest(d18O ~ .,
                 data = select(train, c("d18O", unname(predictors))),
                 controls = cforest_unbiased(ntree = 50, mtry = 3)
                 )

# conditional permutation importance (more-robust)
importance <- varimp(model, conditional=T)

importance <- tibble(
  predictor = names(importance),
  importance = as.numeric(importance)
  ) %>%
  arrange(desc(importance)) %>%
  mutate(
    predictor = factor(predictor, levels = predictor, labels = names(predictors)[match(predictor, predictors)])
  )

train$d18O_predicted <- predict(model, type="response", newdata=select(train, c("d18O", unname(predictors))))[,1]
validation$d18O_predicted <- predict(model, type="response", newdata=select(validation, c("d18O", unname(predictors))))[,1]

rmse <- sqrt(mean((validation$d18O - validation$d18O_predicted)**2, na.rm = T))

fit <- bind_rows(train, validation)

p1 = ggplot()+
  geom_col(data = importance, aes(predictor, importance), fill="#666666", color="#373737")+
  labs(x = "", y = "Variable importance")+
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    plot.margin = margin(l = 3, unit="mm")
  )

p2 = ggplot()+
  geom_point(data = fit, aes(d18O_predicted, d18O, fill=month), color="#373737", pch=21, size=2.5)+
  geom_abline(intercept = 0, slope = 1, color = "#E7315D")+
  annotate("text", x = -Inf, y = Inf, label = paste0("RMSE: ", round(rmse, 2), "‰"), hjust=-0.02, vjust=1.2, size=10/.pt)+
  scale_fill_manual(
    name = "",
    values = c(
      "June" = "#0082C8",
      "July" = "#3CB44B",
      "August" = "#E6194B"
    )
  )+
  labs(x = expression(Predicted ~ delta^{18}*"O (‰)"), y = expression(Observed ~ delta^{18}*"O (‰)"))

p1 + p2
ggsave(file.path("figures", "FigS6.png"), dpi=300, width=190, height=90, units="mm")
