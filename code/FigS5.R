#Fig S5. Random forest model to predict δ18O
library(viridis)
library(party)
library(patchwork)
library(tidyverse)

#Publication theme
source("code/Theme+Settings.R")

#Read: AIMS isotopes RF
df <- as_tibble(read.csv("data/AIMS_isotopes_RF.csv"))
df$watershed <- as.factor(df$watershed)
df$flowing <- as.factor(df$flowing)
df$ysiWaterTemp <- as.numeric(df$ysiWaterTemp)
names(df)

#Define predictors
predictors <- c("watershed", "yday", "Elevation_m", "ContributingArea_ha", "Slope_prc", "TWI", "flowing", "ysiWaterTemp")

#Make data frame with some descriptive variables, the target, and the predictor variables
fit_data_in <-
  df %>%
  select(c("sampleID", "siteID", "long", "lat", "month", "d18OWater", all_of(predictors))) %>%
  subset(complete.cases(.)) #remove any row with an NA value

#Fit model
rfmodel <- cforest(d18OWater ~ .,
                   data = select(fit_data_in, -sampleID, -siteID, -long, -lat, -month), 
                   controls = cforest_unbiased(ntree = 50, mtry = 3))

#Extract more robust variable importance 
rfmodel_varimp <- varimp(rfmodel, conditional = T)
rfmodel_varimp <- sort(rfmodel_varimp, decreasing = T)
var_imp <- tibble(predictor = names(rfmodel_varimp), 
                  IncMSE = rfmodel_varimp)
var_imp$predictor[c(1,2,3,4,5,6,7,8)] <- c("Day of Year", "Flow State", "Water Temperature", "Topographic Wetness Index", "Contributing Area", "Burn Frequency", "Elevation", "Slope")
var_imp <- arrange(var_imp, -IncMSE) %>% #order most to least importance (higher value = greater influence on predictors)
  mutate(predictor = factor(predictor, levels = .$predictor)) #make a factor so they plot in order

#Plot variable importance!
p1 <-
  ggplot(var_imp, aes(x = predictor, y = IncMSE))+
  geom_col(fill = "grey40")+
  labs(x = "Predictor Factor", y = "Variable Importance")+
  theme(axis.text = element_text(angle = 45, hjust = 1))

#Evaluate RF (model diagnostics)
rfmodel_pred <- unlist(treeresponse(rfmodel)) #[c(FALSE,TRUE)]
fit_data_in$predicted <- predict(rfmodel, newdata = fit_data_in, type = "response")

#Plot model fit!
fit_data_in$synopticDate <- factor(fit_data_in$month, levels = c("Jun", "Jul", "Aug"))
p2 <-
  ggplot(fit_data_in, aes(x = predicted, y = d18OWater, fill = synopticDate))+
  geom_point(pch = 21, size = 2)+
  geom_abline(intercept = 0, slope = 1, color = "red")+
  scale_x_continuous(name = "Predicted \U03B4\U00B9\U2078O (‰)")+
  scale_y_continuous(name = "Observed \U03B4\U00B9\U2078O (‰)")+
  scale_fill_manual(values = c('#0082c8', '#3cb44b', '#e6194b'), name = "", labels = c("June", "July", "August"))+
  annotate("text", x = -5.75, y = 0.25, label = "RMSE: 0.764", size = 10/.pt)

sm <- summary(lm(fit_data_in$d18OWater ~ fit_data_in$predicted))

#rmse between obs and predicted, in the same units of obs and predicted, provides the standard deviation of the model prediction error. A smaller value indicates better model performance.
rmse <- function(sm) 
  sqrt(mean(sm$residuals^2))
rmse_lm <- rmse(sm)

p1 + p2 + plot_annotation(tag_levels = "a")
ggsave(path = "figures/", "FigS1.png", dpi=300, width = 190, height = 90, units = "mm")