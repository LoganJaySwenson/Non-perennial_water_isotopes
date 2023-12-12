# Fig 5. Variation in δ18O with d-excess & flow state
library(patchwork)
library(lubridate)
library(tidyverse)

# Publication theme
source("code/Theme+Settings.R")

# Read: AIMS isotopes
AIMS_isotopes <- read_csv("data/AIMS_isotopes_RF.csv") %>%
  select(1:18)

# Plot δ18O & d-excess
AIMS_isotopes$month <- factor(AIMS_isotopes$month, levels = c( "Jun", "Jul", "Aug"))
p1 <-
  ggplot()+
  geom_point(data = AIMS_isotopes, aes(d18OWater, dexcess, fill = month), color= 'black', pch=21, size=2)+
  scale_fill_manual(values = c('#0082c8', '#3cb44b', '#e6194b'), labels = c("June", "July", "August"))+
  labs(x = expression(delta^{18}*"O (‰)"), y = "d-excess (‰)")+
  labs(fill = "", shape = "")

# Unpaired Wilcoxon rank sum test to see if the differences in medians of isotopic signatures between flowing vs. pooled reaches is equal to 0.
# p-value < 0.05 indicates a significant difference.
wilcox.test(AIMS_isotopes %>% filter(month == "Jul" & flowing == "y") %>% pull(d18OWater),
            AIMS_isotopes %>% filter(month == "Jul" & flowing == "n") %>% pull(d18OWater), exact = F)

wilcox.test(AIMS_isotopes %>% filter(month == "Aug" & flowing == "y") %>% pull(d18OWater),
            AIMS_isotopes %>% filter(month == "Aug" & flowing == "n") %>% pull(d18OWater), exact = F)

# Plot δ18O & flow state
AIMS_isotopes$flowing <- factor(AIMS_isotopes$flowing, levels = c("y", "n"))
p2 <-
  ggplot(AIMS_isotopes, aes(month, d18OWater, fill = flowing))+
  geom_boxplot(width = 0.5, outlier.shape = 21, outlier.size = 2)+
  annotate("text", x = 1.04, y = -5.50, label = paste("n =", count(subset(AIMS_isotopes, month == "Jun"))), size = 6/.pt)+
  annotate("text", x = 1.88, y = -5.20, label = paste("n =", count(subset(AIMS_isotopes, month == "Jul" & flowing == "y"))), size = 6/.pt)+
  annotate("text", x = 2.16, y = -4.80, label = paste("n =", count(subset(AIMS_isotopes, month == "Jul" & flowing == "n"))), size = 6/.pt)+
  annotate("text", x = 2.88, y = -5.30, label = paste("n =", count(subset(AIMS_isotopes, month == "Aug" & flowing == "y"))), size = 6/.pt)+
  annotate("text", x = 3.16, y = 0.5, label = paste("n =", count(subset(AIMS_isotopes, month == "Aug" & flowing == "n"))), size = 6/.pt)+
  labs(x = "", y = expression(delta^{18}*"O (‰)"))+
  scale_fill_manual(values = c("#619CFF", "#F8766D"), name = "", labels = c("Flowing", "Pooled"))+
  scale_x_discrete(labels = c("June", "July", "August"))

p1 + p2 + plot_annotation(tag_levels = "a")
ggsave(path = "figures/", "Fig5.png", dpi=300, width = 190, height = 90, units = "mm")
ggsave(path = "figures/", "Figure5.pdf", dpi=300, width = 190, height = 90, units = "mm")