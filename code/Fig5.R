# Fig 5. Variation in δ18O with d-excess & flow state 
library(patchwork)
library(tidyverse)

source(file.path("code", "theme.R"))

AIMS_isotopes <- file.path("data", "AIMS", "AIMS_isotopes.csv") %>% 
  read_csv(show_col_types = F) %>% 
  mutate(
    month = factor(month, levels = 6:8, labels = c("June", "July", "August")),
    flowing = factor(flowing, levels = c(T,F), labels = c("Flowing", "Pooled"))
  )

p1 = ggplot()+
  geom_point(data = AIMS_isotopes, aes(d18O, dexcess, fill = month), color = "#373737", pch=21, size=2.5)+
  scale_fill_manual(
    name = "",
    values = c(
      "June" = "#0082C8",
      "July" = "#3CB44B",
      "August" = "#E6194B"
    )
  )+
  labs(x = expression(delta^{18}*"O (‰)"), y = "d-excess (‰)")

# Unpaired Wilcoxon rank sum test to see if the differences in medians of isotopic signatures between flowing vs. pooled reaches is equal to 0.
wilcox <- AIMS_isotopes %>%
  group_by(month) %>%
  summarise(
    num_flowing = sum(flowing == "Flowing"),
    num_pooled  = sum(flowing == "Pooled"),
    median_flowing = median(d18O[flowing == "Flowing"], na.rm=T),
    median_pooled  = median(d18O[flowing == "Pooled"],  na.rm=T),
    p_value = if (num_flowing > 0 & num_pooled > 0) {
      wilcox.test(
        d18O[flowing == "Flowing"],
        d18O[flowing == "Pooled"],
        exact=F
      )$p.value
    } else {
      NA
    },
    .groups = "drop"
  )

labels <- AIMS_isotopes %>%
  group_by(month, flowing) %>%
  summarise(
    n = n(),
    y = max(d18O, na.rm = T) + 0.2,
    .groups = "drop"
  ) %>%
  mutate(label = paste("n =", n))

p2 = ggplot()+
  geom_boxplot(data = AIMS_isotopes, aes(month, d18O, fill = flowing), width=0.5, outlier.shape=21, outlier.size=2.5, position = position_dodge(width = 0.7))+
  geom_text(data = labels, aes(month, y, label = label, group = flowing), position = position_dodge(width = 0.7), size = 8/.pt, vjust = 0)+
  labs(x = "", y = expression(delta^{18}*"O (‰)"))+
  scale_fill_manual(
    name = "",
    values = c(
      "Flowing" = "#619CFF",
      "Pooled" = "#F8766D"
    )
  )

pp = p1 + p2
pp
ggsave(file.path("figures", "Fig5.png"), dpi=300, width=190, height=90, units="mm")
