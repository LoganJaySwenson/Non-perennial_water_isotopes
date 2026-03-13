# Fig 7. Bayesian-estimated young water fractions at catchment outlet with different priors
library(tidyverse)

source(file.path("code", "theme.R"))

temporal_priors <- file.path("data", "temporal_priors.csv") %>% 
  read_csv(show_col_types = F) %>%
  mutate(
    prior = factor(prior, 
                   levels = c("c(0.5, 0.5)", "c(0.145, 0.855)", "c(0.038, 0.962)")),
    s1 = s1*100,
    s2 = s2*100
  )  

priors <- levels(temporal_priors$prior)

temporal_summary <- temporal_priors %>%
  group_by(prior) %>%
  summarise(
    mean_s1 = mean(s1, na.rm=T),
    sd_s1 = sd(s1, na.rm=T),
    min_s1 = min(s1, na.rm=T),
    max_s1 = max(s1, na.rm=T),
    label = paste0(sprintf("%.1f", mean_s1), "% ± ", sprintf("%.1f", sd_s1)),
    .groups = "drop"
  ) 

titles <- c(
  "c(0.5, 0.5)"       = "'Uninformed'~'(0.5/0.5)'",
  "c(0.145, 0.855)"   = "'Flow-weighted'~F[yw]~'(0.15/0.85)'",
  "c(0.038, 0.962)" = "'Time-weighted'~F[yw]~'(0.04/0.96)'"
)

ggplot()+
  geom_point(data = temporal_priors, aes(date, s1, fill = prior), color = "#373737", pch=21, size=2.5)+
  geom_hline(data = temporal_summary, aes(yintercept = mean_s1, color = prior), linetype = "dashed", lwd=0.8)+
  geom_text(data = temporal_summary, aes(x = as.Date("2020-10-01"), y = mean_s1, label = label), hjust=0, vjust=-1.05, size = 10/.pt)+
  labs(x = "", y = expression(F["yw"]* " (%)"), fill = "Prior")+
  scale_fill_manual(
    values = c(
      "c(0.5, 0.5)" = "#9ecae1",
      "c(0.145, 0.855)" = "#4292c6",
      "c(0.038, 0.962)" = "#084594"
    ),
    labels = c(
      "c(0.5, 0.5)" = "0.5/0.5",
      "c(0.145, 0.855)" = "0.15/0.85",
      "c(0.038, 0.962)" = "0.04/0.96"
    )
  ) +
  scale_color_manual(
    values = c(
      "c(0.5, 0.5)" = "#9ecae1",
      "c(0.145, 0.855)" = "#4292c6",
      "c(0.038, 0.962)" = "#084594"
    ),
    guide = "none" 
  ) +
  scale_x_date(
    breaks = scales::date_breaks("months"), 
    labels = scales::date_format("%b %Y")
    )+
  scale_y_continuous(
    breaks = scales::pretty_breaks(n=4)
  )+
  facet_wrap(~prior, ncol = 1, nrow = 3, scales = "free_y", labeller = as_labeller(titles, label_parsed))+
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, hjust = 0),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
  )
ggsave(file.path("figures", "Fig7.png"), dpi=300, width=190, height=190, units="mm")
