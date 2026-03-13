# Fig 9. Bayesian-estimated young water fractions in headwaters with different priors 
library(ggtext)
library(patchwork)
library(tidyverse)

source(file.path("code", "theme.R"))

crs <- "epsg:5070"

spatial_priors <- file.path("data", "spatial_priors.csv") %>% 
  read_csv(show_col_types = F) %>%
  mutate(
    month = factor(month, levels = 6:8, labels = c("June", "July", "August")),
    prior = factor(prior, 
                   levels = c("c(0.5, 0.5)", "c(0.145, 0.855)", "c(0.038, 0.962)")),
    s1 = s1*100,
    s2 = s2*100
  )

priors <- levels(spatial_priors$prior)

titles <- list(
  expression("Uninformed (0.5/0.5)"),
  expression("Flow-weighted" ~ F["yw"] ~ "(0.15/0.85)"),
  expression("Time-weighted" ~ F["yw"] ~ "(0.04/0.96)")
)

format_p <- function(pval, sig) {
  pval_str <- format(pval, scientific=T, digits=2)
  if_else(
    sig,
    paste0("<b>", pval_str, "*</b>"),
    paste0(pval_str)
  )
}

# Pairwise t-test to test if differences in means of young water between months per prior are equal to 0.
pairwise_t <- spatial_priors %>%
  group_by(prior) %>%
  summarise(
    pt = list(pairwise.t.test(
      s1,
      month,
      p.adjust.method = "bonferroni"
    )),
    .groups = "drop"
  ) %>%
  mutate(
    t1_pval = map_dbl(pt, ~ .x$p.value["July", "June"]),
    t2_pval = map_dbl(pt, ~ .x$p.value["August", "July"]),
    t3_pval = map_dbl(pt, ~ .x$p.value["August", "June"]),
    t1_significance = t1_pval < 0.05,
    t2_significance = t2_pval < 0.05,
    t3_significance = t3_pval < 0.05,
    label = paste0(
      "June to July: ", format_p(t1_pval, t1_significance), "<br>",
      "July to August: ", format_p(t2_pval, t2_significance), "<br>",
      "June to August: ", format_p(t3_pval, t3_significance)
    )
  ) 

pp <- list()

for (i in seq_along(priors)) {
  
  spatial_prior <- filter(spatial_priors, prior == priors[i])
  t <- filter(pairwise_t, prior == priors[i])
  label <- t$label[1]
  
  p = ggplot() +
    geom_boxplot(data = spatial_prior, aes(x = month, y = s1), fill = "#BBBBBB", width = 0.5, outlier.shape = 21, outlier.size = 2) +
    stat_summary(data = spatial_prior, aes(x = month, y = s1), 
                 fun = mean, geom = "point", shape = 23, size = 3, fill = "#AA3377") +
    annotate("richtext", x = 0.4, y = -Inf, label = label,
             hjust = 0, vjust = 0, size = 9/.pt,
             fill = "transparent", label.color = NA,
             color = "#373737", label.padding = unit(0.5, "lines")) +
    labs(
      x = "",
      y = "",
      title = titles[[i]]
    ) +
    scale_x_discrete(labels = c("June", "July", "August"))+
    scale_y_continuous(breaks = scales::pretty_breaks(n=4))+
    theme(
      axis.title = element_text(size=12)
    )
  
  if (i == 1)
    p <- p + labs(y = expression(F["yw"]* " (%)"))
  
  pp[[i]] <- p
}

pp <- wrap_plots(pp, ncol=3)
pp
ggsave(file.path("figures", "Fig9.png"), dpi=300, width=210, height=90, units="mm")